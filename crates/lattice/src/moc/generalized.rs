//! 一般化随伴の固定源の問題（ADR-0247 の段 3。#0279）。
//!
//! 転置した物質の固定源の問題 (L̃ − F̃/k) Γ = S を、前向きの固有値計算と同じスイープで解く。L̃ − F̃/k は k で特異で、右の
//! 零ベクトルは随伴の固有値計算の解 φ†、左の零ベクトルは前向きの解 φ（トラックを両向きに同じ重みで通るスイープの相反性。
//! `solve_adjoint`）。解は ⟨S, φ⟩ = 0 のときだけあり、φ† の倍だけ不定なので、⟨Γ, Fφ⟩ = 0 で 1 つに決める。反復の状態（代表の
//! 領域のモーメントとトラックの端の角度中性子束）から、スイープごとに φ† の状態の倍を除く。除かなければ、源の丸めで入る φ† の
//! 成分が固有値 1 のモードとして残る。
//!
//! 源の確かめ ⟨S, φ⟩ と直交化 ⟨Γ, Fφ⟩ の内積は、領域のスカラー中性子束の体積の重みの和で、前向きと随伴の角度中性子束の
//! 積（`Importance` の全断面積の項）を使わない。どちらも領域ごとの等方な掛け算（源、核分裂）だけを含み、平らな線源から領域の
//! 中性子束への作用素の摂動を含まないので、領域の値の内積が離散の問題について厳密になる（ADR-0247 の「一般化随伴」）

use std::f64::consts::PI;

use rayon::prelude::*;

use super::{
    CoarseMesh, Detected, Material, MocError, Options, Remaining, Representatives, Start, Sweeper,
    Symmetry, equilibrium_moments, planar_moments, sum_tallies, transposed,
};
use crate::cmfd::{CmfdOptions, Coarse};
use crate::tracking::Tracking;

/// ⟨S, φ⟩ を ⟨|S|, φ⟩ と比べる相対の許容値。応答の源 h₁/⟨h₁, φ⟩ − h₂/⟨h₂, φ⟩ を解いた φ から作れば、残るのは足す順の丸め
/// だけなので、それより十分大きく、反復の判定（既定の中性子束 1e-6）より十分小さい値にする
const ORTHOGONALITY_TOLERANCE: f64 = 1.0e-10;

/// 固定源の反復の収束の方法
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GeneralizedMethod {
    /// 射影つきの源の反復（加速なし）
    SourceIteration,
    /// 射影つきの源の反復を、スイープごとに粗メッシュの固定源の問題の解で補正する
    Cmfd(CmfdOptions),
    /// 射影した作用素の GMRES（`restart` 回ごとに再始動する）
    Gmres { restart: usize },
}

/// 一般化随伴の固定源の反復の方法と打ち切り
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeneralizedOptions {
    pub method: GeneralizedMethod,
    /// 中性子束の変化から見積もった残りの誤差の上限（群ごとの、領域の |Γ| の最大値に対する比）
    pub tolerance: f64,
    /// スイープの回数の上限
    pub max_iterations: usize,
    /// すべての源を一緒にスイープする（減衰 1 − exp(−τ) を源の間で共有する）。`false` なら源ごとに解く
    pub batch: bool,
}

/// 1 つの源の一般化随伴
#[derive(Debug, Clone, PartialEq)]
pub struct Generalized {
    /// [領域][群] のスカラー中性子束（源 S と同じ規格化。⟨Γ, Fφ⟩ = 0）
    pub flux: Vec<Vec<f64>>,
    /// 判定を満たすまでのスイープの回数（一緒にスイープすれば、すべての源が満たすまでスイープを続ける）
    pub iterations: usize,
}

/// 物質 `materials`（前向き）の問題の一般化随伴を、源 `sources`（[源][領域][群] の等方な放出率の密度）ごとに解く（ADR-0247）。
/// `forward` と `adjoint` は、同じトラッキング、対称操作 `symmetry`、`options` の散乱の次数で解いた前向きと随伴の固有値計算の
/// 収束した状態（`solve_symmetric_from` と `solve_adjoint_from`）で、k は前向きの値を使う。源の値が欠ける、有限でない、対称操作で
/// 移り合う領域で違う、⟨S, φ⟩ が 0 でない（解がない）、状態が別の問題のものなら、エラー
#[allow(clippy::too_many_arguments)]
pub fn solve_generalized(
    tracking: &Tracking,
    materials: &[&Material],
    options: &Options,
    symmetry: Option<&Detected>,
    forward: &Start,
    adjoint: &Start,
    sources: &[Vec<Vec<f64>>],
    settings: &GeneralizedOptions,
) -> Result<Vec<Generalized>, MocError> {
    let tolerance = settings.tolerance;
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        return Err(MocError::GeneralizedTolerance(tolerance));
    }
    let regions = materials.len();
    let ng = materials.first().ok_or(MocError::NoRegions)?.groups();
    if materials.iter().any(|m| m.groups() != ng) {
        return Err(MocError::MixedGroups);
    }
    let (owned, index) = transposed(materials)?;
    let adjoint_materials: Vec<&Material> = index.iter().map(|&i| &owned[i]).collect();
    let mesh = CoarseMesh::new(tracking);
    let pairs = planar_moments(options.scattering_order);
    let representatives = match symmetry {
        Some(verified) if verified.symmetry() != Symmetry::NONE => {
            if !verified.fits(tracking) {
                return Err(MocError::ForeignSymmetry);
            }
            Representatives::by_images(
                tracking,
                &adjoint_materials,
                &mesh,
                &pairs,
                verified.images(),
            )?
        }
        _ => Representatives::all(tracking, &mesh, pairs.len()),
    };
    let copies = if settings.batch { sources.len() } else { 1 };
    let batched = (copies > 1).then(|| {
        Sweeper::with_copies(
            tracking,
            &adjoint_materials,
            ng,
            copies,
            &mesh,
            &pairs,
            representatives.clone(),
        )
    });
    let sweeper = Sweeper::new(
        tracking,
        &adjoint_materials,
        ng,
        &mesh,
        &pairs,
        representatives,
    );
    let sizes = [regions, tracking.tracks.len()];
    if !forward.fits(sizes, &sweeper, ng) || !adjoint.fits(sizes, &sweeper, ng) {
        return Err(MocError::ForeignStart);
    }
    let reduced = &sweeper.representatives;
    let kept = &reduced.regions.representatives;
    let volumes = tracking.volumes(regions);
    let orbit_volumes = reduced.regions.sums(&volumes);
    // [代表][群] の Fφ（前向きの物質の核分裂の源 χ (νΣf · φ)）に軌道の体積を掛けたもの。⟨Γ, Fφ⟩ = Σ weights · Γ
    let weights: Vec<Vec<f64>> = kept
        .iter()
        .zip(&forward.flux)
        .zip(&orbit_volumes)
        .map(|((&r, phi), v)| {
            let mut w = vec![0.0; ng];
            for t in &materials[r].fission_terms {
                let rate: f64 = t.nu_fission.iter().zip(phi).map(|(n, f)| n * f).sum();
                for (x, chi) in w.iter_mut().zip(&t.chi) {
                    *x += v * chi * rate;
                }
            }
            w
        })
        .collect();
    let adjoint_production = dot(&weights, &adjoint.flux);
    if !(adjoint_production > 0.0 && adjoint_production.is_finite()) {
        return Err(MocError::NoFission);
    }
    // [源][代表][群] の源（確かめたもの）
    let mut reduced_sources = Vec::with_capacity(sources.len());
    for (s, source) in sources.iter().enumerate() {
        if source.len() != regions || source.iter().any(|q| q.len() != ng) {
            return Err(MocError::SourceShape { index: s });
        }
        for (r, q) in source.iter().enumerate() {
            if q.iter().any(|v| !v.is_finite()) {
                return Err(MocError::InvalidSource {
                    index: s,
                    region: r,
                });
            }
            if *q != source[kept[reduced.regions.index[r]]] {
                return Err(MocError::AsymmetricSource {
                    index: s,
                    region: r,
                });
            }
        }
        let q: Vec<Vec<f64>> = kept.iter().map(|&r| source[r].clone()).collect();
        // ⟨S, φ⟩ と、その項の絶対値の和（φ は正）
        let (mut overlap, mut size) = (0.0, 0.0);
        for ((q, phi), v) in q.iter().zip(&forward.flux).zip(&orbit_volumes) {
            for (a, b) in q.iter().zip(phi) {
                overlap += v * a * b;
                size += v * a.abs() * b;
            }
        }
        if overlap.abs() > ORTHOGONALITY_TOLERANCE * size {
            return Err(MocError::NonOrthogonalSource {
                index: s,
                ratio: overlap / size,
            });
        }
        reduced_sources.push(q);
    }
    let problem = Problem {
        single: &sweeper,
        batched: batched.as_ref(),
        mesh: &mesh,
        materials: &adjoint_materials,
        volumes: &volumes,
        pairs: &pairs,
        ng,
        k: forward.k,
        weights: &weights,
        adjoint,
        adjoint_production,
    };
    let coarse = match settings.method {
        GeneralizedMethod::Cmfd(cmfd) => Some(problem.coarse(cmfd)?),
        _ => None,
    };
    let solve = |batch: &[Vec<Vec<f64>>]| -> Result<Vec<(State, usize)>, MocError> {
        match settings.method {
            GeneralizedMethod::Gmres { restart } => problem.gmres(batch, restart, settings),
            _ => problem.source_iteration(batch, coarse.as_ref(), settings),
        }
    };
    let solved = if settings.batch {
        solve(&reduced_sources)?
    } else {
        let mut solved = Vec::with_capacity(reduced_sources.len());
        for q in &reduced_sources {
            solved.extend(solve(std::slice::from_ref(q))?);
        }
        solved
    };
    Ok(solved
        .into_iter()
        .map(|(state, iterations)| Generalized {
            flux: reduced.expand(&state.flux),
            iterations,
        })
        .collect())
}

/// Σ_代表 Σ_群 a b
fn dot(a: &[Vec<f64>], b: &[Vec<f64>]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| x.iter().zip(y).map(|(p, q)| p * q).sum::<f64>())
        .sum()
}

/// 反復の状態
#[derive(Clone)]
struct State {
    /// [代表][群] のスカラー中性子束
    flux: Vec<Vec<f64>>,
    /// [代表][モーメント (1 番目から)][群] の高次のモーメント
    higher: Vec<Vec<f64>>,
    /// [代表のトラック][向き][極角][群] の出ていく角度中性子束
    outgoing: Vec<f64>,
}

impl State {
    /// 並べた値（中性子束、高次のモーメント、角度中性子束の順）
    fn pack(&self) -> Vec<f64> {
        self.flux
            .iter()
            .chain(&self.higher)
            .flatten()
            .chain(&self.outgoing)
            .copied()
            .collect()
    }

    /// `pack` の逆。形は `like` と同じ
    fn unpack(values: &[f64], like: &State) -> State {
        let mut rest = values;
        let mut take = |n: usize| {
            let (head, tail) = rest.split_at(n);
            rest = tail;
            head.to_vec()
        };
        let flux = like.flux.iter().map(|f| take(f.len())).collect();
        let higher = like.higher.iter().map(|h| take(h.len())).collect();
        let outgoing = take(like.outgoing.len());
        State {
            flux,
            higher,
            outgoing,
        }
    }
}

/// 粗メッシュの加速の問題（固有値 k で均質化した、転置した物質の問題）
struct CoarseAcceleration {
    coarse: Coarse,
    /// [代表のセル][群] の W（`Coarse::solve_projected`）。粗メッシュの補正 e を延ばした値の ⟨Πe, Fφ⟩ が W · e になる
    weights: Vec<f64>,
    options: CmfdOptions,
}

/// 解く問題（転置した物質のスイープと、φ† を除くための値）
struct Problem<'a> {
    /// 1 つの源のスイープ
    single: &'a Sweeper<'a>,
    /// すべての源を一緒にスイープするスイープ（`GeneralizedOptions::batch` で源が 2 つ以上のとき）
    batched: Option<&'a Sweeper<'a>>,
    mesh: &'a CoarseMesh,
    /// [領域] の転置した物質
    materials: &'a [&'a Material],
    /// [領域] のトラックから推定した体積
    volumes: &'a [f64],
    pairs: &'a [(usize, isize)],
    ng: usize,
    k: f64,
    /// [代表][群] の軌道の体積 × Fφ
    weights: &'a [Vec<f64>],
    /// 随伴の固有値計算の収束した状態（φ† のモーメントと端の角度中性子束）
    adjoint: &'a Start,
    /// ⟨φ†, Fφ⟩
    adjoint_production: f64,
}

/// 源ごとの値 `parts` を `chunk` 個ずつ交互に並べる（一緒にスイープする問題の [..][問題][群] の並び）
fn interleave(parts: &[&[f64]], chunk: usize) -> Vec<f64> {
    let len = parts.first().map_or(0, |p| p.len());
    let mut out = Vec::with_capacity(len * parts.len());
    for at in (0..len).step_by(chunk) {
        for p in parts {
            out.extend_from_slice(&p[at..][..chunk]);
        }
    }
    out
}

/// `interleave` の逆: 問題 `copy` の値
fn strand(values: &[f64], copies: usize, chunk: usize, copy: usize) -> Vec<f64> {
    values
        .chunks(copies * chunk)
        .flat_map(|c| &c[copy * chunk..][..chunk])
        .copied()
        .collect()
}

impl Problem<'_> {
    fn kept(&self) -> &[usize] {
        &self.single.representatives.regions.representatives
    }

    /// `copies` 個の源を一緒にスイープするスイープ
    fn sweeper(&self, copies: usize) -> &Sweeper<'_> {
        match self.batched {
            Some(b) if b.copies == copies => b,
            _ => {
                assert_eq!(copies, 1, "no sweep for {copies} sources");
                self.single
            }
        }
    }

    /// 0 の状態
    fn zero(&self) -> State {
        let states = self.kept().len();
        State {
            flux: vec![vec![0.0; self.ng]; states],
            higher: vec![vec![0.0; (self.pairs.len() - 1) * self.ng]; states],
            outgoing: vec![0.0; self.single.outgoing_len()],
        }
    }

    /// 状態（[代表][群] の中性子束と高次のモーメント）と源 `source`（[代表][群]。`None` なら 0）の、[代表][モーメント][群] の
    /// 角度中性子束の平衡値の展開係数
    fn equilibrium(&self, state: &State, source: Option<&[Vec<f64>]>) -> Vec<Vec<f64>> {
        let kept = self.kept();
        (0..kept.len())
            .into_par_iter()
            .map(|i| {
                let m = self.materials[kept[i]];
                let mut e =
                    equilibrium_moments(m, &state.flux[i], &state.higher[i], self.k, self.pairs);
                // 真空の領域の線源は Σt で割らずに置く（`Sweeper::sweep`）
                if let Some(q) = source {
                    let void = m.is_void();
                    for ((e, q), t) in e.iter_mut().zip(&q[i]).zip(&m.total) {
                        *e += if void {
                            q / (4.0 * PI)
                        } else {
                            q / (4.0 * PI * t)
                        };
                    }
                }
                e
            })
            .collect()
    }

    /// 状態 `states` から、源 `sources`（[問題]。`None` なら 0）を足して一緒に 1 回スイープした状態（φ† の成分を除いたもの）と、
    /// [粗メッシュの面][群] の正味の中性子流
    fn step(
        &self,
        states: &[State],
        sources: &[Option<&[Vec<f64>]>],
        tallies: &mut [super::Tally],
    ) -> Result<Vec<(State, Vec<f64>)>, MocError> {
        let copies = states.len();
        let sweeper = self.sweeper(copies);
        let (ng, pairs) = (self.ng, self.pairs);
        let equilibria: Vec<Vec<Vec<f64>>> = states
            .iter()
            .zip(sources)
            .map(|(state, source)| self.equilibrium(state, *source))
            .collect();
        let equilibrium: Vec<Vec<f64>> = (0..self.kept().len())
            .map(|i| {
                let parts: Vec<&[f64]> = equilibria.iter().map(|e| e[i].as_slice()).collect();
                interleave(&parts, ng)
            })
            .collect();
        let incoming = {
            let parts: Vec<&[f64]> = states.iter().map(|s| s.outgoing.as_slice()).collect();
            interleave(&parts, ng)
        };
        let mut outgoing = vec![0.0; incoming.len()];
        sweeper.sweep(&equilibrium, &incoming, &mut outgoing, tallies)?;
        let mut sums = sum_tallies(sweeper, tallies);
        sweeper
            .representatives
            .symmetrize(&mut sums, &sweeper.void, copies * ng);
        let kept = self.kept();
        let mut next = Vec::with_capacity(copies);
        for (copy, e) in equilibria.iter().enumerate() {
            // [代表][モーメント][問題][群] の和の、この問題の分
            let delta = |i: usize, k: usize| {
                &sums.delta[((i * pairs.len() + k) * copies + copy) * ng..][..ng]
            };
            let (flux, higher) = (0..kept.len())
                .into_par_iter()
                .map(|i| {
                    let region = kept[i];
                    let volume = self.volumes[region];
                    let mut higher = vec![0.0; (pairs.len() - 1) * ng];
                    if let Some(v) = sweeper.void[region] {
                        let flux = sums.void_flux[(v * copies + copy) * ng..][..ng]
                            .iter()
                            .map(|f| f / volume)
                            .collect();
                        return (flux, higher);
                    }
                    let total = &self.materials[region].total;
                    let flux = (0..ng)
                        .map(|g| 4.0 * PI * e[i][g] + delta(i, 0)[g] / (total[g] * volume))
                        .collect();
                    for (k, &(l, _)) in pairs.iter().enumerate().skip(1) {
                        let e = &e[i][k * ng..][..ng];
                        for (g, h) in higher[(k - 1) * ng..][..ng].iter_mut().enumerate() {
                            *h = 4.0 * PI * e[g] / (2 * l + 1) as f64
                                + delta(i, k)[g] / (total[g] * volume);
                        }
                    }
                    (flux, higher)
                })
                .unzip();
            let mut state = State {
                flux,
                higher,
                outgoing: strand(&outgoing, copies, ng, copy),
            };
            self.project(&mut state);
            next.push((state, strand(&sums.current, copies, ng, copy)));
        }
        Ok(next)
    }

    /// 状態から φ† の状態の倍を除き、⟨Γ, Fφ⟩ = 0 にする
    fn project(&self, state: &mut State) {
        let alpha = dot(self.weights, &state.flux) / self.adjoint_production;
        let a = self.adjoint;
        state
            .flux
            .par_iter_mut()
            .zip(&a.flux)
            .chain(state.higher.par_iter_mut().zip(&a.higher))
            .for_each(|(x, y)| x.iter_mut().zip(y).for_each(|(x, y)| *x -= alpha * y));
        state
            .outgoing
            .par_iter_mut()
            .zip(&a.outgoing)
            .for_each(|(x, y)| *x -= alpha * y);
    }

    /// 射影つきの源の反復: 0 の状態から、源と前の状態の散乱と核分裂の平衡値でスイープし、φ† の成分を除くことを繰り返す。
    /// `coarse` があれば、スイープごとに粗メッシュの補正を足す（`correct`）。源 `sources` を一緒にスイープし、すべてが判定を
    /// 満たすまで続ける。[源] の状態と、判定を満たしたスイープの回数を返す
    fn source_iteration(
        &self,
        sources: &[Vec<Vec<f64>>],
        coarse: Option<&CoarseAcceleration>,
        settings: &GeneralizedOptions,
    ) -> Result<Vec<(State, usize)>, MocError> {
        let copies = sources.len();
        let mut states = vec![self.zero(); copies];
        let mut tallies = self.sweeper(copies).tallies(self.materials.len());
        let mut histories: Vec<Remaining> = (0..copies).map(|_| Remaining::default()).collect();
        let mut done: Vec<Option<usize>> = vec![None; copies];
        let given: Vec<Option<&[Vec<f64>]>> = sources.iter().map(|q| Some(q.as_slice())).collect();
        for iteration in 1..=settings.max_iterations {
            let swept = self.step(&states, &given, &mut tallies)?;
            for (copy, ((mut next, _), state)) in swept.into_iter().zip(&mut states).enumerate() {
                if let Some(c) = coarse {
                    self.correct(c, state, &mut next)?;
                }
                let change = largest_change(&next.flux, &state.flux);
                tracing::debug!(iteration, copy, change, "generalized adjoint");
                *state = next;
                if histories[copy].update(change) < settings.tolerance && done[copy].is_none() {
                    done[copy] = Some(iteration);
                }
            }
            if done.iter().all(Option::is_some) {
                return Ok(states.into_iter().zip(done.into_iter().flatten()).collect());
            }
        }
        Err(MocError::GeneralizedNotConverged {
            iterations: settings.max_iterations,
        })
    }

    /// φ† の状態を 1 回スイープした中性子束と面の中性子流（生成率 1 に規格化）で、転置した物質の問題を粗メッシュのセルごとに
    /// 均質化する（固有値計算の CMFD と同じ式。補正の D̂ は φ† の中性子流から作るので、均質化した問題は k で特異で、均質化した
    /// φ† が零ベクトル）
    fn coarse(&self, options: CmfdOptions) -> Result<CoarseAcceleration, MocError> {
        let a = self.adjoint;
        let mut tallies = self.single.tallies(self.materials.len());
        let start = State {
            flux: a.flux.clone(),
            higher: a.higher.clone(),
            outgoing: a.outgoing.clone(),
        };
        // φ† の成分を除かずに、スイープした値そのものを使う
        let equilibrium = self.equilibrium(&start, None);
        let mut outgoing = vec![0.0; a.outgoing.len()];
        self.single
            .sweep(&equilibrium, &a.outgoing, &mut outgoing, &mut tallies)?;
        let mut sums = sum_tallies(self.single, &tallies);
        self.single
            .representatives
            .symmetrize(&mut sums, &self.single.void, self.ng);
        let reduced = &self.single.representatives;
        let orbit_volumes = reduced.regions.sums(self.volumes);
        let kept = self.kept();
        let ng = self.ng;
        let swept: Vec<Vec<f64>> = (0..kept.len())
            .map(|i| {
                let region = kept[i];
                let volume = self.volumes[region];
                match self.single.void[region] {
                    Some(v) => sums.void_flux[v * ng..][..ng]
                        .iter()
                        .map(|f| f / volume)
                        .collect(),
                    None => {
                        let total = &self.materials[region].total;
                        (0..ng)
                            .map(|g| {
                                4.0 * PI * equilibrium[i][g]
                                    + sums.delta[i * self.pairs.len() * ng + g]
                                        / (total[g] * volume)
                            })
                            .collect()
                    }
                }
            })
            .collect();
        let rate: f64 = swept
            .iter()
            .zip(&orbit_volumes)
            .enumerate()
            .map(|(i, (phi, v))| {
                v * phi
                    .iter()
                    .zip(&self.materials[kept[i]].nu_fission)
                    .map(|(f, n)| f * n)
                    .sum::<f64>()
            })
            .sum();
        let flux: Vec<Vec<f64>> = swept
            .iter()
            .map(|f| f.iter().map(|x| x / rate).collect())
            .collect();
        let current: Vec<f64> = sums.current.iter().map(|c| c / rate).collect();
        let cells = &reduced.cells;
        let coarse = Coarse::homogenize(
            self.mesh,
            cells,
            self.materials,
            self.volumes,
            &flux,
            &reduced.regions.index,
            &current,
        )?;
        // W_c = Σ_{代表の領域 i ⊂ c} w_i φ†_i / φ̄†_c（延ばす形は φ† の状態なので、⟨Πe, Fφ⟩ = W · e）
        let average = coarse.homogenized_flux();
        let mut weights = vec![0.0; cells.representatives.len() * ng];
        for (i, &region) in kept.iter().enumerate() {
            let c = cells.index[self.mesh.cell(region)];
            for g in 0..ng {
                weights[c * ng + g] += self.weights[i][g] * a.flux[i][g] / average[c * ng + g];
            }
        }
        Ok(CoarseAcceleration {
            coarse,
            weights,
            options,
        })
    }

    /// 粗メッシュの補正: スイープの前の状態 `before` と後の状態 `after` の中性子束の差 δ の散乱と核分裂の源 M̃δ を、粗メッシュの
    /// セルで体積の積分にした残差で、誤差の問題 (A − F/k) e = R M̃δ を解き（`Coarse::solve_projected`）、e を φ† の状態の形で
    /// 延ばして `after` に足し、φ† の成分を除く。収束した解では δ = 0 なので、補正は解を変えない
    fn correct(
        &self,
        c: &CoarseAcceleration,
        before: &State,
        after: &mut State,
    ) -> Result<(), MocError> {
        let (ng, k) = (self.ng, self.k);
        let reduced = &self.single.representatives;
        let cells = &reduced.cells;
        let kept = self.kept();
        // [代表][群] の M̃δ（P0 の散乱と核分裂）
        let residual: Vec<Vec<f64>> = (0..kept.len())
            .into_par_iter()
            .map(|i| {
                let m = self.materials[kept[i]];
                let delta: Vec<f64> = after.flux[i]
                    .iter()
                    .zip(&before.flux[i])
                    .map(|(a, b)| a - b)
                    .collect();
                let mut r = vec![0.0; ng];
                for (from, row) in m.scatter.iter().enumerate() {
                    for (to, s) in row.iter().enumerate() {
                        r[to] += s * delta[from];
                    }
                }
                for t in &m.fission_terms {
                    let rate: f64 = t.nu_fission.iter().zip(&delta).map(|(n, d)| n * d).sum();
                    for (x, chi) in r.iter_mut().zip(&t.chi) {
                        *x += chi * rate / k;
                    }
                }
                r
            })
            .collect();
        // 代表のセルの領域（`Coarse::homogenize` と同じ集め方）
        let mut rhs = vec![0.0; cells.representatives.len() * ng];
        for (region, volume) in self.volumes.iter().enumerate() {
            if let Some(cell) = cells.position(self.mesh.cell(region)) {
                let r = &residual[reduced.regions.index[region]];
                for (x, y) in rhs[cell * ng..][..ng].iter_mut().zip(r) {
                    *x += volume * y;
                }
            }
        }
        let e = c.coarse.solve_projected(k, &rhs, &c.weights, &c.options)?;
        let average = c.coarse.homogenized_flux();
        // [代表のセル][群] の e / φ̄†
        let ratio: Vec<f64> = e.iter().zip(average).map(|(x, f)| x / f).collect();
        let a = self.adjoint;
        after
            .flux
            .par_iter_mut()
            .zip(after.higher.par_iter_mut())
            .enumerate()
            .for_each(|(i, (f, h))| {
                let cell = cells.index[self.mesh.cell(kept[i])];
                for (j, (x, y)) in f.iter_mut().zip(&a.flux[i]).enumerate() {
                    *x += ratio[cell * ng + j] * y;
                }
                for (j, (x, y)) in h.iter_mut().zip(&a.higher[i]).enumerate() {
                    *x += ratio[cell * ng + j % ng] * y;
                }
            });
        self.single
            .add_outgoing(&mut after.outgoing, &a.outgoing, |cell, g| {
                ratio[cells.index[cell] * ng + g]
            });
        self.project(after);
        Ok(())
    }

    /// 射影した作用素の GMRES: x = P T₀ x + P T(0)（T₀ は源のないスイープ、P は φ† の成分を除く射影）を (I − P T₀) x = P T(0) と
    /// 書き、`restart` 回ごとに再始動する GMRES で解く。内積は並べた状態の値の和で、残差のノルムが右辺のノルムの判定倍に
    /// なったら止める。源 `sources` ごとの Arnoldi の過程を、作用素を掛けるスイープを一緒にして進める。[源] の状態と、判定を
    /// 満たしたスイープの回数を返す
    fn gmres(
        &self,
        sources: &[Vec<Vec<f64>>],
        restart: usize,
        settings: &GeneralizedOptions,
    ) -> Result<Vec<(State, usize)>, MocError> {
        let copies = sources.len();
        let zero = self.zero();
        let mut tallies = self.sweeper(copies).tallies(self.materials.len());
        let given: Vec<Option<&[Vec<f64>]>> = sources.iter().map(|q| Some(q.as_slice())).collect();
        let first = self.step(&vec![zero.clone(); copies], &given, &mut tallies)?;
        let mut krylov: Vec<Krylov> = first
            .into_iter()
            .map(|(b, _)| Krylov::new(b.pack(), settings.tolerance, restart))
            .collect();
        let mut sweeps = 1;
        let none: Vec<Option<&[Vec<f64>]>> = vec![None; copies];
        while krylov.iter().any(|k| !k.done) {
            if sweeps >= settings.max_iterations {
                return Err(MocError::GeneralizedNotConverged { iterations: sweeps });
            }
            // 判定を満たした源は 0 をスイープする（一緒にスイープするので、回数は変わらない）
            let vectors: Vec<State> = krylov
                .iter()
                .map(|k| {
                    k.request()
                        .map_or_else(|| zero.clone(), |v| State::unpack(v, &zero))
                })
                .collect();
            let images = self.step(&vectors, &none, &mut tallies)?;
            sweeps += 1;
            for ((k, v), (image, _)) in krylov.iter_mut().zip(&vectors).zip(images) {
                if k.done {
                    continue;
                }
                // (I − P T₀) v
                let image = image.pack();
                let applied: Vec<f64> = v
                    .pack()
                    .par_iter()
                    .zip(&image)
                    .map(|(a, b)| a - b)
                    .collect();
                k.receive(applied, sweeps);
            }
        }
        Ok(krylov
            .into_iter()
            .map(|k| (State::unpack(&k.x, &zero), k.sweeps))
            .collect())
    }
}

/// 1 つの源の、再始動つきの GMRES の過程（`Problem::gmres`）。作用素を掛ける値を `request` で出し、掛けた値を `receive` で受ける
struct Krylov {
    /// 解の近似
    x: Vec<f64>,
    /// 右辺
    b: Vec<f64>,
    /// 止める残差のノルム（判定 × 右辺のノルム）
    target: f64,
    restart: usize,
    /// 正規直交の基底
    basis: Vec<Vec<f64>>,
    /// Givens の回転で上三角にした Hessenberg 行列の列
    columns: Vec<Vec<f64>>,
    rotations: Vec<(f64, f64)>,
    /// 最小二乗の右辺
    g: Vec<f64>,
    /// 次に掛ける値が解の近似 x（再始動の残差）か
    residual: bool,
    done: bool,
    /// 判定を満たしたスイープの回数
    sweeps: usize,
}

impl Krylov {
    fn new(b: Vec<f64>, tolerance: f64, restart: usize) -> Self {
        let target = tolerance * norm(&b);
        let mut k = Self {
            x: vec![0.0; b.len()],
            b,
            target,
            restart,
            basis: Vec::new(),
            columns: Vec::new(),
            rotations: Vec::new(),
            g: Vec::new(),
            residual: false,
            done: false,
            sweeps: 1,
        };
        let r = k.b.clone();
        k.begin(r);
        k
    }

    /// 残差 `r` から Arnoldi の過程を始める（判定を満たしていれば終える）
    fn begin(&mut self, r: Vec<f64>) {
        let beta = norm(&r);
        self.residual = false;
        if beta <= self.target {
            self.done = true;
            return;
        }
        self.basis = vec![r.par_iter().map(|v| v / beta).collect()];
        self.columns.clear();
        self.rotations.clear();
        self.g = vec![beta];
    }

    /// 作用素を掛ける値（終えていれば `None`）
    fn request(&self) -> Option<&[f64]> {
        if self.done {
            None
        } else if self.residual {
            Some(&self.x)
        } else {
            self.basis.last().map(Vec::as_slice)
        }
    }

    /// `request` の値に作用素を掛けた値 `applied` を受ける（`sweeps` はそこまでのスイープの回数）
    fn receive(&mut self, applied: Vec<f64>, sweeps: usize) {
        if self.residual {
            let r = self
                .b
                .par_iter()
                .zip(&applied)
                .map(|(a, c)| a - c)
                .collect();
            self.begin(r);
            self.sweeps = sweeps;
            return;
        }
        let j = self.columns.len();
        let mut w = applied;
        let mut h = Vec::with_capacity(j + 2);
        for v in &self.basis {
            let c = inner(&w, v);
            w.par_iter_mut().zip(v).for_each(|(a, b)| *a -= c * b);
            h.push(c);
        }
        let next = norm(&w);
        h.push(next);
        for (i, &(c, s)) in self.rotations.iter().enumerate() {
            let (a, b) = (h[i], h[i + 1]);
            h[i] = c * a + s * b;
            h[i + 1] = -s * a + c * b;
        }
        let (a, b) = (h[j], h[j + 1]);
        let radius = a.hypot(b);
        let (c, s) = (a / radius, b / radius);
        h[j] = radius;
        h.truncate(j + 1);
        self.rotations.push((c, s));
        self.g.push(-s * self.g[j]);
        self.g[j] *= c;
        self.columns.push(h);
        let residual = self.g[j + 1].abs();
        tracing::debug!(
            sweeps,
            residual,
            target = self.target,
            "generalized adjoint GMRES"
        );
        if residual <= self.target || next == 0.0 || j + 1 == self.restart {
            // 上三角の系を後退代入で解き、x に基底の和を足す
            let n = self.columns.len();
            let mut y = vec![0.0; n];
            for i in (0..n).rev() {
                let sum: f64 = ((i + 1)..n).map(|l| self.columns[l][i] * y[l]).sum();
                y[i] = (self.g[i] - sum) / self.columns[i][i];
            }
            for (v, coefficient) in self.basis.iter().zip(&y) {
                self.x
                    .par_iter_mut()
                    .zip(v)
                    .for_each(|(a, b)| *a += coefficient * b);
            }
            self.sweeps = sweeps;
            if residual <= self.target || next == 0.0 {
                self.done = true;
            } else {
                self.residual = true;
            }
            return;
        }
        self.basis.push(w.par_iter().map(|v| v / next).collect());
    }
}

fn norm(v: &[f64]) -> f64 {
    inner(v, v).sqrt()
}

fn inner(a: &[f64], b: &[f64]) -> f64 {
    a.par_iter().zip(b).map(|(x, y)| x * y).sum::<f64>()
}

/// 群ごとの、領域の |Γ| の最大値に対する変化の最大（Γ は符号を変えるので、領域ごとの相対の変化は 0 の近くで定まらない）
fn largest_change(new: &[Vec<f64>], old: &[Vec<f64>]) -> f64 {
    let ng = new.first().map_or(0, Vec::len);
    let mut scale = vec![0.0_f64; ng];
    let mut change = vec![0.0_f64; ng];
    for (a, b) in new.iter().zip(old) {
        for g in 0..ng {
            scale[g] = scale[g].max(a[g].abs());
            change[g] = change[g].max((a[g] - b[g]).abs());
        }
    }
    // 群の中性子束が 0 になったのに変化があれば、残りの誤差を見積もれない
    change
        .iter()
        .zip(&scale)
        .map(|(&c, &s)| {
            if c == 0.0 {
                0.0
            } else if s > 0.0 {
                c / s
            } else {
                f64::INFINITY
            }
        })
        .fold(0.0, f64::max)
}

#[cfg(test)]
mod tests;
