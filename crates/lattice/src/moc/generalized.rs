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
}

/// 1 つの源の一般化随伴
#[derive(Debug, Clone, PartialEq)]
pub struct Generalized {
    /// [領域][群] のスカラー中性子束（源 S と同じ規格化。⟨Γ, Fφ⟩ = 0）
    pub flux: Vec<Vec<f64>>,
    /// スイープの回数
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
        sweeper: &sweeper,
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
    reduced_sources
        .iter()
        .map(|q| {
            let (state, iterations) = match (settings.method, &coarse) {
                (GeneralizedMethod::Gmres { restart }, _) => problem.gmres(q, restart, settings)?,
                (_, coarse) => problem.source_iteration(q, coarse.as_ref(), settings)?,
            };
            Ok(Generalized {
                flux: reduced.expand(&state.flux),
                iterations,
            })
        })
        .collect()
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
    sweeper: &'a Sweeper<'a>,
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

impl Problem<'_> {
    fn kept(&self) -> &[usize] {
        &self.sweeper.representatives.regions.representatives
    }

    /// 0 の状態
    fn zero(&self) -> State {
        let states = self.kept().len();
        State {
            flux: vec![vec![0.0; self.ng]; states],
            higher: vec![vec![0.0; (self.pairs.len() - 1) * self.ng]; states],
            outgoing: vec![0.0; self.sweeper.outgoing_len()],
        }
    }

    /// 状態（[代表][群] の中性子束と高次のモーメント）と源 `source`（[代表][群]。`None` なら 0）の、[代表][モーメント][群] の
    /// 角度中性子束の平衡値の展開係数
    fn equilibrium(
        &self,
        flux: &[Vec<f64>],
        higher: &[Vec<f64>],
        source: Option<&[Vec<f64>]>,
    ) -> Vec<Vec<f64>> {
        let kept = self.kept();
        (0..kept.len())
            .into_par_iter()
            .map(|i| {
                let m = self.materials[kept[i]];
                let mut e = equilibrium_moments(m, &flux[i], &higher[i], self.k, self.pairs);
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

    /// 平衡値 `equilibrium` でスイープした和 `tallies` から、[代表][群] の中性子束と [代表][モーメント (1 番目から)][群] の高次の
    /// モーメント（`solve_with` と同じ式）と、[粗メッシュの面][群] の正味の中性子流
    fn moments(&self, equilibrium: &[Vec<f64>], tallies: &[super::Tally]) -> (State, Vec<f64>) {
        let (sweeper, ng, pairs) = (self.sweeper, self.ng, self.pairs);
        let mut sums = sum_tallies(sweeper, tallies);
        sweeper
            .representatives
            .symmetrize(&mut sums, &sweeper.void, ng);
        let kept = self.kept();
        let (flux, higher) = (0..kept.len())
            .into_par_iter()
            .map(|i| {
                let region = kept[i];
                let volume = self.volumes[region];
                let mut higher = vec![0.0; (pairs.len() - 1) * ng];
                if let Some(v) = sweeper.void[region] {
                    let flux = sums.void_flux[v * ng..][..ng]
                        .iter()
                        .map(|f| f / volume)
                        .collect();
                    return (flux, higher);
                }
                let total = &self.materials[region].total;
                let flux = (0..ng)
                    .map(|g| {
                        4.0 * PI * equilibrium[i][g]
                            + sums.delta[i * pairs.len() * ng + g] / (total[g] * volume)
                    })
                    .collect();
                for (k, &(l, _)) in pairs.iter().enumerate().skip(1) {
                    let delta = &sums.delta[(i * pairs.len() + k) * ng..][..ng];
                    let e = &equilibrium[i][k * ng..][..ng];
                    for (g, h) in higher[(k - 1) * ng..][..ng].iter_mut().enumerate() {
                        *h = 4.0 * PI * e[g] / (2 * l + 1) as f64 + delta[g] / (total[g] * volume);
                    }
                }
                (flux, higher)
            })
            .unzip();
        (
            State {
                flux,
                higher,
                outgoing: Vec::new(),
            },
            sums.current,
        )
    }

    /// 状態 `state` から、源 `source` を足して 1 回スイープした状態（φ† の成分を除いたもの）と、面の正味の中性子流
    fn step(
        &self,
        state: &State,
        source: Option<&[Vec<f64>]>,
        tallies: &mut [super::Tally],
    ) -> Result<(State, Vec<f64>), MocError> {
        let equilibrium = self.equilibrium(&state.flux, &state.higher, source);
        let mut outgoing = vec![0.0; state.outgoing.len()];
        self.sweeper
            .sweep(&equilibrium, &state.outgoing, &mut outgoing, tallies)?;
        let (mut next, current) = self.moments(&equilibrium, tallies);
        next.outgoing = outgoing;
        self.project(&mut next);
        Ok((next, current))
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
    /// `coarse` があれば、スイープごとに粗メッシュの補正を足す（`correct`）。状態とスイープの回数を返す
    fn source_iteration(
        &self,
        source: &[Vec<f64>],
        coarse: Option<&CoarseAcceleration>,
        settings: &GeneralizedOptions,
    ) -> Result<(State, usize), MocError> {
        let mut state = self.zero();
        let mut tallies = self.sweeper.tallies(self.materials.len());
        let mut history = Remaining::default();
        for iteration in 1..=settings.max_iterations {
            let (mut next, _) = self.step(&state, Some(source), &mut tallies)?;
            if let Some(c) = coarse {
                self.correct(c, &state, &mut next)?;
            }
            let change = largest_change(&next.flux, &state.flux);
            tracing::debug!(iteration, change, "generalized adjoint");
            state = next;
            if history.update(change) < settings.tolerance {
                return Ok((state, iteration));
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
        let mut tallies = self.sweeper.tallies(self.materials.len());
        let equilibrium = self.equilibrium(&a.flux, &a.higher, None);
        let mut outgoing = vec![0.0; a.outgoing.len()];
        self.sweeper
            .sweep(&equilibrium, &a.outgoing, &mut outgoing, &mut tallies)?;
        let (swept, mut current) = self.moments(&equilibrium, &tallies);
        let reduced = &self.sweeper.representatives;
        let orbit_volumes = reduced.regions.sums(self.volumes);
        let kept = self.kept();
        let rate: f64 = swept
            .flux
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
            .flux
            .iter()
            .map(|f| f.iter().map(|x| x / rate).collect())
            .collect();
        for c in &mut current {
            *c /= rate;
        }
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
        let ng = self.ng;
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
        let reduced = &self.sweeper.representatives;
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
        self.sweeper
            .add_outgoing(&mut after.outgoing, &a.outgoing, |cell, g| {
                ratio[cells.index[cell] * ng + g]
            });
        self.project(after);
        Ok(())
    }

    /// 射影した作用素の GMRES: x = P T₀ x + P T(0)（T₀ は源のないスイープ、P は φ† の成分を除く射影）を (I − P T₀) x = P T(0) と
    /// 書き、`restart` 回ごとに再始動する GMRES で解く。内積は並べた状態の値の和で、残差のノルムが右辺のノルムの判定倍に
    /// なったら止める。状態とスイープの回数を返す
    fn gmres(
        &self,
        source: &[Vec<f64>],
        restart: usize,
        settings: &GeneralizedOptions,
    ) -> Result<(State, usize), MocError> {
        let zero = self.zero();
        let mut tallies = self.sweeper.tallies(self.materials.len());
        let (first, _) = self.step(&zero, Some(source), &mut tallies)?;
        let mut sweeps = 1;
        let b = first.pack();
        let norm = |v: &[f64]| v.par_iter().map(|x| x * x).sum::<f64>().sqrt();
        let inner = |a: &[f64], b: &[f64]| a.par_iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
        let target = settings.tolerance * norm(&b);
        let mut x = vec![0.0; b.len()];
        // (I − P T₀) v
        let mut apply = |v: &[f64], sweeps: &mut usize| -> Result<Vec<f64>, MocError> {
            *sweeps += 1;
            let (image, _) = self.step(&State::unpack(v, &zero), None, &mut tallies)?;
            let image = image.pack();
            Ok(v.par_iter().zip(&image).map(|(a, b)| a - b).collect())
        };
        if target == 0.0 {
            return Ok((zero, sweeps));
        }
        loop {
            let r: Vec<f64> = if sweeps == 1 {
                b.clone()
            } else {
                let ax = apply(&x, &mut sweeps)?;
                b.par_iter().zip(&ax).map(|(a, c)| a - c).collect()
            };
            let beta = norm(&r);
            if beta <= target {
                return Ok((State::unpack(&x, &zero), sweeps));
            }
            let mut basis: Vec<Vec<f64>> = vec![r.par_iter().map(|v| v / beta).collect()];
            // 上三角にした Hessenberg 行列の列と、Givens の回転、最小二乗の右辺
            let mut columns: Vec<Vec<f64>> = Vec::new();
            let mut rotations: Vec<(f64, f64)> = Vec::new();
            let mut g = vec![beta];
            for j in 0..restart {
                if sweeps >= settings.max_iterations {
                    return Err(MocError::GeneralizedNotConverged { iterations: sweeps });
                }
                let mut w = apply(&basis[j], &mut sweeps)?;
                let mut h = Vec::with_capacity(j + 2);
                for v in &basis {
                    let c = inner(&w, v);
                    w.par_iter_mut().zip(v).for_each(|(a, b)| *a -= c * b);
                    h.push(c);
                }
                let next = norm(&w);
                h.push(next);
                for (i, &(c, s)) in rotations.iter().enumerate() {
                    let (a, b) = (h[i], h[i + 1]);
                    h[i] = c * a + s * b;
                    h[i + 1] = -s * a + c * b;
                }
                let (a, b) = (h[j], h[j + 1]);
                let radius = a.hypot(b);
                let (c, s) = (a / radius, b / radius);
                h[j] = radius;
                h[j + 1] = 0.0;
                rotations.push((c, s));
                g.push(-s * g[j]);
                g[j] *= c;
                h.truncate(j + 1);
                columns.push(h);
                let converged = g[j + 1].abs() <= target;
                tracing::debug!(
                    sweeps,
                    residual = g[j + 1].abs() / (target / settings.tolerance),
                    "generalized adjoint GMRES"
                );
                if converged || next == 0.0 || j + 1 == restart {
                    // 上三角の系を後退代入で解き、x に基底の和を足す
                    let n = columns.len();
                    let mut y = vec![0.0; n];
                    for i in (0..n).rev() {
                        let sum: f64 = ((i + 1)..n).map(|l| columns[l][i] * y[l]).sum();
                        y[i] = (g[i] - sum) / columns[i][i];
                    }
                    for (v, coefficient) in basis.iter().zip(&y) {
                        x.par_iter_mut()
                            .zip(v)
                            .for_each(|(a, b)| *a += coefficient * b);
                    }
                    break;
                }
                basis.push(w.par_iter().map(|v| v / next).collect());
            }
        }
    }
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
