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
use crate::cmfd::CmfdOptions;
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
/// 収束した状態（`solve_symmetric_from` と `solve_adjoint_from`）で、k は前向きの値を使う。源の値が欠ける、有限でない、真空で
/// 0 でない、対称操作で移り合う領域で違う、⟨S, φ⟩ が 0 でない（解がない）、状態が別の問題のものなら、エラー
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
        for (r, (q, m)) in source.iter().zip(materials).enumerate() {
            if q.iter().any(|v| !v.is_finite()) || (m.is_void() && q.iter().any(|&v| v != 0.0)) {
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
        materials: &adjoint_materials,
        volumes: &volumes,
        pairs: &pairs,
        ng,
        k: forward.k,
        weights: &weights,
        adjoint,
        adjoint_production,
    };
    reduced_sources
        .iter()
        .map(|q| {
            let (flux, iterations) = match settings.method {
                GeneralizedMethod::SourceIteration => problem.source_iteration(q, settings)?,
                GeneralizedMethod::Cmfd(_) | GeneralizedMethod::Gmres { .. } => {
                    return Err(MocError::GeneralizedNotConverged { iterations: 0 });
                }
            };
            Ok(Generalized {
                flux: reduced.expand(&flux),
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

/// 解く問題（転置した物質のスイープと、φ† を除くための値）
struct Problem<'a> {
    sweeper: &'a Sweeper<'a>,
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
                if let Some(q) = source {
                    for ((e, q), t) in e.iter_mut().zip(&q[i]).zip(&m.total) {
                        if *t > 0.0 {
                            *e += q / (4.0 * PI * t);
                        }
                    }
                }
                e
            })
            .collect()
    }

    /// 平衡値 `equilibrium` でスイープした和 `tallies` から、[代表][群] の中性子束と [代表][モーメント (1 番目から)][群] の高次の
    /// モーメント（`solve_with` と同じ式）
    fn moments(
        &self,
        equilibrium: &[Vec<f64>],
        tallies: &[super::Tally],
    ) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
        let (sweeper, ng, pairs) = (self.sweeper, self.ng, self.pairs);
        let mut sums = sum_tallies(sweeper, tallies);
        sweeper
            .representatives
            .symmetrize(&mut sums, &sweeper.void, ng);
        let kept = self.kept();
        (0..kept.len())
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
            .unzip()
    }

    /// 状態から φ† の状態の倍を除き、⟨Γ, Fφ⟩ = 0 にする
    fn project(&self, flux: &mut [Vec<f64>], higher: &mut [Vec<f64>], outgoing: &mut [f64]) {
        let alpha = dot(self.weights, flux) / self.adjoint_production;
        let a = self.adjoint;
        flux.par_iter_mut()
            .zip(&a.flux)
            .chain(higher.par_iter_mut().zip(&a.higher))
            .for_each(|(x, y)| x.iter_mut().zip(y).for_each(|(x, y)| *x -= alpha * y));
        outgoing
            .par_iter_mut()
            .zip(&a.outgoing)
            .for_each(|(x, y)| *x -= alpha * y);
    }

    /// 射影つきの源の反復: 0 の状態から、源と前の状態の散乱と核分裂の平衡値でスイープし、φ† の成分を除くことを繰り返す。
    /// [代表][群] の中性子束とスイープの回数を返す
    fn source_iteration(
        &self,
        source: &[Vec<f64>],
        settings: &GeneralizedOptions,
    ) -> Result<(Vec<Vec<f64>>, usize), MocError> {
        let (sweeper, ng) = (self.sweeper, self.ng);
        let states = self.kept().len();
        let mut flux = vec![vec![0.0; ng]; states];
        let mut higher = vec![vec![0.0; (self.pairs.len() - 1) * ng]; states];
        let mut outgoing = vec![0.0; sweeper.outgoing_len()];
        let mut next = outgoing.clone();
        let mut tallies = sweeper.tallies(self.materials.len());
        let mut history = Remaining::default();
        for iteration in 1..=settings.max_iterations {
            let equilibrium = self.equilibrium(&flux, &higher, Some(source));
            sweeper.sweep(&equilibrium, &outgoing, &mut next, &mut tallies)?;
            std::mem::swap(&mut outgoing, &mut next);
            let (mut new_flux, mut new_higher) = self.moments(&equilibrium, &tallies);
            self.project(&mut new_flux, &mut new_higher, &mut outgoing);
            let change = largest_change(&new_flux, &flux);
            tracing::debug!(iteration, change, "generalized adjoint");
            flux = new_flux;
            higher = new_higher;
            if history.update(change) < settings.tolerance {
                return Ok((flux, iteration));
            }
        }
        Err(MocError::GeneralizedNotConverged {
            iterations: settings.max_iterations,
        })
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
