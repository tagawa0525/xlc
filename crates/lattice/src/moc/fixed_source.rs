//! 散乱のない固定源の MOC（ESSM のバックグラウンド断面積を求めるのに使う）。
//!
//! Ω·∇ψ + Σt ψ = Q / 4π を、領域・群ごとの全断面積 Σt と等方線源 Q で解く。散乱がないので、反復するのは
//! 境界（鏡面反射の辺）の角度中性子束だけである。スカラー中性子束は、Q = Σt の無限媒質で 1 になる規格化。
//! 全断面積が 0 の領域は真空、無限大の領域は黒体（中性子束は 0）として扱う。

use std::f64::consts::PI;

use super::{Representatives, Sweeper, sum_tallies};
use crate::cmfd::{CmfdError, CoarseMesh};
use crate::material::Material;
use crate::tracking::Tracking;

/// 反復の打ち切り条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FixedSourceOptions {
    /// 境界の角度中性子束の変化の上限（群ごとの、そのときの角度中性子束の最大値に対する比）
    pub tolerance: f64,
    pub max_iterations: usize,
}

impl Default for FixedSourceOptions {
    fn default() -> Self {
        Self {
            tolerance: 1.0e-10,
            max_iterations: 10_000,
        }
    }
}

/// 固定源の計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum FixedSourceError {
    #[error("{given} values for {regions} regions")]
    Regions { given: usize, regions: usize },
    #[error("region {0}: totals and sources need the same number of groups")]
    Groups(usize),
    #[error("region {0}: total cross sections must be positive (or zero in every group for void)")]
    Total(usize),
    #[error(
        "region {0}: sources must be finite and not negative, and zero in void and black regions"
    )]
    Source(usize),
    #[error("no track crosses region {0}; use a finer track spacing")]
    Untracked(usize),
    #[error("boundary fluxes not converged in {0} iterations")]
    NotConverged(usize),
    #[error(transparent)]
    Cmfd(#[from] CmfdError),
}

/// [領域][群] のスカラー中性子束。`totals` と `sources` は [領域][群] の全断面積 (1/cm) と等方線源。
pub fn fixed_source(
    tracking: &Tracking,
    totals: &[Vec<f64>],
    sources: &[Vec<f64>],
    options: &FixedSourceOptions,
) -> Result<Vec<Vec<f64>>, FixedSourceError> {
    let regions = tracking.lattice.volumes().len();
    for given in [totals.len(), sources.len()] {
        if given != regions {
            return Err(FixedSourceError::Regions { given, regions });
        }
    }
    let ng = totals.first().map_or(0, Vec::len);
    for (region, (t, q)) in totals.iter().zip(sources).enumerate() {
        if t.len() != ng || q.len() != ng {
            return Err(FixedSourceError::Groups(region));
        }
        let void = t.iter().all(|&v| v == 0.0);
        if !void && !t.iter().all(|&v| v > 0.0 && !v.is_nan()) {
            return Err(FixedSourceError::Total(region));
        }
        // 真空と黒体には線源を置けない（平衡値 Q/(4πΣt) が定まらないか 0 になる）
        let sourceless = |(tv, qv): (&f64, &f64)| (*tv == 0.0 || tv.is_infinite()) && *qv != 0.0;
        if q.iter().any(|v| !(*v >= 0.0 && v.is_finite())) || t.iter().zip(q).any(sourceless) {
            return Err(FixedSourceError::Source(region));
        }
    }

    let materials: Vec<Material> = totals
        .iter()
        .map(|t| Material {
            total: t.clone(),
            ..Material::void(ng)
        })
        .collect();
    let references: Vec<&Material> = materials.iter().collect();
    let mesh = CoarseMesh::new(tracking);
    let sweeper = Sweeper::new(
        tracking,
        &references,
        ng,
        &mesh,
        &[(0, 0)],
        Representatives::all(tracking, &mesh, 1),
    );
    let volumes = tracking.volumes(regions);
    if let Some(region) = volumes.iter().position(|&v| v <= 0.0) {
        return Err(FixedSourceError::Untracked(region));
    }
    // 角度中性子束の平衡値 Q / (4πΣt)
    let equilibrium: Vec<Vec<f64>> = totals
        .iter()
        .zip(sources)
        .map(|(t, q)| {
            t.iter()
                .zip(q)
                .map(|(tv, qv)| {
                    if *qv == 0.0 {
                        0.0
                    } else {
                        qv / (4.0 * PI * tv)
                    }
                })
                .collect()
        })
        .collect();
    // 最初は端の領域の平衡値が等方に入ってくるとする
    let start: Vec<Vec<f64>> = equilibrium
        .iter()
        .map(|q| q.iter().map(|v| 4.0 * PI * v).collect())
        .collect();
    let mut outgoing = sweeper.isotropic_outgoing(&start);
    let mut next = outgoing.clone();
    let mut tallies = sweeper.tallies(regions);
    for _ in 0..options.max_iterations {
        sweeper.sweep(&equilibrium, &outgoing, &mut next, &mut tallies)?;
        // 群ごとに、その群の角度中性子束の最大値に対する変化で判定する（群の間で大きさが何桁も
        // 違っても、小さな群を収束しないまま打ち切らないように）。並びは群がいちばん速く変わる
        let mut scale = vec![0.0_f64; ng];
        let mut change = vec![0.0_f64; ng];
        for (i, (a, b)) in outgoing.iter().zip(&next).enumerate() {
            scale[i % ng] = scale[i % ng].max(b.abs());
            change[i % ng] = change[i % ng].max((a - b).abs());
        }
        std::mem::swap(&mut outgoing, &mut next);
        if change
            .iter()
            .zip(&scale)
            .all(|(c, s)| *c <= options.tolerance * s)
        {
            let sums = sum_tallies(&sweeper, &tallies);
            return Ok((0..regions)
                .map(|i| match sweeper.void[i] {
                    Some(v) => sums.void_flux[v * ng..][..ng]
                        .iter()
                        .map(|f| f / volumes[i])
                        .collect(),
                    None => (0..ng)
                        .map(|g| {
                            let total = totals[i][g];
                            if total.is_infinite() {
                                return 0.0;
                            }
                            4.0 * PI * equilibrium[i][g]
                                + sums.delta[i * ng + g] / (total * volumes[i])
                        })
                        .collect(),
                })
                .collect());
        }
    }
    Err(FixedSourceError::NotConverged(options.max_iterations))
}

#[cfg(test)]
mod tests;
