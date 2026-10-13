//! 等価なセルの SPH 因子。
//!
//! 超微細群の減速計算（[`crate::ufg`]）のリングごとの群の中性子束 φ_ref を、同じセルを多群の群定数で解いても
//! 再現するよう、リング r・群 g の反応に掛ける因子 μ を求める（Hébert の SPH 等価）。μ を掛けた群定数で解いた
//! 中性子束 φ̃ が φ_ref / μ に比例するよう、μ = φ_ref / φ̃ を反復する。群の釣り合い（反応率の和 = 上の群からの
//! 線源）は μ では変えられず、多群の散乱行列は UFG と同じ群間の移行を与えないので、再現できるのは群ごとの
//! リング間の形である。大きさは、群ごとに因子を掛けるリングの Σ V φ_ref / μ = Σ V φ_ref となるよう規格化して
//! 決める（それらのリングの平均の中性子束を保つ）。真空のリングと、因子を掛ける先のないリング（格子に使わない
//! 等価なセルの枠のリング）は μ = 1 とし、規格化にも入れない。群定数にすると、群の中のスペクトルの空間分布を
//! 平坦にした誤差（文献で約 −200 pcm）を直す。

use std::ops::Range;

use crate::cp::{CpError, Cylinder};
use crate::material::{Material, MaterialError};
use xlc_numeric::lu::solve;

/// SPH 因子の反復の打ち切り（μ の変化）と回数の上限
const TOLERANCE: f64 = 1.0e-10;
const ITERATIONS: usize = 1000;

/// 範囲の中に上方散乱があるときの、群の掃き直しの打ち切り（中性子束の相対変化）
const UPSCATTER_TOLERANCE: f64 = 1.0e-12;

/// SPH 因子の計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum SphError {
    #[error("{given} materials for {rings} rings")]
    Rings { given: usize, rings: usize },
    #[error("{given} fixed flags for {rings} rings")]
    Fixed { given: usize, rings: usize },
    #[error("groups {groups:?} need {needed} groups above and in the materials, got {above} above")]
    Groups {
        groups: Range<usize>,
        needed: usize,
        above: usize,
    },
    #[error("the reference flux must be {rings} rings of {groups} groups")]
    Reference { rings: usize, groups: usize },
    #[error("ring {ring}: the reference flux in group {group} must be positive and finite")]
    ReferenceValue { ring: usize, group: usize },
    #[error("ring {ring}: no multigroup flux in group {group} to match")]
    Flux { ring: usize, group: usize },
    #[error("the multigroup equations of group {0} are singular")]
    Singular(usize),
    #[error("not converged in {0} iterations")]
    NotConverged(usize),
    #[error(transparent)]
    Cp(#[from] CpError),
    #[error(transparent)]
    Material(#[from] MaterialError),
}

/// リングと群の数を確かめる。
fn check(
    cell: &Cylinder,
    materials: &[Material],
    groups: &Range<usize>,
    above: &[f64],
) -> Result<(), SphError> {
    let rings = cell.volumes().len();
    if materials.len() != rings {
        return Err(SphError::Rings {
            given: materials.len(),
            rings,
        });
    }
    let short = materials.iter().any(|m| m.groups() < groups.end);
    if above.len() != groups.start || groups.is_empty() || short {
        return Err(SphError::Groups {
            groups: groups.clone(),
            needed: groups.end,
            above: above.len(),
        });
    }
    Ok(())
}

/// 円筒のセル（白色境界）の、群 `groups` の多群の固定源問題の中性子束 [リング][群（`groups` の中の順）]。
/// `materials` はリングごとの物質、`above` は範囲の上の群（0..groups.start）の中性子束（空間的に平坦）で、
/// その群からの散乱が線源になる。
pub fn multigroup_flux(
    cell: &Cylinder,
    materials: &[Material],
    groups: Range<usize>,
    above: &[f64],
) -> Result<Vec<Vec<f64>>, SphError> {
    check(cell, materials, &groups, above)?;
    let volumes = cell.volumes();
    let n = volumes.len();
    let ng = groups.len();
    // 範囲の上の群からの散乱源 [リング][群]
    let external: Vec<Vec<f64>> = materials
        .iter()
        .map(|m| {
            groups
                .clone()
                .map(|g| (0..groups.start).map(|h| m.scatter[h][g] * above[h]).sum())
                .collect()
        })
        .collect();
    let tracks = groups
        .clone()
        .map(|g| {
            let totals: Vec<f64> = materials.iter().map(|m| m.total[g]).collect();
            cell.tracks(&totals)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let upscatter = materials.iter().any(|m| {
        groups
            .clone()
            .any(|h| groups.clone().any(|g| g < h && m.scatter[h][g] != 0.0))
    });
    let mut flux = vec![vec![0.0; ng]; n];
    for _ in 0..ITERATIONS {
        // 中性子束の相対変化が許容以下か（NaN は収束とみなさない）
        let mut converged = true;
        for (k, g) in groups.clone().enumerate() {
            // 他の群からの散乱源と、自群散乱の係数
            let source: Vec<f64> = (0..n)
                .map(|i| {
                    external[i][k]
                        + groups
                            .clone()
                            .enumerate()
                            .filter(|&(j, _)| j != k)
                            .map(|(j, h)| materials[i].scatter[h][g] * flux[i][j])
                            .sum::<f64>()
                })
                .collect();
            let w = &tracks[k];
            // φ_j V_j − Σ_i W_ij Σs,i(g→g) φ_i = Σ_i W_ij Q_i（W は対称）
            let matrix: Vec<f64> = (0..n * n)
                .map(|k| {
                    let (j, i) = (k / n, k % n);
                    let diagonal = if i == j { volumes[j] } else { 0.0 };
                    diagonal - w[j][i] * materials[i].scatter[g][g]
                })
                .collect();
            let rhs: Vec<f64> = (0..n)
                .map(|j| (0..n).map(|i| w[j][i] * source[i]).sum())
                .collect();
            let phi = solve(matrix, rhs).ok_or(SphError::Singular(g))?;
            for (i, p) in phi.into_iter().enumerate() {
                converged &= (p - flux[i][k]).abs() <= UPSCATTER_TOLERANCE * p.abs();
                flux[i][k] = p;
            }
        }
        if !upscatter || converged {
            return Ok(flux);
        }
    }
    Err(SphError::NotConverged(ITERATIONS))
}

/// 参照の中性子束 `reference` [リング][群] を再現する SPH 因子 [リング][群（`groups` の中の順）]。`fixed` の
/// リングは μ = 1 のままにする。
pub fn sph_factors(
    cell: &Cylinder,
    materials: &[Material],
    fixed: &[bool],
    groups: Range<usize>,
    above: &[f64],
    reference: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, SphError> {
    check(cell, materials, &groups, above)?;
    let n = materials.len();
    let ng = groups.len();
    if reference.len() != n || reference.iter().any(|r| r.len() != ng) {
        return Err(SphError::Reference {
            rings: n,
            groups: ng,
        });
    }
    if fixed.len() != n {
        return Err(SphError::Fixed {
            given: fixed.len(),
            rings: n,
        });
    }
    let volumes = cell.volumes();
    let active: Vec<usize> = (0..n)
        .filter(|&r| !(materials[r].is_void() || fixed[r]))
        .collect();
    for &r in &active {
        if let Some(k) = reference[r]
            .iter()
            .position(|&f| !(f > 0.0 && f.is_finite()))
        {
            return Err(SphError::ReferenceValue {
                ring: r,
                group: groups.start + k,
            });
        }
    }
    let mut factors = vec![vec![1.0; ng]; n];
    for _ in 0..ITERATIONS {
        let scaled = materials
            .iter()
            .zip(&factors)
            .map(|(m, mu)| {
                let mut all = vec![1.0; m.groups()];
                all[groups.clone()].copy_from_slice(mu);
                m.scaled(&all)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let flux = multigroup_flux(cell, &scaled, groups.clone(), above)?;
        let mut next = vec![vec![1.0; ng]; n];
        for k in 0..ng {
            for &r in &active {
                if !(flux[r][k] > 0.0 && flux[r][k].is_finite()) {
                    return Err(SphError::Flux {
                        ring: r,
                        group: groups.start + k,
                    });
                }
                next[r][k] = reference[r][k] / flux[r][k];
            }
            // Σ V φ_ref / μ = Σ V φ_ref に規格化する
            let preserved: f64 = active.iter().map(|&r| volumes[r] * reference[r][k]).sum();
            let current: f64 = active
                .iter()
                .map(|&r| volumes[r] * reference[r][k] / next[r][k])
                .sum();
            for &r in &active {
                next[r][k] *= current / preserved;
            }
        }
        // NaN は収束とみなさない（f64::max は NaN を無視する）
        let converged = next
            .iter()
            .zip(&factors)
            .flat_map(|(a, b)| a.iter().zip(b))
            .all(|(x, y)| (x - y).abs() <= TOLERANCE);
        factors = next;
        if converged {
            return Ok(factors);
        }
    }
    Err(SphError::NotConverged(ITERATIONS))
}

#[cfg(test)]
mod tests;
