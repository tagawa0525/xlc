//! 等価なセルの SPH 因子。
//!
//! 超微細群の減速計算（[`crate::ufg`]）のリングごとの群の中性子束 φ_ref を、同じセルを多群の群定数で解いても
//! 再現するよう、リング r・群 g の反応に掛ける因子 μ を求める（Hébert の SPH 等価）。μ を掛けた群定数で解いた
//! 中性子束 φ̃ が φ_ref / μ に比例するよう、μ = φ_ref / φ̃ を反復する。群のバランス（反応率の和 = 上の群からの
//! 線源）は μ では変えられず、多群の散乱行列は UFG と同じ群間の移行を与えないので、再現できるのは群ごとの
//! リング間の形である。大きさは、群ごとに真空でないリングの Σ V φ_ref / μ = Σ V φ_ref となるよう規格化して
//! 決める（セルの平均の中性子束を保つ）。真空のリングには掛ける断面積がないので μ = 1 とする。群定数にすると、
//! 群の中のスペクトルの空間分布を平坦にした誤差（文献で約 −200 pcm）を直す。

use std::ops::Range;

use crate::cp::{CpError, Cylinder};
use crate::material::{Material, MaterialError};
use crate::ufg::solve;

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
    #[error("groups {groups:?} need {needed} groups above and in the materials, got {above} above")]
    Groups {
        groups: Range<usize>,
        needed: usize,
        above: usize,
    },
    #[error("the reference flux must be {rings} rings of {groups} groups")]
    Reference { rings: usize, groups: usize },
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
        let mut change: f64 = 0.0;
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
            let matrix: Vec<Vec<f64>> = (0..n)
                .map(|j| {
                    (0..n)
                        .map(|i| {
                            let diagonal = if i == j { volumes[j] } else { 0.0 };
                            diagonal - w[j][i] * materials[i].scatter[g][g]
                        })
                        .collect()
                })
                .collect();
            let rhs: Vec<f64> = (0..n)
                .map(|j| (0..n).map(|i| w[j][i] * source[i]).sum())
                .collect();
            let phi = solve(matrix, rhs).ok_or(SphError::Singular(g))?;
            for (i, p) in phi.into_iter().enumerate() {
                change = change.max(((p - flux[i][k]) / p).abs());
                flux[i][k] = p;
            }
        }
        if !upscatter || change <= UPSCATTER_TOLERANCE {
            return Ok(flux);
        }
    }
    Err(SphError::NotConverged(ITERATIONS))
}

/// 参照の中性子束 `reference` [リング][群] を再現する SPH 因子 [リング][群（`groups` の中の順）]。
pub fn sph_factors(
    cell: &Cylinder,
    materials: &[Material],
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
    let volumes = cell.volumes();
    let active: Vec<usize> = (0..n).filter(|&r| !materials[r].is_void()).collect();
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
        let change = next
            .iter()
            .zip(&factors)
            .flat_map(|(a, b)| a.iter().zip(b).map(|(x, y)| (x - y).abs()))
            .fold(0.0, f64::max);
        factors = next;
        if change <= TOLERANCE {
            return Ok(factors);
        }
    }
    Err(SphError::NotConverged(ITERATIONS))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::Material;

    /// 群 0 が範囲の上、群 1〜3 が範囲の中
    const GROUPS: std::ops::Range<usize> = 1..4;
    const RADII: [f64; 4] = [0.3, 0.4096, 0.418, 0.7];

    fn material(total: [f64; 4], scattering: f64) -> Material {
        // 散乱は自群と次の群へ（下方散乱だけ）
        let scatter = (0..4)
            .map(|g| {
                let mut row = vec![0.0; 4];
                row[g] = 0.6 * scattering * total[g];
                if g + 1 < 4 {
                    row[g + 1] = 0.4 * scattering * total[g];
                }
                row
            })
            .collect();
        Material {
            total: total.to_vec(),
            scatter,
            fission: vec![0.0; 4],
            nu_fission: vec![0.0; 4],
            fission_terms: Vec::new(),
        }
    }

    fn materials(fuel_total: f64) -> Vec<Material> {
        let fuel = material([0.5, fuel_total, 1.0, 0.8], 0.3);
        let water = material([0.6, 1.0, 1.1, 1.2], 0.99);
        vec![fuel.clone(), fuel, Material::void(4), water]
    }

    #[test]
    fn consistent_reference_needs_no_correction() {
        let cell = Cylinder::new(RADII.to_vec()).unwrap();
        let reference = multigroup_flux(&cell, &materials(2.0), GROUPS, &[1.0]).unwrap();
        let factors = sph_factors(&cell, &materials(2.0), GROUPS, &[1.0], &reference).unwrap();
        for row in &factors {
            for mu in row {
                assert!((mu - 1.0).abs() < 1e-10, "{mu}");
            }
        }
    }

    #[test]
    fn factors_reproduce_the_reference_shape() {
        // 燃料の全断面積を変えた解を参照にする。群のバランス（反応率の和 = 上の群からの線源）は μ では変えられない
        // ので、群ごとのリング間の形 μ φ̃ ∝ φ_ref を再現し、大きさは Σ V φ_ref / μ = Σ V φ_ref で決める
        let cell = Cylinder::new(RADII.to_vec()).unwrap();
        let reference = multigroup_flux(&cell, &materials(3.0), GROUPS, &[1.0]).unwrap();
        let original = materials(2.0);
        let factors = sph_factors(&cell, &original, GROUPS, &[1.0], &reference).unwrap();
        let scaled: Vec<Material> = original
            .iter()
            .zip(&factors)
            .map(|(m, mu)| {
                let mut all = vec![1.0; 4];
                all[GROUPS].copy_from_slice(mu);
                m.scaled(&all).unwrap()
            })
            .collect();
        let flux = multigroup_flux(&cell, &scaled, GROUPS, &[1.0]).unwrap();
        let volumes = cell.volumes();
        for k in 0..GROUPS.len() {
            // 真空でないリングで μ φ̃ / φ_ref が等しい
            let ratio = |r: usize| factors[r][k] * flux[r][k] / reference[r][k];
            for r in [1, 3] {
                assert!(
                    (ratio(r) - ratio(0)).abs() < 1e-8 * ratio(0),
                    "{r} {k}: {} vs {}",
                    ratio(r),
                    ratio(0)
                );
            }
            let preserved: f64 = [0, 1, 3]
                .iter()
                .map(|&r| volumes[r] * reference[r][k])
                .sum();
            let normalized: f64 = [0, 1, 3]
                .iter()
                .map(|&r| volumes[r] * reference[r][k] / factors[r][k])
                .sum();
            assert!((normalized - preserved).abs() < 1e-10 * preserved);
            // 真空のリングには掛ける断面積がないので 1
            assert_eq!(factors[2][k], 1.0);
        }
        // 燃料の全断面積を大きくした参照なので、燃料の μ は 1 と違う
        assert!((factors[0][0] - 1.0).abs() > 1e-3, "{:?}", factors[0]);
    }

    #[test]
    fn rejects_rings_without_flux() {
        // 範囲の上から散乱してくる中性子がなければ、多群の中性子束が 0 になり、μ が決まらない
        let cell = Cylinder::new(RADII.to_vec()).unwrap();
        let reference = multigroup_flux(&cell, &materials(2.0), GROUPS, &[1.0]).unwrap();
        let mut isolated = materials(2.0);
        for m in &mut isolated {
            m.scatter[0] = vec![0.0; 4];
        }
        assert!(matches!(
            sph_factors(&cell, &isolated, GROUPS, &[1.0], &reference),
            Err(SphError::Flux { .. })
        ));
    }

    #[test]
    fn rejects_mismatched_sizes() {
        let cell = Cylinder::new(RADII.to_vec()).unwrap();
        assert!(multigroup_flux(&cell, &materials(2.0)[..3], GROUPS, &[1.0]).is_err());
        assert!(multigroup_flux(&cell, &materials(2.0), GROUPS, &[]).is_err());
        let reference = vec![vec![1.0; 2]; 4];
        assert!(sph_factors(&cell, &materials(2.0), GROUPS, &[1.0], &reference).is_err());
        // 物質の群より先まで求める群の範囲
        let reference = vec![vec![1.0; 5]; 4];
        assert!(sph_factors(&cell, &materials(2.0), 1..6, &[1.0], &reference).is_err());
    }
}
