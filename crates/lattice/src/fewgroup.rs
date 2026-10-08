//! MOC の解から、少数群定数（`xlc_format::fewgroup::FewGroup`）を作る。

use xlc_format::fewgroup::FewGroup;
use xlc_format::group::GroupStructure;

use crate::material::Material;
use crate::moc::Solution;

/// 少数群への縮約のエラー。
#[derive(Debug, thiserror::Error)]
pub enum HomogenizeError {
    #[error("{materials} materials for {regions} regions")]
    RegionCount { materials: usize, regions: usize },
    #[error("few-group bounds must be bounds of the library groups")]
    NotLibraryBounds,
    #[error("few-group bounds must be descending and span the library groups")]
    Span,
    #[error("a few group has no flux")]
    NoFlux,
}

/// MOC の解を使い、領域・群ごとの反応率を保存するようにセル全体で均質化し、
/// `coarse_bounds`（`fine` の境界の一部、降順）の少数群に縮約する。
/// 拡散係数は D = 1/(3Σt) とする（M0 の近似）。
pub fn homogenize(
    solution: &Solution,
    materials: &[&Material],
    fine: &GroupStructure,
    coarse_bounds: &[f64],
) -> Result<FewGroup, HomogenizeError> {
    let ng = fine.len();
    if materials.len() != solution.flux.len() {
        return Err(HomogenizeError::RegionCount {
            materials: materials.len(),
            regions: solution.flux.len(),
        });
    }
    // 少数群の境界を細かい群の境界の位置に対応させる
    let positions = coarse_bounds
        .iter()
        .map(|b| fine.bounds().iter().position(|f| f == b))
        .collect::<Option<Vec<usize>>>()
        .ok_or(HomogenizeError::NotLibraryBounds)?;
    let covers = positions.first() == Some(&0) && positions.last() == Some(&ng);
    if !covers || positions.windows(2).any(|w| w[0] >= w[1]) {
        return Err(HomogenizeError::Span);
    }
    let coarse = positions.len() - 1;
    let mut group_of = vec![0; ng];
    for (big, w) in positions.windows(2).enumerate() {
        group_of[w[0]..w[1]].fill(big);
    }

    let mut flux = vec![0.0; coarse];
    let mut total = vec![0.0; coarse];
    let mut nu_fission = vec![0.0; coarse];
    let mut fission = vec![0.0; coarse];
    let mut scatter = vec![vec![0.0; coarse]; coarse];
    let mut chi = vec![0.0; coarse];
    for ((m, phi), volume) in materials.iter().zip(&solution.flux).zip(&solution.volumes) {
        for g in 0..ng {
            let w = volume * phi[g];
            let big = group_of[g];
            flux[big] += w;
            total[big] += w * m.total[g];
            nu_fission[big] += w * m.nu_fission[g];
            fission[big] += w * m.fission[g];
            for (to, s) in m.scatter[g].iter().enumerate() {
                scatter[big][group_of[to]] += w * s;
            }
        }
        for term in &m.fission_terms {
            let rate: f64 = term
                .nu_fission
                .iter()
                .zip(phi)
                .map(|(n, f)| volume * n * f)
                .sum();
            for (g, c) in term.chi.iter().enumerate() {
                chi[group_of[g]] += c * rate;
            }
        }
    }
    if flux.iter().any(|&f| f <= 0.0) {
        return Err(HomogenizeError::NoFlux);
    }
    let per_flux =
        |values: Vec<f64>| -> Vec<f64> { values.iter().zip(&flux).map(|(v, f)| v / f).collect() };
    let total = per_flux(total);
    let scatter: Vec<Vec<f64>> = scatter
        .into_iter()
        .zip(&flux)
        .map(|(row, f)| row.iter().map(|v| v / f).collect())
        .collect();
    let produced: f64 = chi.iter().sum();
    Ok(FewGroup {
        bounds: coarse_bounds.to_vec(),
        k_inf: solution.k,
        absorption: total
            .iter()
            .zip(&scatter)
            .map(|(t, row)| t - row.iter().sum::<f64>())
            .collect(),
        diffusion: total.iter().map(|t| 1.0 / (3.0 * t)).collect(),
        total,
        nu_fission: per_flux(nu_fission),
        fission: per_flux(fission),
        chi: if produced > 0.0 {
            chi.iter().map(|c| c / produced).collect()
        } else {
            chi
        },
        scatter,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Lattice, PinCell};
    use crate::material::FissionTerm;
    use crate::moc::{Options, solve};
    use crate::tracking::Tracking;

    fn two_group(total: [f64; 2], scatter: [[f64; 2]; 2], nu_fission: [f64; 2]) -> Material {
        Material {
            total: total.to_vec(),
            scatter: scatter.iter().map(|r| r.to_vec()).collect(),
            anisotropic: Vec::new(),
            fission: nu_fission.iter().map(|v| v / 2.5).collect(),
            nu_fission: nu_fission.to_vec(),
            fission_terms: if nu_fission.iter().any(|&v| v > 0.0) {
                vec![FissionTerm {
                    chi: vec![1.0, 0.0],
                    nu_fission: nu_fission.to_vec(),
                }]
            } else {
                vec![]
            },
        }
    }

    fn solved() -> (Solution, Vec<Material>, GroupStructure) {
        let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
        let moderator = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
        let materials = vec![fuel.clone(), moderator.clone(), moderator];
        let cell = PinCell::new(1.26, vec![0.4096, 0.55]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 16, 0.05).unwrap();
        let refs: Vec<&Material> = materials.iter().collect();
        let solution = solve(&tracking, &refs, &Options::default()).unwrap();
        (
            solution,
            materials,
            GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap(),
        )
    }

    #[test]
    fn one_group_collapse_preserves_k() {
        // 反応率を保存するので、1群の νΣf/Σa はセルの k∞ に一致する
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let one = homogenize(&solution, &refs, &fine, &[1.0e7, 1.0e-5]).unwrap();
        assert_eq!(one.groups(), 1);
        let k = one.nu_fission[0] / one.absorption[0];
        assert!(
            (k - solution.k).abs() < 1e-6 * solution.k,
            "{k} vs {}",
            solution.k
        );
        assert_eq!(one.k_inf, solution.k);
        assert!((one.chi[0] - 1.0).abs() < 1e-12);
        assert!((one.diffusion[0] - 1.0 / (3.0 * one.total[0])).abs() < 1e-15);
    }

    #[test]
    fn two_group_constants_reproduce_k_inf() {
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let fg = homogenize(&solution, &refs, &fine, fine.bounds()).unwrap();
        // 均質な2群の無限媒質の k = νΣfᵀ (diag Σt - Σsᵀ)⁻¹ χ（χ = (1, 0)）
        let a = [
            [fg.total[0] - fg.scatter[0][0], -fg.scatter[1][0]],
            [-fg.scatter[0][1], fg.total[1] - fg.scatter[1][1]],
        ];
        let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        let phi = [a[1][1] / det, -a[1][0] / det];
        let k = fg.nu_fission[0] * phi[0] + fg.nu_fission[1] * phi[1];
        assert!((k - solution.k).abs() < 1e-6 * solution.k);
        for g in 0..2 {
            let removal = fg.total[g] - fg.scatter[g].iter().sum::<f64>();
            assert!((fg.absorption[g] - removal).abs() < 1e-14);
        }
    }

    #[test]
    fn coarse_bounds_must_be_fine_bounds() {
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        assert!(homogenize(&solution, &refs, &fine, &[1.0e7, 1.0, 1.0e-5]).is_err());
        assert!(homogenize(&solution, &refs, &fine, &[1.0e7, 0.625]).is_err());
    }
}
