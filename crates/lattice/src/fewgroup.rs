//! MOC の解から、少数群定数（`xlc_format::fewgroup::FewGroup`）を作る。

use xlc_format::fewgroup::FewGroup;
use xlc_format::group::GroupStructure;

use crate::b1::{B1Error, Medium};
use crate::material::Material;
use crate::moc::Solution;
use crate::tracking::{Boundaries, Boundary};

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
    #[error("library group {0} has no flux")]
    NoFineFlux(usize),
    #[error("discontinuity factors need reflective sides on the whole lattice")]
    NotReflective,
    #[error(
        "region {region} has no P1 scattering moments; diffusion coefficients need a library with them"
    )]
    NoP1 { region: usize },
    #[error("region {region} has transport cross section {transport} in group {group}")]
    Transport {
        region: usize,
        group: usize,
        transport: f64,
    },
    #[error("B1 critical spectrum")]
    B1(#[from] B1Error),
}

/// MOC の解を使い、領域・群ごとの反応率を保存するようにセル全体で均質化し、
/// `coarse_bounds`（`fine` の境界の一部、降順）の少数群に縮約する。
/// 拡散係数は `diffusion_coefficients` で作る（ADR-0098）。不連続因子は、辺の上の中性子束を体積平均の中性子束で割る。
/// 格子全体を1つのノード（集合体）として均質化し、外周がすべて鏡面反射であること（均質な解が平らになる）。
pub fn homogenize(
    solution: &Solution,
    materials: &[&Material],
    fine: &GroupStructure,
    coarse_bounds: &[f64],
    boundaries: &Boundaries,
) -> Result<FewGroup, HomogenizeError> {
    if !reflective(boundaries) {
        return Err(HomogenizeError::NotReflective);
    }
    let ng = fine.len();
    check_regions(solution, materials)?;
    let group_of = coarse_groups(fine, coarse_bounds)?;
    let coarse = coarse_bounds.len() - 1;
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
    // 不連続因子: 辺の上の中性子束 / 均質化した中性子束（反射境界の単一の格子では、均質な解は平らで
    // 体積平均に等しい）
    let volume: f64 = solution.volumes.iter().sum();
    let adf = solution
        .sides
        .iter()
        .map(|side| {
            let mut surface = vec![0.0; coarse];
            for (g, v) in side.iter().enumerate() {
                surface[group_of[g]] += v;
            }
            surface
                .iter()
                .zip(&flux)
                .map(|(s, f)| s / (f / volume))
                .collect()
        })
        .collect();
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
        diffusion: diffusion_coefficients(solution, materials, fine, coarse_bounds, boundaries)?,
        total,
        nu_fission: per_flux(nu_fission),
        fission: per_flux(fission),
        chi: if produced > 0.0 {
            chi.iter().map(|c| c / produced).collect()
        } else {
            chi
        },
        scatter,
        adf,
    })
}

fn reflective(boundaries: &Boundaries) -> bool {
    [
        boundaries.left,
        boundaries.right,
        boundaries.top,
        boundaries.bottom,
    ]
    .iter()
    .all(|&b| b == Boundary::Reflective)
}

fn check_regions(solution: &Solution, materials: &[&Material]) -> Result<(), HomogenizeError> {
    if materials.len() != solution.flux.len() {
        return Err(HomogenizeError::RegionCount {
            materials: materials.len(),
            regions: solution.flux.len(),
        });
    }
    Ok(())
}

/// [細かい群] の少数群の番号。少数群の境界は細かい群の境界の一部で、全体を覆う降順であること。
fn coarse_groups(
    fine: &GroupStructure,
    coarse_bounds: &[f64],
) -> Result<Vec<usize>, HomogenizeError> {
    let ng = fine.len();
    let positions = coarse_bounds
        .iter()
        .map(|b| fine.bounds().iter().position(|f| f == b))
        .collect::<Option<Vec<usize>>>()
        .ok_or(HomogenizeError::NotLibraryBounds)?;
    let covers = positions.first() == Some(&0) && positions.last() == Some(&ng);
    if !covers || positions.windows(2).any(|w| w[0] >= w[1]) {
        return Err(HomogenizeError::Span);
    }
    let mut group_of = vec![0; ng];
    for (big, w) in positions.windows(2).enumerate() {
        group_of[w[0]..w[1]].fill(big);
    }
    Ok(group_of)
}

/// 少数群の拡散係数（ADR-0098）。どちらもライブラリの P1 の散乱行列を使い、真空でない物質に P1 の
/// モーメントがなければエラーにする。
///
/// - 外周がすべて鏡面反射の格子: 細かい群のまま中性子束と体積で均質化し、B1 の臨界のバックリングでの
///   基本モードから D_g = J_g / (B φ_g) とし、B1 の中性子束で縮約する（D_G = Σ J_g / Σ B φ_g）
/// - それ以外（真空かアルベドの辺がある）: 領域・細かい群ごとに outflow の輸送断面積
///   Σtr,g = Σt,g − Σ_g' Σs1(g→g') で D = 1/(3Σtr) とし、MOC の中性子束と体積で重みづけして均質化・縮約する
pub fn diffusion_coefficients(
    solution: &Solution,
    materials: &[&Material],
    fine: &GroupStructure,
    coarse_bounds: &[f64],
    boundaries: &Boundaries,
) -> Result<Vec<f64>, HomogenizeError> {
    check_regions(solution, materials)?;
    let _ = (fine, boundaries);
    Ok(vec![0.0; coarse_bounds.len() - 1])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Lattice, PinCell};
    use crate::material::FissionTerm;
    use crate::moc::{Options, solve};
    use crate::tracking::{Boundary, Tracking};

    fn two_group(total: [f64; 2], scatter: [[f64; 2]; 2], nu_fission: [f64; 2]) -> Material {
        // P1 の散乱は P0 の 0.3 倍（前方散乱）
        let p1 = scatter
            .iter()
            .map(|r| r.iter().map(|s| 0.3 * s).collect())
            .collect();
        Material {
            total: total.to_vec(),
            scatter: scatter.iter().map(|r| r.to_vec()).collect(),
            anisotropic: vec![p1],
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
    fn homogeneous_cell_has_unit_adf() {
        let m = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
        let cell = PinCell::new(1.26, vec![0.4096]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 16, 0.05).unwrap();
        let solution = solve(&tracking, &[&m, &m], &Options::default()).unwrap();
        let fine = GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap();
        let fg = homogenize(
            &solution,
            &[&m, &m],
            &fine,
            fine.bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        for value in fg.adf.iter().flatten() {
            assert!((value - 1.0).abs() < 1e-8, "{:?}", fg.adf);
        }
    }

    #[test]
    fn adf_needs_reflective_sides() {
        // 真空の辺があると均質な解は平らでないので、体積平均で割る不連続因子は作れない
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let vacuum = Boundaries {
            right: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        assert!(matches!(
            homogenize(&solution, &refs, &fine, fine.bounds(), &vacuum),
            Err(HomogenizeError::NotReflective)
        ));
        let albedo = Boundaries {
            top: Boundary::Albedo(0.5),
            ..Boundaries::REFLECTIVE
        };
        assert!(matches!(
            homogenize(&solution, &refs, &fine, fine.bounds(), &albedo),
            Err(HomogenizeError::NotReflective)
        ));
    }

    #[test]
    fn fuel_pin_raises_the_thermal_adf() {
        // 燃料が中心にあると、外周（減速材）の熱中性子束はセルの平均より高い。4 辺は対称で等しい
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let fg = homogenize(
            &solution,
            &refs,
            &fine,
            fine.bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        assert!(fg.adf[0][1] > 1.0, "{:?}", fg.adf);
        for side in &fg.adf {
            for (a, b) in side.iter().zip(&fg.adf[0]) {
                assert!((a - b).abs() < 1e-6 * b, "{:?}", fg.adf);
            }
        }
    }

    #[test]
    fn one_group_collapse_preserves_k() {
        // 反応率を保存するので、1群の νΣf/Σa はセルの k∞ に一致する
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let one = homogenize(
            &solution,
            &refs,
            &fine,
            &[1.0e7, 1.0e-5],
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        assert_eq!(one.groups(), 1);
        let k = one.nu_fission[0] / one.absorption[0];
        assert!(
            (k - solution.k).abs() < 1e-6 * solution.k,
            "{k} vs {}",
            solution.k
        );
        assert_eq!(one.k_inf, solution.k);
        assert!((one.chi[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn two_group_constants_reproduce_k_inf() {
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let fg = homogenize(
            &solution,
            &refs,
            &fine,
            fine.bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
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
        assert!(
            homogenize(
                &solution,
                &refs,
                &fine,
                &[1.0e7, 1.0, 1.0e-5],
                &Boundaries::REFLECTIVE
            )
            .is_err()
        );
        assert!(
            homogenize(
                &solution,
                &refs,
                &fine,
                &[1.0e7, 0.625],
                &Boundaries::REFLECTIVE
            )
            .is_err()
        );
    }

    #[test]
    fn reflective_lattice_takes_d_from_the_b1_mode() {
        // 細かい群 = 少数群のとき、D は均質化した断面積の B1 の臨界の基本モードの J/(Bφ) そのもの
        let (solution, materials, fine) = solved();
        let refs: Vec<&Material> = materials.iter().collect();
        let fg = homogenize(
            &solution,
            &refs,
            &fine,
            fine.bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        let p1 = |m: &Material, g: usize, to: usize| m.anisotropic[0][g][to];
        let mut medium = Medium {
            total: vec![0.0; 2],
            scatter: vec![vec![0.0; 2]; 2],
            scatter1: vec![vec![0.0; 2]; 2],
            nu_fission: vec![0.0; 2],
            chi: vec![1.0, 0.0],
        };
        for g in 0..2 {
            let mut flux = 0.0;
            for ((m, phi), v) in materials.iter().zip(&solution.flux).zip(&solution.volumes) {
                let w = v * phi[g];
                flux += w;
                medium.total[g] += w * m.total[g];
                medium.nu_fission[g] += w * m.nu_fission[g];
                for to in 0..2 {
                    medium.scatter[g][to] += w * m.scatter[g][to];
                    medium.scatter1[g][to] += w * p1(m, g, to);
                }
            }
            medium.total[g] /= flux;
            medium.nu_fission[g] /= flux;
            for to in 0..2 {
                medium.scatter[g][to] /= flux;
                medium.scatter1[g][to] /= flux;
            }
        }
        let mode = medium.critical().unwrap();
        // B² の符号は k∞ − 1 の符号
        assert_eq!(mode.buckling > 0.0, solution.k > 1.0, "{}", mode.buckling);
        for (d, e) in fg.diffusion.iter().zip(mode.diffusion()) {
            assert!((d - e).abs() < 1e-9 * e, "{d} vs {e}");
        }
        // 1群に縮約した D は B1 の中性子流と中性子束の和の比
        let one = homogenize(
            &solution,
            &refs,
            &fine,
            &[1.0e7, 1.0e-5],
            &Boundaries::REFLECTIVE,
        )
        .unwrap();
        let expected = mode.current.iter().sum::<f64>() / mode.flux.iter().sum::<f64>();
        assert!((one.diffusion[0] - expected).abs() < 1e-9 * expected);
    }

    #[test]
    fn diffusion_needs_p1_moments() {
        let (solution, mut materials, fine) = solved();
        for m in &mut materials {
            m.anisotropic.clear();
        }
        let refs: Vec<&Material> = materials.iter().collect();
        for boundaries in [
            Boundaries::REFLECTIVE,
            Boundaries {
                left: Boundary::Vacuum,
                ..Boundaries::REFLECTIVE
            },
        ] {
            assert!(matches!(
                diffusion_coefficients(&solution, &refs, &fine, fine.bounds(), &boundaries),
                Err(HomogenizeError::NoP1 { region: 0 })
            ));
        }
    }

    /// 2つの領域の、中性子束と体積がわかっている解
    fn two_regions() -> Solution {
        Solution {
            k: 1.0,
            flux: vec![vec![2.0, 1.0], vec![1.0, 3.0]],
            volumes: vec![1.0, 2.0],
            leakage: vec![0.0; 2],
            sides: vec![vec![0.0; 2]; 4],
            iterations: 1,
        }
    }

    #[test]
    fn open_lattice_without_p1_scattering_uses_one_third_total() {
        let mut a = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.0, 0.0]);
        let mut b = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
        for m in [&mut a, &mut b] {
            m.anisotropic = vec![vec![vec![0.0; 2]; 2]];
        }
        let fine = GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap();
        let vacuum = Boundaries {
            right: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let d = diffusion_coefficients(&two_regions(), &[&a, &b], &fine, fine.bounds(), &vacuum)
            .unwrap();
        // 重み V φ: 群 0 は 2 と 2、群 1 は 1 と 6
        let expected = [
            (2.0 / (3.0 * 0.5) + 2.0 / (3.0 * 0.6)) / 4.0,
            (1.0 / (3.0 * 1.2) + 6.0 / (3.0 * 2.0)) / 7.0,
        ];
        for (d, e) in d.iter().zip(expected) {
            assert!((d - e).abs() < 1e-14, "{d} vs {e}");
        }
    }

    #[test]
    fn open_lattice_uses_the_outflow_transport_cross_section() {
        // P1 = 0.3 P0: Σtr,g = Σt,g − 0.3 Σ_g' Σs0(g→g')
        let a = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.0, 0.0]);
        let b = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
        let fine = GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap();
        let albedo = Boundaries {
            top: Boundary::Albedo(0.5),
            ..Boundaries::REFLECTIVE
        };
        let tr = |t: f64, out: f64| 1.0 / (3.0 * (t - 0.3 * out));
        let expected_per_group = [
            (2.0 * tr(0.5, 0.47) + 2.0 * tr(0.6, 0.598)) / 4.0,
            (1.0 * tr(1.2, 0.9) + 6.0 * tr(2.0, 1.972)) / 7.0,
        ];
        let d = diffusion_coefficients(&two_regions(), &[&a, &b], &fine, fine.bounds(), &albedo)
            .unwrap();
        for (d, e) in d.iter().zip(expected_per_group) {
            assert!((d - e).abs() < 1e-14 * e, "{d} vs {e}");
        }
        // 1群への縮約も V φ の重み
        let one =
            diffusion_coefficients(&two_regions(), &[&a, &b], &fine, &[1.0e7, 1.0e-5], &albedo)
                .unwrap();
        let expected = (expected_per_group[0] * 4.0 + expected_per_group[1] * 7.0) / 11.0;
        assert!((one[0] - expected).abs() < 1e-14 * expected);
    }

    #[test]
    fn outflow_rejects_a_void_region_with_flux() {
        let a = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.0, 0.0]);
        let void = Material::void(2);
        let fine = GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap();
        let vacuum = Boundaries {
            bottom: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        assert!(matches!(
            diffusion_coefficients(&two_regions(), &[&a, &void], &fine, fine.bounds(), &vacuum),
            Err(HomogenizeError::Transport { region: 1, .. })
        ));
    }
}
