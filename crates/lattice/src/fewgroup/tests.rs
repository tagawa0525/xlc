use super::*;
use crate::geometry::{Cell, Lattice, PinCell};
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
fn reflective_adfs_divide_by_the_average_flux() {
    // 外周がすべて鏡面反射なら、1ノードの問題の均質な解は平らで、ADF は辺の上の中性子束 / 体積平均
    // の中性子束（ADR-0093 までの定義）
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
    let volume: f64 = solution.volumes.iter().sum();
    for g in 0..2 {
        let average: f64 = solution
            .flux
            .iter()
            .zip(&solution.volumes)
            .map(|(f, v)| v * f[g])
            .sum::<f64>()
            / volume;
        for (side, values) in fg.adf.iter().take(4).enumerate() {
            let expected = solution.sides[side][g] / average;
            assert!(
                (values[g] - expected).abs() < 1e-6 * expected,
                "side {side} group {g}: {} vs {expected}",
                values[g]
            );
        }
        // 2次元の格子の z の辺は 1（ADR-0128）
        assert_eq!((fg.adf[4][g], fg.adf[5][g]), (1.0, 1.0));
    }
}

/// 燃料と水の 2×2 の格子（右が真空）の解
fn open_lattice() -> (Solution, Vec<Material>, GroupStructure, Boundaries) {
    let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
    let water = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let cell = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::new(vec![pin, cell], vec![vec![0, 1], vec![1, 0]]).unwrap();
    let boundaries = vacuum();
    let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, boundaries).unwrap();
    let materials = vec![
        fuel.clone(),
        water.clone(),
        water.clone(),
        water.clone(),
        fuel,
        water,
    ];
    let refs: Vec<&Material> = materials.iter().collect();
    let solution = solve(&tracking, &refs, &Options::default()).unwrap();
    (solution, materials, fine(), boundaries)
}

#[test]
fn open_lattice_adfs_divide_by_the_one_node_flux() {
    // 真空の辺がある格子の ADF は、参照の k、ノード平均の中性子束、辺の正味の流れで解いた1ノードの
    // NEM の均質な面の中性子束で割る（ADR-0103）。辺は（左、右、上、下）で、-x が左、+x が右、
    // -y が下、+y が上
    use xlc_numeric::nem::{Edge, OneNode, surface_fluxes};
    let (solution, materials, fine, boundaries) = open_lattice();
    let refs: Vec<&Material> = materials.iter().collect();
    let fg = homogenize(&solution, &refs, &fine, fine.bounds(), &boundaries).unwrap();
    let volume: f64 = solution.volumes.iter().sum();
    let average: Vec<f64> = (0..2)
        .map(|g| {
            solution
                .flux
                .iter()
                .zip(&solution.volumes)
                .map(|(f, v)| v * f[g])
                .sum::<f64>()
                / volume
        })
        .collect();
    let c = &solution.currents;
    let surface = surface_fluxes(&OneNode {
        diffusion: &fg.diffusion,
        total: &fg.total,
        scatter: &fg.scatter,
        nu_fission: &fg.nu_fission,
        chi: &fg.chi,
        k: solution.k,
        flux: &average,
        widths: solution.size,
        currents: [&c[0], &c[1], &c[3], &c[2]],
        edges: [
            Edge::Reflective,
            Edge::Vacuum,
            Edge::Reflective,
            Edge::Reflective,
        ],
    })
    .unwrap();
    for (side, face) in [(0, 0), (1, 1), (2, 3), (3, 2)] {
        for (g, homogeneous) in surface[face].iter().enumerate() {
            let expected = solution.sides[side][g] / homogeneous;
            let adf = fg.adf[side][g];
            assert!(
                (adf - expected).abs() < 1e-12 * expected,
                "side {side} group {g}: {adf} vs {expected}"
            );
        }
    }
    // 真空の辺では、均質な解が平らでないので、体積平均で割った値と違う
    for (g, mean) in average.iter().enumerate() {
        let flat = solution.sides[1][g] / mean;
        assert!(
            (fg.adf[1][g] - flat).abs() > 1e-2 * flat,
            "{} vs {flat}",
            fg.adf[1][g]
        );
    }
}

#[test]
fn vacuum_sides_carry_the_current_ratio() {
    // 真空の辺は、MOC の辺の外向きの正味の流れ / 辺の上の不均質な中性子束を、少数群ごとに細かい群の
    // 和の比で持つ（ADR-0101）。鏡面反射の辺は持たない
    let (solution, materials, fine, boundaries) = open_lattice();
    let refs: Vec<&Material> = materials.iter().collect();
    let fg = homogenize(&solution, &refs, &fine, fine.bounds(), &boundaries).unwrap();
    for side in [0, 2, 3] {
        assert_eq!(fg.vacuum_ratio[side], None, "side {side}");
    }
    let ratio = fg.vacuum_ratio[1].clone().unwrap();
    for (g, r) in ratio.iter().enumerate() {
        let expected = solution.currents[1][g] / solution.sides[1][g];
        assert!((r - expected).abs() < 1e-14 * expected, "{r} vs {expected}");
        // 出ていく角度中性子束だけなので P1 の 1/2 より大きく、1 より小さい
        assert!(*r > 0.5 && *r < 1.0, "{r}");
    }
    let one = homogenize(&solution, &refs, &fine, &[1.0e7, 1.0e-5], &boundaries).unwrap();
    let collapsed = one.vacuum_ratio[1].clone().unwrap();
    let expected = solution.currents[1].iter().sum::<f64>() / solution.sides[1].iter().sum::<f64>();
    assert!((collapsed[0] - expected).abs() < 1e-14 * expected);
    // 外周がすべて鏡面反射の格子は、どの辺も持たない
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
    assert_eq!(fg.vacuum_ratio, vec![None; 6]);
}

#[test]
fn adf_rejects_albedo_sides() {
    // core にアルベドの境界がないので、アルベドの辺の格子では ADF を作らない（ADR-0103）
    let (solution, materials, fine) = solved();
    let refs: Vec<&Material> = materials.iter().collect();
    let albedo = Boundaries {
        top: Boundary::Albedo(0.5),
        ..Boundaries::REFLECTIVE
    };
    assert!(matches!(
        homogenize(&solution, &refs, &fine, fine.bounds(), &albedo),
        Err(HomogenizeError::Albedo)
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
    for side in &fg.adf[..4] {
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
        currents: vec![vec![0.0; 2]; 4],
        face_fluxes: vec![vec![0.0; 2]; 4],
        face_currents: vec![vec![0.0; 2]; 4],
        // [領域][x, y][群]。体積 (1, 2) を掛けて ê_G へ射影した重みは、群 0 で (0.5, 2)、群 1 で (1.0, −0.2)
        // （群 1 は領域 1 で逆向き）。どちらの群も ê = (0.6, 0.8)
        region_currents: vec![
            [vec![0.3, 0.6], vec![0.4, 0.8]],
            [vec![0.6, -0.06], vec![0.8, -0.08]],
        ],
        size: [1.0, 3.0],
        iterations: 1,
    }
}

fn fine() -> GroupStructure {
    GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]).unwrap()
}

/// どの領域も媒質でない
fn no_media(solution: &Solution) -> Vec<Option<usize>> {
    vec![None; solution.volumes.len()]
}

fn vacuum() -> Boundaries {
    Boundaries {
        right: Boundary::Vacuum,
        ..Boundaries::REFLECTIVE
    }
}

#[test]
fn open_lattice_without_p1_scattering_uses_one_third_total() {
    // P1 が 0 なら Σtr = Σt。重みは中性子流（`two_regions` のコメント）。D = 1/(3⟨Σt⟩)
    let mut a = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.0, 0.0]);
    let mut b = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
    for m in [&mut a, &mut b] {
        m.anisotropic = vec![vec![vec![0.0; 2]; 2]];
    }
    let d = diffusion_coefficients(
        &two_regions(),
        &[&a, &b],
        &fine(),
        fine().bounds(),
        &vacuum(),
    )
    .unwrap();
    let expected = [
        1.0 / (3.0 * (0.5 * 0.5 + 2.0 * 0.6) / 2.5),
        1.0 / (3.0 * (1.0 * 1.2 - 0.2 * 2.0) / 0.8),
    ];
    for (d, e) in d.iter().zip(expected) {
        assert!((d - e).abs() < 1e-14 * e, "{d} vs {e}");
    }
}

#[test]
fn open_lattice_weights_the_outflow_transport_cross_section_by_currents() {
    // P1 = 0.3 P0: Σtr,g = Σt,g − 0.3 Σ_g' Σs0(g→g')。⟨Σtr⟩ を、ノードの正味の中性子流の向きへ射影した
    // 領域の中性子流 × 体積で重みづけし（符号つき）、D = 1/(3⟨Σtr⟩)（ADR-0104）
    let a = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.0, 0.0]);
    let b = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
    let albedo = Boundaries {
        top: Boundary::Albedo(0.5),
        ..Boundaries::REFLECTIVE
    };
    let tr = |t: f64, out: f64| t - 0.3 * out;
    let sums = [
        0.5 * tr(0.5, 0.47) + 2.0 * tr(0.6, 0.598),
        1.0 * tr(1.2, 0.9) - 0.2 * tr(2.0, 1.972),
    ];
    let expected = [1.0 / (3.0 * sums[0] / 2.5), 1.0 / (3.0 * sums[1] / 0.8)];
    let d = diffusion_coefficients(&two_regions(), &[&a, &b], &fine(), fine().bounds(), &albedo)
        .unwrap();
    for (d, e) in d.iter().zip(expected) {
        assert!((d - e).abs() < 1e-14 * e, "{d} vs {e}");
    }
    // 1群への縮約も同じ向き（ê = (0.6, 0.8)）への射影で重みづけする
    let one = diffusion_coefficients(
        &two_regions(),
        &[&a, &b],
        &fine(),
        &[1.0e7, 1.0e-5],
        &albedo,
    )
    .unwrap();
    let expected = 1.0 / (3.0 * (sums[0] + sums[1]) / 3.3);
    assert!(
        (one[0] - expected).abs() < 1e-14 * expected,
        "{one:?} vs {expected}"
    );
}

#[test]
fn current_weighting_counts_a_void_gap_with_its_current() {
    // 真空の隙間（ピンのギャップ）は Σtr = 0 で、中性子流は重みに入る
    let a = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.0, 0.0]);
    let void = Material::void(2);
    let d = diffusion_coefficients(
        &two_regions(),
        &[&a, &void],
        &fine(),
        fine().bounds(),
        &vacuum(),
    )
    .unwrap();
    let expected = [
        1.0 / (3.0 * 0.5 * tr_of(&a, 0) / 2.5),
        1.0 / (3.0 * 1.0 * tr_of(&a, 1) / 0.8),
    ];
    for (d, e) in d.iter().zip(expected) {
        assert!((d - e).abs() < 1e-14 * e, "{d} vs {e}");
    }
}

#[test]
fn current_weighting_needs_a_net_current_and_a_positive_mean() {
    let a = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.0, 0.0]);
    let b = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
    // ノードの正味の中性子流が 0（対称な位置）なら、向きが定まらない
    let mut zero = two_regions();
    zero.region_currents = vec![[vec![0.0; 2], vec![0.0; 2]]; 2];
    assert!(matches!(
        diffusion_coefficients(&zero, &[&a, &b], &fine(), fine().bounds(), &vacuum()),
        Err(HomogenizeError::NoNetCurrent { group: 0 })
    ));
    // 逆向きの重みが大きく、⟨Σtr⟩ が正でない
    let mut opposite = two_regions();
    // 群 0 の重みは (1.0, −0.9)（和 0.1）で、⟨Σtr⟩ = (0.359 − 0.9 × 0.4206)/0.1 < 0
    opposite.region_currents = vec![
        [vec![1.0, 0.6], vec![0.0, 0.8]],
        [vec![-0.45, -0.06], vec![0.0, -0.08]],
    ];
    assert!(matches!(
        diffusion_coefficients(&opposite, &[&a, &b], &fine(), fine().bounds(), &vacuum()),
        Err(HomogenizeError::Transport { group: 0, .. })
    ));
}

#[test]
fn open_lattice_d_uses_the_moc_region_currents() {
    // 格子の解の領域の中性子流で、手で重みづけした値と一致する
    let (solution, materials, fine, boundaries) = open_lattice();
    let refs: Vec<&Material> = materials.iter().collect();
    let fg = homogenize(&solution, &refs, &fine, fine.bounds(), &boundaries).unwrap();
    for g in 0..2 {
        let net: [f64; 2] = std::array::from_fn(|c| {
            solution
                .region_currents
                .iter()
                .zip(&solution.volumes)
                .map(|(j, v)| v * j[c][g])
                .sum()
        });
        let length = net[0].hypot(net[1]);
        let (mut weight, mut sum) = (0.0, 0.0);
        for ((j, v), m) in solution
            .region_currents
            .iter()
            .zip(&solution.volumes)
            .zip(&materials)
        {
            let w = v * (j[0][g] * net[0] + j[1][g] * net[1]) / length;
            weight += w;
            sum += w * tr_of(m, g);
        }
        let expected = 1.0 / (3.0 * sum / weight);
        assert!(
            (fg.diffusion[g] - expected).abs() < 1e-12 * expected,
            "group {g}"
        );
    }
}

fn tr_of(m: &Material, g: usize) -> f64 {
    m.total[g] - m.anisotropic[0][g].iter().sum::<f64>()
}

#[test]
fn reflective_lattice_without_fission_uses_the_p1_current_of_the_moc_flux() {
    // 核分裂がないと臨界のバックリングは定まらない。B² = 0 の中性子流の式
    // 3 ⟨Σt⟩_g d_g − 3 Σ_g' ⟨Σs1⟩(g'→g) d_g' = Φ_g を MOC の中性子束 Φ で解き、D_G = Σ d / Σ Φ
    let a = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.0, 0.0]);
    let b = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
    let flux = [4.0, 7.0];
    let hom = |g: usize, f: &dyn Fn(&Material) -> f64| {
        let w = [[2.0, 1.0], [2.0, 6.0]];
        (w[0][g] * f(&a) + w[1][g] * f(&b)) / flux[g]
    };
    let t0 = hom(0, &|m| m.total[0]);
    let t1 = hom(1, &|m| m.total[1]);
    let s = |from: usize, to: usize| hom(from, &|m| m.anisotropic[0][from][to]);
    // [3(t0 − s00), −3 s10; −3 s01, 3(t1 − s11)] d = Φ
    let m = [
        [3.0 * (t0 - s(0, 0)), -3.0 * s(1, 0)],
        [-3.0 * s(0, 1), 3.0 * (t1 - s(1, 1))],
    ];
    let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    let d = [
        (m[1][1] * flux[0] - m[0][1] * flux[1]) / det,
        (m[0][0] * flux[1] - m[1][0] * flux[0]) / det,
    ];
    let got = diffusion_coefficients(
        &two_regions(),
        &[&a, &b],
        &fine(),
        fine().bounds(),
        &Boundaries::REFLECTIVE,
    )
    .unwrap();
    for g in 0..2 {
        let e = d[g] / flux[g];
        assert!((got[g] - e).abs() < 1e-13 * e, "{} vs {e}", got[g]);
    }
    let one = diffusion_coefficients(
        &two_regions(),
        &[&a, &b],
        &fine(),
        &[1.0e7, 1.0e-5],
        &Boundaries::REFLECTIVE,
    )
    .unwrap();
    let e = (d[0] + d[1]) / (flux[0] + flux[1]);
    assert!((one[0] - e).abs() < 1e-13 * e);
}

/// 2つの核分裂性の物質（χ が違い、熱群から高速群への上方散乱がある）
fn fissile(
    total: [f64; 2],
    scatter: [[f64; 2]; 2],
    p1: [[f64; 2]; 2],
    nu_fission: [f64; 2],
    chi: [f64; 2],
) -> Material {
    Material {
        total: total.to_vec(),
        scatter: scatter.iter().map(|r| r.to_vec()).collect(),
        anisotropic: vec![p1.iter().map(|r| r.to_vec()).collect()],
        fission: nu_fission.iter().map(|v| v / 2.5).collect(),
        nu_fission: nu_fission.to_vec(),
        fission_terms: vec![FissionTerm {
            chi: chi.to_vec(),
            nu_fission: nu_fission.to_vec(),
        }],
    }
}

/// 2群の B1 の基本モード（手で書いた 2 × 2 の式）。(k, φ, d)
fn oracle_mode(
    total: [f64; 2],
    s0: [[f64; 2]; 2],
    s1: [[f64; 2]; 2],
    nu_fission: [f64; 2],
    chi: [f64; 2],
    b2: f64,
) -> (f64, [f64; 2], [f64; 2]) {
    let gamma = |t: f64| {
        let x = b2.abs().sqrt() / t;
        if b2 == 0.0 {
            1.0
        } else if b2 > 0.0 {
            x * x * x.atan() / (3.0 * (x - x.atan()))
        } else {
            x * x * x.atanh() / (3.0 * (x.atanh() - x))
        }
    };
    // d = M⁻¹ φ、M = [3γ0 t0 − 3 s1(0→0), −3 s1(1→0); −3 s1(0→1), 3γ1 t1 − 3 s1(1→1)]
    let m = [
        [
            3.0 * gamma(total[0]) * total[0] - 3.0 * s1[0][0],
            -3.0 * s1[1][0],
        ],
        [
            -3.0 * s1[0][1],
            3.0 * gamma(total[1]) * total[1] - 3.0 * s1[1][1],
        ],
    ];
    let det_m = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    let inv = [
        [m[1][1] / det_m, -m[0][1] / det_m],
        [-m[1][0] / det_m, m[0][0] / det_m],
    ];
    // A = diag(t) − S0ᵀ + B² M⁻¹、A φ = χ
    let a = [
        [
            total[0] - s0[0][0] + b2 * inv[0][0],
            -s0[1][0] + b2 * inv[0][1],
        ],
        [
            -s0[0][1] + b2 * inv[1][0],
            total[1] - s0[1][1] + b2 * inv[1][1],
        ],
    ];
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    let phi = [
        (a[1][1] * chi[0] - a[0][1] * chi[1]) / det,
        (a[0][0] * chi[1] - a[1][0] * chi[0]) / det,
    ];
    let d = [
        inv[0][0] * phi[0] + inv[0][1] * phi[1],
        inv[1][0] * phi[0] + inv[1][1] * phi[1],
    ];
    (nu_fission[0] * phi[0] + nu_fission[1] * phi[1], phi, d)
}

#[test]
fn b1_diffusion_matches_a_hand_computed_two_region_lattice() {
    let a = fissile(
        [0.5, 1.2],
        [[0.42, 0.05], [0.01, 0.95]],
        [[0.15, 0.01], [0.002, 0.3]],
        [0.016, 0.45],
        [1.0, 0.0],
    );
    let b = fissile(
        [0.6, 1.6],
        [[0.5, 0.08], [0.02, 1.4]],
        [[0.2, 0.02], [0.004, 0.5]],
        [0.008, 0.19],
        [0.6, 0.4],
    );
    let solution = Solution {
        k: 1.0,
        flux: vec![vec![2.0, 1.0], vec![1.5, 2.5]],
        volumes: vec![1.0, 0.5],
        leakage: vec![0.0; 2],
        sides: vec![vec![1.0; 2]; 4],
        currents: vec![vec![0.0; 2]; 4],
        face_fluxes: vec![vec![1.0; 2]; 4],
        face_currents: vec![vec![0.0; 2]; 4],
        region_currents: vec![[vec![0.0; 2], vec![0.0; 2]]; 2],
        size: [1.0, 0.5],
        iterations: 1,
    };
    // 重み V φ: A は (2, 1)、B は (0.75, 1.25)。中性子束 (2.75, 2.25)
    let w = [[2.0, 1.0], [0.75, 1.25]];
    let flux = [2.75, 2.25];
    let hom = |g: usize, xa: f64, xb: f64| (w[0][g] * xa + w[1][g] * xb) / flux[g];
    let total = [hom(0, 0.5, 0.6), hom(1, 1.2, 1.6)];
    let s0 = [
        [hom(0, 0.42, 0.5), hom(0, 0.05, 0.08)],
        [hom(1, 0.01, 0.02), hom(1, 0.95, 1.4)],
    ];
    let s1 = [
        [hom(0, 0.15, 0.2), hom(0, 0.01, 0.02)],
        [hom(1, 0.002, 0.004), hom(1, 0.3, 0.5)],
    ];
    let nu_fission = [hom(0, 0.016, 0.008), hom(1, 0.45, 0.19)];
    let rate_a = 2.0 * 0.016 + 1.0 * 0.45;
    let rate_b = 0.75 * 0.008 + 1.25 * 0.19;
    let chi = [
        (rate_a + 0.6 * rate_b) / (rate_a + rate_b),
        0.4 * rate_b / (rate_a + rate_b),
    ];
    let k = |b2: f64| oracle_mode(total, s0, s1, nu_fission, chi, b2);
    // 二分法（k は B² について単調に減る）
    let k_inf = k(0.0).0;
    assert!(k_inf > 1.0, "{k_inf}");
    let (mut low, mut high) = (0.0, 1.0);
    assert!(k(high).0 < 1.0);
    for _ in 0..200 {
        let middle = 0.5 * (low + high);
        if k(middle).0 > 1.0 {
            low = middle;
        } else {
            high = middle;
        }
    }
    let (_, phi, d) = k(0.5 * (low + high));
    let fg = homogenize(
        &solution,
        &[&a, &b],
        &fine(),
        fine().bounds(),
        &Boundaries::REFLECTIVE,
    )
    .unwrap();
    for g in 0..2 {
        let e = d[g] / phi[g];
        assert!(
            (fg.diffusion[g] - e).abs() < 1e-8 * e,
            "{} vs {e}",
            fg.diffusion[g]
        );
    }
    // 1群: B1 の中性子束で縮約する。MOC の中性子束で重みづけした値とは違う
    let one = homogenize(
        &solution,
        &[&a, &b],
        &fine(),
        &[1.0e7, 1.0e-5],
        &Boundaries::REFLECTIVE,
    )
    .unwrap();
    let e = (d[0] + d[1]) / (phi[0] + phi[1]);
    assert!(
        (one.diffusion[0] - e).abs() < 1e-8 * e,
        "{} vs {e}",
        one.diffusion[0]
    );
    let moc = (d[0] / phi[0] * flux[0] + d[1] / phi[1] * flux[1]) / (flux[0] + flux[1]);
    assert!((moc - e).abs() > 1e-3 * e, "{moc} vs {e}");
}

/// 燃料のピンと水のセルを横に並べた帯（1行2列、右が真空）と、その解
fn strip() -> (Lattice, Solution, Vec<Material>, Boundaries) {
    let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
    let water = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let cell = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::rectilinear(
        vec![1.26; 2],
        vec![1.26],
        vec![Cell::Pin(pin), Cell::Pin(cell)],
        vec![vec![0, 1]],
    )
    .unwrap();
    let boundaries = vacuum();
    let tracking = Tracking::with_boundaries(&lattice, 32, 0.02, boundaries).unwrap();
    let materials = vec![fuel, water.clone(), water];
    let refs: Vec<&Material> = materials.iter().collect();
    let solution = solve(&tracking, &refs, &Options::default()).unwrap();
    (lattice, solution, materials, boundaries)
}

fn node(columns: std::ops::Range<usize>, rows: std::ops::Range<usize>) -> Node {
    Node {
        columns,
        rows,
        reflector: false,
    }
}

#[test]
fn nodes_of_an_infinite_medium_are_the_medium() {
    // 1つの物質の格子（外周がすべて鏡面反射）を2つのノードに分けると、どちらも物質の断面積そのもので、
    // 面の正味の流れは 0、ADF は 1。複数のノードの格子の拡散係数は outflow（ADR-0098）で、
    // 1/(3(Σt − Σs1))
    let m = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
    let cell = PinCell::new(1.26, vec![0.4]).unwrap();
    let lattice = Lattice::new(vec![cell], vec![vec![0, 0], vec![0, 0]]).unwrap();
    let tracking = Tracking::new(&lattice, 16, 0.05).unwrap();
    let refs = vec![&m; 8];
    let solution = solve(&tracking, &refs, &Options::default()).unwrap();
    let fine = fine();
    let nodes = [node(0..1, 0..2), node(1..2, 0..2)];
    let fgs = homogenize_nodes(
        &solution,
        &lattice,
        &refs,
        &fine,
        fine.bounds(),
        &Boundaries::REFLECTIVE,
        &nodes,
        &no_media(&solution),
        false,
    )
    .unwrap();
    assert_eq!(fgs.len(), 2);
    for fg in &fgs {
        for g in 0..2 {
            assert!((fg.total[g] - m.total[g]).abs() < 1e-12, "{:?}", fg.total);
            assert!((fg.nu_fission[g] - m.nu_fission[g]).abs() < 1e-12);
            let transport = m.total[g] - m.anisotropic[0][g].iter().sum::<f64>();
            let d = 1.0 / (3.0 * transport);
            assert!(
                (fg.diffusion[g] - d).abs() < 1e-12 * d,
                "{:?}",
                fg.diffusion
            );
        }
        for value in fg.adf.iter().flatten() {
            assert!((value - 1.0).abs() < 1e-6, "{:?}", fg.adf);
        }
        assert_eq!(fg.vacuum_ratio, vec![None; 6]);
        assert_eq!(fg.normalized_adf, vec![false; 6]);
    }
}

#[test]
fn nodes_of_one_medium_take_its_b1_diffusion() {
    // 均質な燃料のセル 2 つと水のセル（右が真空）。燃料のセルの領域が媒質なら、燃料のノードの D は媒質そのものの
    // B1 の D（ADR-0141）で、ノードの中性子流の重みに依らない。媒質でなければ中性子流で重みづけした outflow の D
    let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
    let water = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
    let cell = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::rectilinear(
        vec![1.26; 3],
        vec![1.26],
        vec![Cell::Pin(cell)],
        vec![vec![0, 0, 0]],
    )
    .unwrap();
    let boundaries = vacuum();
    let tracking = Tracking::with_boundaries(&lattice, 32, 0.02, boundaries).unwrap();
    let refs = vec![&fuel, &fuel, &water];
    let solution = solve(&tracking, &refs, &Options::default()).unwrap();
    let fine = fine();
    let nodes = [
        node(0..1, 0..1),
        node(1..2, 0..1),
        Node {
            reflector: true,
            ..node(2..3, 0..1)
        },
    ];
    let solve_nodes = |media: &[Option<usize>]| {
        homogenize_nodes(
            &solution,
            &lattice,
            &refs,
            &fine,
            fine.bounds(),
            &boundaries,
            &nodes,
            media,
            false,
        )
        .unwrap()
    };
    let group_of = coarse_groups(&fine, fine.bounds()).unwrap();
    let b1 = b1_diffusion(&solution, &[&fuel], &[0], &group_of, 2).unwrap();
    let with_media = solve_nodes(&[Some(0), Some(0), None]);
    let without = solve_nodes(&no_media(&solution));
    for n in 0..2 {
        for (g, &b) in b1.iter().enumerate() {
            let d = with_media[n].diffusion[g];
            assert!((d - b).abs() < 1e-12 * b, "node {n}: {d} vs {b}");
            let other = without[n].diffusion[g];
            assert!((other - b).abs() > 1e-3 * b, "node {n}: {other}");
        }
    }
    // 媒質と媒質でない物質が混ざったノードは、中性子流の重みのまま
    let mixed = homogenize_nodes(
        &solution,
        &lattice,
        &refs,
        &fine,
        fine.bounds(),
        &boundaries,
        &[
            node(0..2, 0..1),
            Node {
                reflector: true,
                ..node(2..3, 0..1)
            },
        ],
        &[Some(0), None, None],
        false,
    )
    .unwrap();
    let whole = homogenize_nodes(
        &solution,
        &lattice,
        &refs,
        &fine,
        fine.bounds(),
        &boundaries,
        &[
            node(0..2, 0..1),
            Node {
                reflector: true,
                ..node(2..3, 0..1)
            },
        ],
        &no_media(&solution),
        false,
    )
    .unwrap();
    assert_eq!(mixed[0].diffusion, whole[0].diffusion);
    // 核分裂のない媒質のノードは、媒質の印に依らず ADR-0104 と ADR-0134 の D（Koebke の D を中性子流の重みの D から
    // 探す）
    let water_medium = solve_nodes(&[None, None, Some(0)]);
    assert_eq!(water_medium[2].diffusion, without[2].diffusion);
}

#[test]
fn strip_nodes_make_get_constants_and_normalized_adfs() {
    // 燃料と水の帯を、燃料のノードと反射体のノードに分ける（ADR-0095）
    use xlc_numeric::nem::{Edge, OneNode, surface_fluxes};
    let (lattice, solution, materials, boundaries) = strip();
    let refs: Vec<&Material> = materials.iter().collect();
    let fine = fine();
    let nodes = [
        node(0..1, 0..1),
        Node {
            reflector: true,
            ..node(1..2, 0..1)
        },
    ];
    let fgs = homogenize_nodes(
        &solution,
        &lattice,
        &refs,
        &fine,
        fine.bounds(),
        &boundaries,
        &nodes,
        &no_media(&solution),
        false,
    )
    .unwrap();
    let (fuel, reflector) = (&fgs[0], &fgs[1]);
    // 定数は、ノードの領域だけの反応率の保存（燃料は領域 0, 1、反射体は領域 2）
    let average = |regions: std::ops::Range<usize>, g: usize, x: &dyn Fn(&Material) -> f64| {
        let (mut rate, mut flux) = (0.0, 0.0);
        for r in regions {
            let w = solution.volumes[r] * solution.flux[r][g];
            rate += w * x(&materials[r]);
            flux += w;
        }
        (rate / flux, flux / 1.26_f64.powi(2))
    };
    let mut fluxes = [[0.0; 2]; 2];
    for (node_fluxes, (fg, regions)) in fluxes.iter_mut().zip([(fuel, 0..2), (reflector, 2..3)]) {
        for (g, node_flux) in node_fluxes.iter_mut().enumerate() {
            let (total, flux) = average(regions.clone(), g, &|m| m.total[g]);
            assert!((fg.total[g] - total).abs() < 1e-12 * total);
            // 燃料のノードの D は、領域の中性子流をノードの正味の中性子流の向きへ射影した重みの outflow の
            // Σtr から（ADR-0104）。反射体のノードの D は等価な D（下で x の面の不連続因子がそろうことを見る）
            if regions.start == 0 {
                let net: [f64; 2] = std::array::from_fn(|c| {
                    regions
                        .clone()
                        .map(|r| solution.volumes[r] * solution.region_currents[r][c][g])
                        .sum()
                });
                let length = net[0].hypot(net[1]);
                let (mut weight, mut sum) = (0.0, 0.0);
                for r in regions.clone() {
                    let j = &solution.region_currents[r];
                    let w = solution.volumes[r] * (j[0][g] * net[0] + j[1][g] * net[1]) / length;
                    let m = &materials[r];
                    weight += w;
                    sum += w * (m.total[g] - m.anisotropic[0][g].iter().sum::<f64>());
                }
                let d = 1.0 / (3.0 * sum / weight);
                assert!(
                    (fg.diffusion[g] - d).abs() < 1e-12 * d,
                    "{} vs {d}",
                    fg.diffusion[g]
                );
            }
            *node_flux = flux;
        }
        assert_eq!(fg.k_inf, solution.k);
    }
    // 面の値: 縦の面は 0（左）、1（燃料と反射体の境）、2（右）、横の面は 3, 4（上）と 5, 6（下）
    let het = |face: usize| solution.face_fluxes[face].clone();
    let out = |face: usize, sign: f64| -> Vec<f64> {
        solution.face_currents[face]
            .iter()
            .map(|j| sign * j)
            .collect()
    };
    let one_node =
        |fg: &FewGroup, flux: &[f64], faces: [usize; 4], signs: [f64; 4], edges: [Edge; 4]| {
            let c: Vec<Vec<f64>> = (0..4).map(|i| out(faces[i], signs[i])).collect();
            surface_fluxes(&OneNode {
                diffusion: &fg.diffusion,
                total: &fg.total,
                scatter: &fg.scatter,
                nu_fission: &fg.nu_fission,
                chi: &fg.chi,
                k: solution.k,
                flux,
                widths: [1.26, 1.26],
                currents: [&c[0], &c[1], &c[2], &c[3]],
                edges,
            })
            .unwrap()
        };
    // 隣の横方向（y）の漏れは、鏡面反射の上下の面の正味の流れ（収束の許容差の程度で 0）
    let y_leakage = |top: usize, bottom: usize| -> Vec<f64> {
        (0..2)
            .map(|g| (-solution.face_currents[top][g] + solution.face_currents[bottom][g]) / 1.26)
            .collect()
    };
    let (fuel_y, reflector_y) = (y_leakage(3, 5), y_leakage(4, 6));
    let fuel_surface = one_node(
        fuel,
        &fluxes[0],
        [0, 1, 5, 3],
        [-1.0, 1.0, 1.0, -1.0],
        [
            Edge::Reflective,
            Edge::Neighbor {
                width: 1.26,
                leakage: &reflector_y,
            },
            Edge::Reflective,
            Edge::Reflective,
        ],
    );
    let reflector_surface = one_node(
        reflector,
        &fluxes[1],
        [1, 2, 6, 4],
        [-1.0, 1.0, 1.0, -1.0],
        [
            Edge::Neighbor {
                width: 1.26,
                leakage: &fuel_y,
            },
            Edge::Vacuum,
            Edge::Reflective,
            Edge::Reflective,
        ],
    );
    for g in 0..2 {
        // 燃料の右の辺の ADF は、境の面の不均質な中性子束 / 均質な面の中性子束
        let expected = het(1)[g] / fuel_surface[1][g];
        assert!((fuel.adf[1][g] - expected).abs() < 1e-10 * expected);
        // 反射体の左の辺は f* = f_反射体 / f_燃料
        let f_reflector = het(1)[g] / reflector_surface[0][g];
        let expected = f_reflector / fuel.adf[1][g];
        assert!(
            (reflector.adf[0][g] - expected).abs() < 1e-10 * expected,
            "{} vs {expected}",
            reflector.adf[0][g]
        );
        // 真空の辺の ADF と J/φ
        let expected = het(2)[g] / reflector_surface[1][g];
        assert!((reflector.adf[1][g] - expected).abs() < 1e-10 * expected);
        // 反射体のノードの x の2つの面の不連続因子は等しい（Koebke の等価な D。ADR-0104）
        let left = reflector.adf[0][g] * fuel.adf[1][g];
        assert!(
            (left - reflector.adf[1][g]).abs() < 1e-9 * left,
            "group {g}: {left} vs {}",
            reflector.adf[1][g]
        );
        let ratio = reflector.vacuum_ratio[1].as_ref().unwrap()[g];
        let expected = solution.face_currents[2][g] / het(2)[g];
        assert!((ratio - expected).abs() < 1e-12 * expected);
    }
    assert_eq!(fuel.normalized_adf, vec![false; 6]);
    assert_eq!(
        reflector.normalized_adf,
        vec![true, false, false, false, false, false]
    );
    assert_eq!(fuel.vacuum_ratio, vec![None; 6]);
    assert_eq!(
        reflector
            .vacuum_ratio
            .iter()
            .map(Option::is_some)
            .collect::<Vec<_>>(),
        vec![false, true, false, false, false, false]
    );
}

#[test]
fn reflector_nodes_away_from_fuel_keep_plain_factors() {
    // 燃料と接しない反射体のノード（軸方向のノズル、炉心板、水の層）は、f* でない GET の ADF を持つ（ADR-0128）。
    // 燃料と接する反射体のノードだけが、燃料の側の辺に f* の印を持つ
    let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
    let water = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let cell = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::rectilinear(
        vec![1.26; 3],
        vec![1.26],
        vec![Cell::Pin(pin), Cell::Pin(cell)],
        vec![vec![0, 1, 1]],
    )
    .unwrap();
    let boundaries = vacuum();
    let tracking = Tracking::with_boundaries(&lattice, 32, 0.02, boundaries).unwrap();
    let materials = [fuel, water.clone(), water.clone(), water];
    let refs: Vec<&Material> = materials.iter().collect();
    let solution = solve(&tracking, &refs, &Options::default()).unwrap();
    let fine = fine();
    let run = |reflector: bool| {
        let nodes: Vec<Node> = [0..1, 1..2, 2..3]
            .into_iter()
            .enumerate()
            .map(|(n, columns)| Node {
                reflector: reflector && n > 0,
                ..node(columns, 0..1)
            })
            .collect();
        homogenize_nodes(
            &solution,
            &lattice,
            &refs,
            &fine,
            fine.bounds(),
            &boundaries,
            &nodes,
            &no_media(&solution),
            false,
        )
        .unwrap()
    };
    let (marked, plain) = (run(true), run(false));
    assert_eq!(marked[1].normalized_adf[..4], [true, false, false, false]);
    assert_eq!(marked[2].normalized_adf, vec![false; 6]);
    assert_eq!(marked[2].adf, plain[2].adf);
    for g in 0..2 {
        let back = marked[1].adf[0][g] * plain[0].adf[1][g];
        assert!((back - plain[1].adf[0][g]).abs() < 1e-14 * back);
    }
}

#[test]
fn invalid_nodes_are_errors() {
    let (lattice, solution, materials, boundaries) = strip();
    let refs: Vec<&Material> = materials.iter().collect();
    let fine = fine();
    let run = |nodes: &[Node], boundaries: &Boundaries| {
        homogenize_nodes(
            &solution,
            &lattice,
            &refs,
            &fine,
            fine.bounds(),
            boundaries,
            nodes,
            &no_media(&solution),
            false,
        )
    };
    let cases: Vec<(Vec<Node>, &str)> = vec![
        (vec![], "no nodes"),
        (vec![node(0..1, 0..1)], "not covered"),
        (vec![node(0..2, 0..1), node(1..2, 0..1)], "overlap"),
        (vec![node(0..1, 0..1), node(1..3, 0..1)], "outside"),
        (vec![node(0..1, 0..1), node(1..1, 0..1)], "empty"),
    ];
    for (nodes, what) in &cases {
        assert!(
            matches!(run(nodes, &boundaries), Err(HomogenizeError::Nodes(_))),
            "{what}"
        );
    }
    let albedo = Boundaries {
        top: Boundary::Albedo(0.5),
        ..boundaries
    };
    assert!(matches!(
        run(&[node(0..1, 0..1), node(1..2, 0..1)], &albedo),
        Err(HomogenizeError::Albedo)
    ));
    // 隣が1つのノードで面がそろわない（左のノードの右の辺に2つのノードが接する）
    let m = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
    let cell = PinCell::new(1.26, vec![0.4]).unwrap();
    let square = Lattice::new(vec![cell], vec![vec![0, 0], vec![0, 0]]).unwrap();
    let tracking = Tracking::new(&square, 16, 0.05).unwrap();
    let refs = vec![&m; 8];
    let solution = solve(&tracking, &refs, &Options::default()).unwrap();
    let misaligned = homogenize_nodes(
        &solution,
        &square,
        &refs,
        &fine,
        fine.bounds(),
        &Boundaries::REFLECTIVE,
        &[node(0..1, 0..2), node(1..2, 0..1), node(1..2, 1..2)],
        &no_media(&solution),
        false,
    );
    assert!(
        matches!(misaligned, Err(HomogenizeError::Nodes(_))),
        "{misaligned:?}"
    );
}

/// `nodes`（x が速く、-y の行から並べた nx × ny の格子）の定数と形状（ADR-0131）で、core の形の角の点の釣り合い
/// （未知数が角の点の不均質な中性子束、均質な角の値は Φ / CDF、真空の辺の流れは r × ADF × 均質な角の値）を解くと、格子の角の点の不均質な中性子束（角の点に
/// 接する領域の中性子束の平均）に戻り、再構成したセルの出力（Σf × 均質な中性子束のセルの平均 × 形状関数）が
/// MOC のセルの核分裂率に戻ることを確かめる。ノード平均と面の正味の流れは、core の NEM が GET で再現する参照の値
fn assert_pins_reproduce(
    lattice: &Lattice,
    solution: &Solution,
    materials: &[Material],
    boundaries: &Boundaries,
    nodes: &[Node],
    nx: usize,
    ny: usize,
) {
    use xlc_numeric::pins::FACE_POINTS;
    use xlc_numeric::pins::corners::{CornerNode, corner_balance};
    use xlc_numeric::pins::node::{Conditions, Constants, NodeModel};
    let refs: Vec<&Material> = materials.iter().collect();
    let fine = fine();
    let fgs = homogenize_nodes(
        solution,
        lattice,
        &refs,
        &fine,
        fine.bounds(),
        boundaries,
        nodes,
        &no_media(solution),
        true,
    )
    .unwrap();
    let references = node_references(
        solution,
        lattice,
        &refs,
        &fine,
        fine.bounds(),
        boundaries,
        nodes,
        &no_media(solution),
    )
    .unwrap()
    .references;
    let rule = xlc_numeric::quadrature::gauss_legendre(FACE_POINTS).unwrap();
    let cells_rule = xlc_numeric::quadrature::gauss_legendre(6).unwrap();
    let pins: Vec<&Pins> = fgs.iter().map(|f| f.pins.as_ref().unwrap()).collect();
    // ノード n の角 c の角の点（x が速い (nx + 1) × (ny + 1)）
    let point = |n: usize, c: usize| {
        let (i, j) = (n % nx, n / nx);
        let (ci, cj) = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)][c];
        ci + (nx + 1) * cj
    };
    // 炉心の CDF: 燃料の CDF で割った角は、角の点に接する燃料（核分裂のある）のノードの CDF の調和平均を掛ける
    let refs_fissile = |f: &FewGroup| f.nu_fission.iter().any(|&v| v > 0.0);
    let cdf = |n: usize, c: usize, g: usize| -> f64 {
        if !pins[n].normalized_corners[c] {
            return pins[n].corners[c][g];
        }
        let p = point(n, c);
        let fuel: Vec<f64> = (0..nodes.len())
            .filter(|&m| refs_fissile(&fgs[m]))
            .flat_map(|m| {
                (0..4)
                    .filter(move |&d| point(m, d) == p)
                    .map(move |d| (m, d))
            })
            .map(|(m, d)| pins[m].corners[d][g])
            .collect();
        assert!(!fuel.is_empty(), "node {n} corner {c}");
        let harmonic = fuel.len() as f64 / fuel.iter().map(|f| 1.0 / f).sum::<f64>();
        pins[n].corners[c][g] * harmonic
    };
    let models: Vec<NodeModel> = references
        .iter()
        .zip(&fgs)
        .map(|(r, f)| {
            let constants = Constants {
                total: &f.total,
                scatter: &f.scatter,
                nu_fission: &f.nu_fission,
                chi: &f.chi,
                diffusion: &f.diffusion,
            };
            NodeModel::new(
                &constants,
                solution.k,
                [0.5 * r.size[0], 0.5 * r.size[1]],
                &rule,
            )
            .unwrap()
        })
        .collect();
    let conditions: Vec<Conditions> = references
        .iter()
        .map(|r| {
            // 辺（左、右、上、下）の外向きの流れを、面（-x、+x、-y、+y）の座標の正の向きの流れにする
            let current =
                |side: usize, sign: f64| [sign * r.currents[side][0], sign * r.currents[side][1]];
            Conditions {
                half: [0.5 * r.size[0], 0.5 * r.size[1]],
                average: [r.average[0], r.average[1]],
                corners: [[0.0; 2]; 4],
                currents: [
                    current(0, -1.0),
                    current(1, 1.0),
                    current(3, -1.0),
                    current(2, 1.0),
                ],
                axial_leakage: [0.0; 2],
            }
        })
        .collect();
    let corner_nodes: Vec<Option<CornerNode>> = (0..nodes.len())
        .map(|n| {
            Some(CornerNode {
                model: &models[n],
                conditions: &conditions[n],
                diffusion: [fgs[n].diffusion[0], fgs[n].diffusion[1]],
                weights: std::array::from_fn(|c| std::array::from_fn(|g| 1.0 / cdf(n, c, g))),
                outflow: std::array::from_fn(|face| {
                    let side = [0, 1, 3, 2][face];
                    fgs[n].vacuum_ratio[side]
                        .as_ref()
                        .map(|r| std::array::from_fn(|g| r[g] * fgs[n].adf[side][g]))
                }),
            })
        })
        .collect();
    let zero = vec![false; (nx + 1) * (ny + 1)];
    let phi = corner_balance(nx, ny, &corner_nodes, &zero).unwrap();
    // 角の点の不均質な中性子束: 角の点に接する格子の中の領域の中性子束の平均（少数群は 2 群のまま）
    let (width, height): (f64, f64) = (
        lattice.widths().iter().sum(),
        lattice.heights().iter().sum(),
    );
    let x_edges: Vec<f64> = (0..=nx)
        .map(|i| nodes[..i].iter().map(|m| r_width(lattice, m)).sum())
        .collect();
    let y_edges: Vec<f64> = (0..=ny)
        .map(|j| (0..j).map(|jj| r_height(lattice, &nodes[jj * nx])).sum())
        .collect();
    for cj in 0..=ny {
        for ci in 0..=nx {
            let (x, y) = (x_edges[ci], y_edges[cj]);
            let e = 1e-9 * width.max(height);
            let mut want = [0.0; 2];
            let mut count = 0.0;
            for (sx, sy) in [(1.0, 1.0), (-1.0, 1.0), (-1.0, -1.0), (1.0, -1.0)] {
                let (px, py) = (x + sx * e, y + sy * e);
                if px <= 0.0 || px >= width || py <= 0.0 || py >= height {
                    continue;
                }
                let region = lattice.region_at(px - 0.5 * width, py - 0.5 * height);
                for (w, f) in want.iter_mut().zip(&solution.flux[region]) {
                    *w += f;
                }
                count += 1.0;
            }
            let got = phi[ci + (nx + 1) * cj];
            for g in 0..2 {
                let want = want[g] / count;
                assert!(
                    (got[g] - want).abs() < 1e-9 * want,
                    "corner point ({ci}, {cj}) group {g}: {} vs {want}",
                    got[g]
                );
            }
        }
    }
    // セルの出力が MOC のセルの核分裂率（解析的な体積で平均）に戻る
    let volumes = lattice.volumes();
    for (n, node) in nodes.iter().enumerate() {
        let mut c = conditions[n].clone();
        c.corners = std::array::from_fn(|k| {
            let p = phi[point(n, k)];
            std::array::from_fn(|g| p[g] / cdf(n, k, g))
        });
        let flux = models[n].solve(&c);
        let [a, b] = c.half;
        let mut y0 = -b;
        for (j, dy) in pins[n].dy.iter().enumerate() {
            let mut x0 = -a;
            for (i, dx) in pins[n].dx.iter().enumerate() {
                let average = flux.average(x0, x0 + dx, y0, y0 + dy, &cells_rule);
                // 行は -y から（格子の行は上から）
                let row = node.rows.end - 1 - j;
                let cell = row * lattice.columns() + node.columns.start + i;
                let (mut rate, mut area) = ([0.0; 2], 0.0);
                for region in lattice.cell_regions(cell) {
                    area += volumes[region];
                    for (g, r) in rate.iter_mut().enumerate() {
                        *r += volumes[region] * solution.flux[region][g] * refs[region].fission[g];
                    }
                }
                for g in 0..2 {
                    let want = rate[g] / area;
                    let got = fgs[n].fission[g] * average[g] * pins[n].form[g][j][i];
                    assert!(
                        (got - want).abs() < 1e-9 * want.max(1e-12),
                        "node {n} cell ({j}, {i}) group {g}: {got} vs {want}"
                    );
                }
                x0 += dx;
            }
            y0 += dy;
        }
    }
}

/// ノードの幅と高さ
fn r_width(lattice: &Lattice, node: &Node) -> f64 {
    lattice.widths()[node.columns.clone()].iter().sum()
}

fn r_height(lattice: &Lattice, node: &Node) -> f64 {
    lattice.heights()[node.rows.clone()].iter().sum()
}

#[test]
fn strip_pins_reproduce_the_lattice() {
    // 燃料と水の帯（右が真空）を燃料と反射体のノードに分けた定数と形状は、定数を作った体系の角の中性子束と
    // セルの核分裂率を再現する（ADR-0131）
    let (lattice, solution, materials, boundaries) = strip();
    let nodes = [
        node(0..1, 0..1),
        Node {
            reflector: true,
            ..node(1..2, 0..1)
        },
    ];
    assert_pins_reproduce(&lattice, &solution, &materials, &boundaries, &nodes, 2, 1);
}

/// 燃料のピン 2 × 2 の右と上に水のセルを 1 列ずつ並べた 3 × 3 と、外周の境界 `boundaries` の解。燃料は右と上の
/// 辺に接しない
fn l_shaped(boundaries: Boundaries) -> (Lattice, Solution, Vec<Material>, Boundaries) {
    let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
    let water = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let cell = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::rectilinear(
        vec![1.26; 3],
        vec![1.26; 3],
        vec![Cell::Pin(pin), Cell::Pin(cell)],
        vec![vec![1, 1, 1], vec![0, 0, 1], vec![0, 0, 1]],
    )
    .unwrap();
    let tracking = Tracking::with_boundaries(&lattice, 32, 0.02, boundaries).unwrap();
    // ピンの中心の領域だけ燃料、ほかは水
    let materials: Vec<Material> = (0..lattice.cells())
        .flat_map(|cell| {
            let fuelled = lattice
                .cell_type(cell)
                .pin()
                .is_some_and(|p| !p.radii().is_empty());
            lattice
                .cell_regions(cell)
                .enumerate()
                .map(|(k, _)| {
                    if fuelled && k == 0 {
                        fuel.clone()
                    } else {
                        water.clone()
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect();
    let refs: Vec<&Material> = materials.iter().collect();
    let solution = solve(&tracking, &refs, &Options::default()).unwrap();
    (lattice, solution, materials, boundaries)
}

#[test]
fn l_shaped_pins_reproduce_the_lattice() {
    // 燃料のノード（2 × 2 のピン）の上と右に反射体のノード、右上に水のノード（反射体の印なし。燃料に角だけで
    // 接する）。反射体の角は燃料の面と、右上のノードに接する（ADR-0131）
    let (lattice, solution, materials, boundaries) = l_shaped(vacuum());
    let reflector = |n: Node| Node {
        reflector: true,
        ..n
    };
    // x が速く、-y の行（格子の行 1、2）から
    let nodes = [
        node(0..2, 1..3),
        reflector(node(2..3, 1..3)),
        reflector(node(0..2, 0..1)),
        node(2..3, 0..1),
    ];
    assert_pins_reproduce(&lattice, &solution, &materials, &boundaries, &nodes, 2, 2);
}

#[test]
fn two_vacuum_sides_pins_reproduce_the_lattice() {
    // 右と上の 2 辺が真空の格子でも、真空の辺の流れを均質な角の値に課す角の点の釣り合いで、定数と形状は定数を
    // 作った体系を再現する（ADR-0131。#0065）
    let (lattice, solution, materials, boundaries) = l_shaped(Boundaries {
        top: Boundary::Vacuum,
        ..vacuum()
    });
    let reflector = |n: Node| Node {
        reflector: true,
        ..n
    };
    let nodes = [
        node(0..2, 1..3),
        reflector(node(2..3, 1..3)),
        reflector(node(0..2, 0..1)),
        node(2..3, 0..1),
    ];
    assert_pins_reproduce(&lattice, &solution, &materials, &boundaries, &nodes, 2, 2);
}

#[test]
fn diagonal_reflector_nodes_normalize_only_the_fuel_corner() {
    // 角型の格子（#0029）: 燃料に角の点だけで接する右上のノードも反射体とする。反射体の f* の印は燃料と接する
    // 辺だけ（反射体どうしの辺は印なし）、C* の印は燃料のノードに接する角の点だけ（ADR-0095、ADR-0131）。
    // 定数と形状は、定数を作った体系を再現する
    let (lattice, solution, materials, boundaries) = l_shaped(Boundaries {
        top: Boundary::Vacuum,
        ..vacuum()
    });
    let reflector = |n: Node| Node {
        reflector: true,
        ..n
    };
    // x が速く、-y の行（格子の行 1、2）から
    let nodes = [
        node(0..2, 1..3),
        reflector(node(2..3, 1..3)),
        reflector(node(0..2, 0..1)),
        reflector(node(2..3, 0..1)),
    ];
    let refs: Vec<&Material> = materials.iter().collect();
    let fine = fine();
    let fgs = homogenize_nodes(
        &solution,
        &lattice,
        &refs,
        &fine,
        fine.bounds(),
        &boundaries,
        &nodes,
        &no_media(&solution),
        true,
    )
    .unwrap();
    // 辺は（左、右、上、下、−z、+z）、角は（-x -y、+x -y、+x +y、-x +y）
    assert_eq!(
        fgs[1].normalized_adf,
        vec![true, false, false, false, false, false]
    );
    assert_eq!(
        fgs[2].normalized_adf,
        vec![false, false, false, true, false, false]
    );
    assert_eq!(fgs[3].normalized_adf, vec![false; 6]);
    let corners = |n: usize| fgs[n].pins.as_ref().unwrap().normalized_corners.clone();
    assert_eq!(corners(1), vec![true, false, false, true]);
    assert_eq!(corners(2), vec![true, true, false, false]);
    assert_eq!(corners(3), vec![true, false, false, false]);
    assert_pins_reproduce(&lattice, &solution, &materials, &boundaries, &nodes, 2, 2);
}

#[test]
fn non_positive_corner_scales_are_errors() {
    // 角の点の大きさ s_p が正でない（0、負、NaN）と CDF を作れないのでエラー（ADR-0131）
    assert!(positive_scales(&[[1.0, 0.9], [1.1, 1.0]]).is_ok());
    for bad in [0.0, -0.5, f64::NAN] {
        let err = positive_scales(&[[1.0, 0.9], [1.1, bad]]).unwrap_err();
        assert!(err.contains("corner point 1"), "{err}");
    }
}

#[test]
fn pins_need_reflector_nodes_without_fission() {
    // 燃料の CDF で割る角は、核分裂のあるノードを燃料とする（ADR-0131）。反射体の印のノードに核分裂があると、
    // 自分の CDF で割ることになり、core で掛け戻す値と食い違うので、形状関数を作るときはエラー
    let (lattice, solution, materials, boundaries) = strip();
    let refs: Vec<&Material> = materials.iter().collect();
    let fine = fine();
    let nodes = [
        Node {
            reflector: true,
            ..node(0..1, 0..1)
        },
        node(1..2, 0..1),
    ];
    let err = homogenize_nodes(
        &solution,
        &lattice,
        &refs,
        &fine,
        fine.bounds(),
        &boundaries,
        &nodes,
        &no_media(&solution),
        true,
    )
    .unwrap_err();
    assert!(
        matches!(&err, HomogenizeError::Nodes(m) if m.contains("fission")),
        "{err:?}"
    );
}

#[test]
fn strip_corner_factors_follow_the_face_discontinuity() {
    // 帯（上と下は鏡面反射）では均質な解が x にしか変わらず、境の角の均質な値の見込みは境の面の均質な中性子束に
    // なる。反射体の燃料に接する角の C* は、境の面の f* に一致する（ADR-0131）
    let (lattice, solution, materials, boundaries) = strip();
    let refs: Vec<&Material> = materials.iter().collect();
    let fine = fine();
    let nodes = [
        node(0..1, 0..1),
        Node {
            reflector: true,
            ..node(1..2, 0..1)
        },
    ];
    let fgs = homogenize_nodes(
        &solution,
        &lattice,
        &refs,
        &fine,
        fine.bounds(),
        &boundaries,
        &nodes,
        &no_media(&solution),
        true,
    )
    .unwrap();
    let reflector = &fgs[1];
    let pins = reflector.pins.as_ref().unwrap();
    // 燃料に接する角は -x -y と -x +y
    assert_eq!(pins.normalized_corners, vec![true, false, false, true]);
    assert!(
        !fgs[0]
            .pins
            .as_ref()
            .unwrap()
            .normalized_corners
            .contains(&true)
    );
    for c in [0, 3] {
        for g in 0..2 {
            let (got, want) = (pins.corners[c][g], reflector.adf[0][g]);
            assert!(
                (got - want).abs() < 1e-6 * want,
                "corner {c} group {g}: {got} vs {want}"
            );
        }
    }
    // 反射体の形状関数は 0
    assert!(pins.form.iter().flatten().flatten().all(|&f| f == 0.0));
}

#[test]
fn one_reflective_node_keeps_the_assembly_pins() {
    // 外周がすべて鏡面反射の格子を1つのノードにすると、均質な解は平らで、CDF は角の不均質な中性子束 /
    // ノード平均（`pin_shapes`）に一致する。形状関数の分母は Σf × 均質な中性子束（トラックの体積で平均した
    // ノード平均）で、`pin_shapes` の分母（解析的な体積で平均した核分裂率）との比だけ違う
    let (lattice, solution, materials) = margined(0.04);
    let refs: Vec<&Material> = materials.iter().collect();
    let fine = fine();
    let assembly = pin_shapes(
        &solution,
        &lattice,
        &refs,
        &fine,
        fine.bounds(),
        &Boundaries::REFLECTIVE,
    )
    .unwrap();
    let fgs = homogenize_nodes(
        &solution,
        &lattice,
        &refs,
        &fine,
        fine.bounds(),
        &Boundaries::REFLECTIVE,
        &[node(0..4, 0..4)],
        &no_media(&solution),
        true,
    )
    .unwrap();
    let pins = fgs[0].pins.as_ref().unwrap();
    assert_eq!((&pins.dx, &pins.dy), (&assembly.dx, &assembly.dy));
    // 平らさは MOC の鏡面反射の辺の正味の流れ（収束の許容差の程度で 0）の程度
    let flat = |got: f64, want: f64| (got - want).abs() < 1e-6 * want;
    for (got, want) in pins
        .corners
        .iter()
        .flatten()
        .zip(assembly.corners.iter().flatten())
    {
        assert!(flat(*got, *want), "{got} vs {want}");
    }
    let volumes = lattice.volumes();
    let area: f64 = volumes.iter().sum();
    let (track_volume, track) = (solution.volumes.iter().sum::<f64>(), &solution.volumes);
    for g in 0..2 {
        let analytic = (0..refs.len())
            .map(|r| volumes[r] * solution.flux[r][g] * refs[r].fission[g])
            .sum::<f64>()
            / area;
        let average = (0..refs.len())
            .map(|r| track[r] * solution.flux[r][g])
            .sum::<f64>()
            / track_volume;
        let ratio = analytic / (fgs[0].fission[g] * average);
        for (got, want) in pins.form[g]
            .iter()
            .flatten()
            .zip(assembly.form[g].iter().flatten())
        {
            assert!(
                (got - want * ratio).abs() < 1e-6 * want.max(1e-3),
                "group {g}: {got} vs {want} × {ratio}"
            );
        }
    }
    assert_eq!(pins.normalized_corners, vec![false; 4]);
}

/// 2×2 のピン（対角に燃料）の周りに幅 0.04 の余白を付けた格子、領域ごとの物質、鏡面反射の解
fn margined(margin: f64) -> (Lattice, Solution, Vec<Material>) {
    let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
    let water = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let cell = PinCell::new(1.26, vec![]).unwrap();
    let lattice =
        Lattice::with_margin(vec![pin, cell], vec![vec![0, 1], vec![1, 0]], margin).unwrap();
    // ピンの中心の領域だけ燃料、ほかは水
    let materials: Vec<Material> = (0..lattice.cells())
        .flat_map(|cell| {
            let fuelled = lattice
                .cell_type(cell)
                .pin()
                .is_some_and(|p| !p.radii().is_empty());
            lattice
                .cell_regions(cell)
                .enumerate()
                .map(|(k, _)| {
                    if fuelled && k == 0 {
                        fuel.clone()
                    } else {
                        water.clone()
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect();
    let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
    let refs: Vec<&Material> = materials.iter().collect();
    let solution = solve(&tracking, &refs, &Options::default()).unwrap();
    (lattice, solution, materials)
}

#[test]
fn pin_shapes_average_to_one_and_take_corner_fluxes() {
    // 形状関数はセルの群の核分裂率の平均 / ノードの群の核分裂率の平均で、面積で重みを付けたノードの平均は 1。
    // 角の不連続因子は、余白の角のセルの中性子束 / ノード平均の中性子束（ADR-0094、ADR-0096）。形状関数は
    // 解析的な体積で、角の不連続因子はトラックの体積で平均する
    let (lattice, solution, materials) = margined(0.04);
    let refs: Vec<&Material> = materials.iter().collect();
    let pins = pin_shapes(
        &solution,
        &lattice,
        &refs,
        &fine(),
        fine().bounds(),
        &Boundaries::REFLECTIVE,
    )
    .unwrap();
    let close = |a: f64, b: f64| (a - b).abs() < 1e-12 * b.abs().max(1.0);
    // 対角の対称は、MOC の収束の程度で成り立つ
    let symmetric = |a: f64, b: f64| (a - b).abs() < 1e-6 * b.abs();
    // セルの大きさは格子の列の幅と行の高さ（行は下から）
    assert_eq!(pins.dx, vec![0.04, 1.26, 1.26, 0.04]);
    assert_eq!(pins.dy, vec![0.04, 1.26, 1.26, 0.04]);
    let area = |j: usize, i: usize| pins.dx[i] * pins.dy[j];
    for g in 0..2 {
        let form = &pins.form[g];
        assert_eq!((form.len(), form[0].len()), (4, 4));
        let mean: f64 = (0..4)
            .flat_map(|j| (0..4).map(move |i| (j, i)))
            .map(|(j, i)| area(j, i) * form[j][i])
            .sum::<f64>()
            / (2.6 * 2.6);
        assert!(close(mean, 1.0), "group {g}: {mean}");
        // 燃料は上の行（下から 2 行目）の左と、下の行（下から 1 行目）の右。対称なので等しい
        assert!(
            form[2][1] > 1.0 && symmetric(form[2][1], form[1][2]),
            "{form:?}"
        );
        for (j, i) in [
            (2, 2),
            (1, 1),
            (0, 0),
            (0, 3),
            (3, 0),
            (3, 3),
            (0, 1),
            (2, 0),
        ] {
            assert_eq!(form[j][i], 0.0, "cell ({j}, {i})");
        }
    }
    // 角のセル（格子の行は上から）: -x -y は左下、+x -y は右下、+x +y は右上、-x +y は左上
    let corner_cells = [(3, 0), (3, 3), (0, 3), (0, 0)];
    // 面の ADF と同じく、トラックから推定した体積で平均する
    let volumes = &solution.volumes;
    let volume: f64 = volumes.iter().sum();
    for g in 0..2 {
        let average = solution
            .flux
            .iter()
            .zip(volumes)
            .map(|(f, v)| v * f[g])
            .sum::<f64>()
            / volume;
        for (c, (r, col)) in corner_cells.into_iter().enumerate() {
            let cell = r * 4 + col;
            let regions = lattice.cell_regions(cell);
            let (flux, v) = regions.fold((0.0, 0.0), |(f, v), i| {
                (f + volumes[i] * solution.flux[i][g], v + volumes[i])
            });
            let expected = flux / v / average;
            assert!(close(pins.corners[c][g], expected), "corner {c} group {g}");
        }
        // 左下と右上は水のピンの角、左上と右下は燃料のピンの角で、対称
        assert!(symmetric(pins.corners[0][g], pins.corners[2][g]));
        assert!(symmetric(pins.corners[1][g], pins.corners[3][g]));
    }
}

#[test]
fn pin_shapes_need_a_reflective_lattice_with_corner_margins() {
    let (lattice, solution, materials) = margined(0.04);
    let refs: Vec<&Material> = materials.iter().collect();
    let shapes = |boundaries: &Boundaries| {
        pin_shapes(
            &solution,
            &lattice,
            &refs,
            &fine(),
            fine().bounds(),
            boundaries,
        )
    };
    // 格子全体の形状関数は外周がすべて鏡面反射の格子でだけ作る（ほかの格子は `nodes`。ADR-0131）
    let err = shapes(&vacuum()).unwrap_err().to_string();
    assert!(err.contains("reflective"), "{err}");
    // 角の中性子束は余白の角のセルで取る
    let (lattice, solution, materials) = margined(0.0);
    let refs: Vec<&Material> = materials.iter().collect();
    let err = pin_shapes(
        &solution,
        &lattice,
        &refs,
        &fine(),
        fine().bounds(),
        &Boundaries::REFLECTIVE,
    )
    .unwrap_err()
    .to_string();
    assert!(err.contains("margin"), "{err}");
}

/// `layout` のピン（0 が燃料、1 が水）の周りに幅 0.04 の余白を付けた格子、物質、鏡面反射の解と、集合体の定数
fn assembly(layout: Vec<Vec<usize>>) -> (Lattice, Solution, Vec<Material>, FewGroup) {
    let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
    let water = two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]);
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let cell = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::with_margin(vec![pin, cell], layout, 0.04).unwrap();
    let materials: Vec<Material> = (0..lattice.cells())
        .flat_map(|cell| {
            let fuelled = lattice
                .cell_type(cell)
                .pin()
                .is_some_and(|p| !p.radii().is_empty());
            lattice
                .cell_regions(cell)
                .enumerate()
                .map(|(k, _)| {
                    if fuelled && k == 0 {
                        fuel.clone()
                    } else {
                        water.clone()
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect();
    let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
    let refs: Vec<&Material> = materials.iter().collect();
    let solution = solve(&tracking, &refs, &Options::default()).unwrap();
    let mut few = homogenize(
        &solution,
        &refs,
        &fine(),
        fine().bounds(),
        &Boundaries::REFLECTIVE,
    )
    .unwrap();
    few.pins = Some(
        pin_shapes(
            &solution,
            &lattice,
            &refs,
            &fine(),
            fine().bounds(),
            &Boundaries::REFLECTIVE,
        )
        .unwrap(),
    );
    (lattice, solution, materials, few)
}

fn split_of(
    (lattice, solution, materials, few): &(Lattice, Solution, Vec<Material>, FewGroup),
    n: usize,
) -> Result<Vec<FewGroup>, HomogenizeError> {
    let refs: Vec<&Material> = materials.iter().collect();
    split_constants(
        few,
        solution,
        lattice,
        &refs,
        &fine(),
        fine().bounds(),
        &Boundaries::REFLECTIVE,
        n,
    )
}

#[test]
fn one_by_one_split_is_the_assembly() {
    // 分けなければ、集合体の定数そのもの（ADR-0108）
    let a = assembly(vec![vec![0, 1], vec![1, 0]]);
    let subs = split_of(&a, 1).unwrap();
    assert_eq!(subs.len(), 1);
    // 集合体の ADF は、MOC の外周の正味の流れ（鏡面反射なので収束の程度で 0）を条件にした1ノードの NEM で作る。
    // 分けた定数は均質な解を平らとするので、MOC の収束の程度で一致する
    let close = |x: f64, y: f64| (x - y).abs() < 1e-8 * y.abs();
    for (s, w) in subs[0].adf.iter().flatten().zip(a.3.adf.iter().flatten()) {
        assert!(close(*s, *w), "{:?} vs {:?}", subs[0].adf, a.3.adf);
    }
    assert_eq!(subs[0].pins, a.3.pins);
    assert_eq!(subs[0].total, a.3.total);
}

#[test]
fn two_by_two_split_takes_segment_fluxes() {
    // 2 × 2 のピン（左上と右下が燃料）を 2 × 2 に分ける。分ける線はピンのセルの境に乗る。外周の辺の ADF は、
    // その辺の半分の区間の面の中性子束 / ノード平均（トラックの体積）。内側の辺は 1。断面積は集合体のもの
    let a = assembly(vec![vec![0, 1], vec![1, 0]]);
    let (lattice, solution, _, few) = &a;
    let subs = split_of(&a, 2).unwrap();
    assert_eq!(subs.len(), 4);
    let volume: f64 = solution.volumes.iter().sum();
    let average: Vec<f64> = (0..2)
        .map(|g| {
            solution
                .flux
                .iter()
                .zip(&solution.volumes)
                .map(|(f, v)| v * f[g])
                .sum::<f64>()
                / volume
        })
        .collect();
    // 左の辺（縦の面の線 0）の下の区間は行 2、3（行は上から。行 3 は余白）
    let (rows, columns) = (lattice.rows(), lattice.columns());
    let heights = lattice.heights();
    let face = |row: usize| crate::cmfd::vertical_face(columns, row, 0);
    for (g, avg) in average.iter().enumerate() {
        let lower = (solution.face_fluxes[face(2)][g] * heights[2]
            + solution.face_fluxes[face(3)][g] * heights[3])
            / (heights[2] + heights[3]);
        let got = subs[0].adf[0][g];
        assert!((got - lower / avg).abs() < 1e-12 * got, "group {g}");
        // 内側の辺（部分ノード 0 の右と上）は 1
        assert_eq!((subs[0].adf[1][g], subs[0].adf[2][g]), (1.0, 1.0));
    }
    assert!(rows == 4 && columns == 4);
    for s in &subs {
        assert_eq!(
            (&s.total, &s.diffusion, s.k_inf),
            (&few.total, &few.diffusion, few.k_inf)
        );
    }
    // 部分ノード 0（左下）の形状関数は、集合体の左下の 2 × 2 のセル。分ける線はセルの境なので半分のセルはない
    let pins = subs[0].pins.as_ref().unwrap();
    let whole = few.pins.as_ref().unwrap();
    assert_eq!(pins.dx, vec![0.04, 1.26]);
    assert_eq!(pins.dy, vec![0.04, 1.26]);
    assert_eq!(pins.halves, vec![false; 4]);
    for g in 0..2 {
        assert_eq!(
            pins.form[g],
            vec![
                whole.form[g][0][..2].to_vec(),
                whole.form[g][1][..2].to_vec()
            ]
        );
    }
    // 集合体の外周の角（-x -y）は集合体の CDF（足し合わせの順の丸めまで）
    for (got, want) in pins.corners[0].iter().zip(&whole.corners[0]) {
        assert!((got - want).abs() < 1e-12 * want);
    }
    // 集合体の中心（部分ノード 0 の +x +y の角）は、4 つのピンのセルの角の領域の中性子束の平均
    let center: Vec<f64> = (0..2)
        .map(|g| {
            let e = 1e-7;
            [(e, e), (-e, e), (-e, -e), (e, -e)]
                .iter()
                .map(|&(x, y)| solution.flux[lattice.region_at(x, y)][g])
                .sum::<f64>()
                / 4.0
                / average[g]
        })
        .collect();
    for (got, want) in pins.corners[2].iter().zip(&center) {
        assert!((got - want).abs() < 1e-12 * want);
    }
}

#[test]
fn odd_lattices_split_through_the_middle_pins() {
    // 3 × 3 のピンを 2 × 2 に分けると、中央の行と列のピンのセルを半分ずつに分ける
    let a = assembly(vec![vec![0, 1, 0], vec![1, 0, 1], vec![0, 1, 0]]);
    let subs = split_of(&a, 2).unwrap();
    let pins = subs[0].pins.as_ref().unwrap();
    let close = |got: &[f64], want: &[f64]| {
        got.len() == want.len() && got.iter().zip(want).all(|(a, b)| (a - b).abs() < 1e-12)
    };
    assert!(close(&pins.dx, &[0.04, 1.26, 0.63]), "{:?}", pins.dx);
    assert!(close(&pins.dy, &[0.04, 1.26, 0.63]), "{:?}", pins.dy);
    assert_eq!(pins.halves, vec![false, true, true, false]);
    let whole = a.3.pins.as_ref().unwrap();
    // 半分のセルは元のセルの形状関数を持つ（中央のピンは行 2（下から）、列 2）
    assert_eq!(pins.form[0][2][2], whole.form[0][2][2]);
    // 右上の部分ノードは中央のピンの残りの半分から始まる
    let top_right = subs[3].pins.as_ref().unwrap();
    assert!(
        close(&top_right.dx, &[0.63, 1.26, 0.04]),
        "{:?}",
        top_right.dx
    );
    assert_eq!(top_right.halves, vec![true, false, false, true]);
}

#[test]
fn splits_need_a_reflective_lattice() {
    let a = assembly(vec![vec![0, 1], vec![1, 0]]);
    let refs: Vec<&Material> = a.2.iter().collect();
    let err = split_constants(
        &a.3,
        &a.1,
        &a.0,
        &refs,
        &fine(),
        fine().bounds(),
        &vacuum(),
        2,
    )
    .unwrap_err()
    .to_string();
    assert!(err.contains("reflective"), "{err}");
    assert!(split_of(&a, 0).is_err());
}
