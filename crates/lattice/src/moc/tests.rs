use super::*;
use crate::geometry::{Cell, Lattice, PinCell};
use crate::material::FissionTerm;
use crate::symmetry::{detect, require};
use crate::tracking::{Boundaries, Boundary};

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

fn fuel() -> Material {
    two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25])
}

fn moderator() -> Material {
    two_group([0.6, 2.0], [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0])
}

/// 無限媒質の k = νΣfᵀ (diag Σt - Σsᵀ)⁻¹ χ（χ = (1, 0)）
fn infinite_medium_k(m: &Material) -> f64 {
    let a = [
        [m.total[0] - m.scatter[0][0], -m.scatter[1][0]],
        [-m.scatter[0][1], m.total[1] - m.scatter[1][1]],
    ];
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    let phi = [a[1][1] / det, -a[1][0] / det];
    m.nu_fission[0] * phi[0] + m.nu_fission[1] * phi[1]
}

#[test]
fn sides_of_an_infinite_medium_see_the_mean_flux() {
    // 均質な無限媒質では中性子束は平らで等方なので、外周の辺の上の中性子束は領域の中性子束と同じ
    let cell = PinCell::new(1.26, vec![0.3, 0.45]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.1).unwrap();
    let m = fuel();
    let solution = solve(&tracking, &[&m, &m, &m], &Options::default()).unwrap();
    assert_eq!(solution.sides.len(), 4);
    for side in &solution.sides {
        for (s, f) in side.iter().zip(&solution.flux[0]) {
            assert!(
                (s - f).abs() < 1e-8 * f,
                "{side:?} vs {:?}",
                solution.flux[0]
            );
        }
    }
}

#[test]
fn faces_of_an_infinite_medium_see_the_mean_flux() {
    // 均質な無限媒質では、内側の面（セルの境）の上の中性子束も領域の中性子束と同じで、正味の中性子流は 0
    let cell = PinCell::new(1.26, vec![0.3, 0.45]).unwrap();
    let lattice = Lattice::new(vec![cell], vec![vec![0, 0], vec![0, 0]]).unwrap();
    let tracking = Tracking::new(&lattice, 8, 0.1).unwrap();
    let m = fuel();
    let solution = solve(&tracking, &[&m; 12], &Options::default()).unwrap();
    // 2×2 のセルの面は 縦 2×3 と横 3×2
    assert_eq!(solution.face_fluxes.len(), 12);
    assert_eq!(solution.face_currents.len(), 12);
    let mean = &solution.flux[0];
    for (face, values) in solution.face_fluxes.iter().enumerate() {
        for (g, (v, f)) in values.iter().zip(mean).enumerate() {
            assert!(
                (v - f).abs() < 1e-8 * f,
                "face {face} group {g}: {v} vs {f}"
            );
            let j = solution.face_currents[face][g];
            assert!(j.abs() < 1e-8 * f, "face {face} group {g}: current {j}");
        }
    }
}

#[test]
fn faces_agree_with_sides_and_mirror_lattices() {
    // 外周の面を辺ごとに長さで重みづけして平均すると、辺の値（`sides`、`currents`）に一致する。
    // 真ん中が対称面の格子の、真ん中の内側の面の中性子束は、半分の格子の鏡面反射の辺の値に近い
    // （方位角がセルの比で変わる cyclic tracking の差だけ違う）
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let open = Boundaries {
        left: Boundary::Vacuum,
        right: Boundary::Vacuum,
        ..Boundaries::REFLECTIVE
    };
    let full = Lattice::new(
        vec![pin.clone(), water.clone()],
        vec![
            vec![1, 0, 0, 1],
            vec![1, 0, 0, 1],
            vec![1, 0, 0, 1],
            vec![1, 0, 0, 1],
        ],
    )
    .unwrap();
    let (f, m) = (fuel(), moderator());
    // 行ごとの領域: 水、ピン（燃料、減速材）、ピン、水
    let row = [&m, &f, &m, &f, &m, &m];
    let materials: Vec<&Material> = (0..4).flat_map(|_| row).collect();
    let tracking = Tracking::with_boundaries(&full, 32, 0.02, open).unwrap();
    let s = solve(&tracking, &materials, &Options::default()).unwrap();
    // 縦の面 (r, i) は r × 5 + i、横の面 (j, c) は 20 + j × 4 + c
    let vertical = |r: usize, i: usize| r * 5 + i;
    let horizontal = |j: usize, c: usize| 20 + j * 4 + c;
    let sides: [Vec<(usize, f64)>; 4] = [
        (0..4).map(|r| (vertical(r, 0), -1.0)).collect(),
        (0..4).map(|r| (vertical(r, 4), 1.0)).collect(),
        (0..4).map(|c| (horizontal(0, c), -1.0)).collect(),
        (0..4).map(|c| (horizontal(4, c), 1.0)).collect(),
    ];
    for (side, faces) in sides.iter().enumerate() {
        for g in 0..2 {
            let flux: f64 = faces.iter().map(|&(i, _)| s.face_fluxes[i][g]).sum::<f64>() / 4.0;
            let out: f64 = faces
                .iter()
                .map(|&(i, o)| o * s.face_currents[i][g])
                .sum::<f64>()
                / 4.0;
            let (expected, current) = (s.sides[side][g], s.currents[side][g]);
            assert!(
                (flux - expected).abs() < 1e-12 * expected,
                "{flux} vs {expected}"
            );
            assert!(
                (out - current).abs() < 1e-12 * expected,
                "{out} vs {current}"
            );
        }
    }
    let half = Lattice::rectilinear(
        vec![1.26; 2],
        vec![1.26; 4],
        vec![Cell::Pin(pin), Cell::Pin(water)],
        vec![vec![1, 0]; 4],
    )
    .unwrap();
    let tracking = Tracking::with_boundaries(
        &half,
        32,
        0.02,
        Boundaries {
            left: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        },
    )
    .unwrap();
    let materials: Vec<&Material> = (0..4).flat_map(|_| [&m, &f, &m]).collect();
    let h = solve(&tracking, &materials, &Options::default()).unwrap();
    for g in 0..2 {
        let center: f64 = (0..4)
            .map(|r| s.face_fluxes[vertical(r, 2)][g])
            .sum::<f64>()
            / 4.0;
        // 全体の中性子束の規格化（生成率 1）は半分の格子の 2 倍の体積にわたるので、半分の値の 1/2
        let expected = 0.5 * h.sides[1][g];
        assert!(
            (center - expected).abs() < 2e-3 * expected,
            "group {g}: {center} vs {expected}"
        );
        let net: f64 = (0..4)
            .map(|r| s.face_currents[vertical(r, 2)][g])
            .sum::<f64>();
        assert!(net.abs() < 1e-6 * center, "{net}");
    }
}

#[test]
fn side_currents_add_up_to_the_leakage() {
    // 真空の辺（右、下）の外向きの正味の中性子流 × 辺の長さの和は漏れに等しく、鏡面反射の辺（左、上）
    // では 0
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap();
    let boundaries = Boundaries {
        right: Boundary::Vacuum,
        bottom: Boundary::Vacuum,
        ..Boundaries::REFLECTIVE
    };
    let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, boundaries).unwrap();
    let (f, m) = (fuel(), moderator());
    let s = solve(&tracking, &[&f, &m, &m, &m, &f, &m], &Options::default()).unwrap();
    assert_eq!(s.size, [2.52, 2.52]);
    assert_eq!(s.currents.len(), 4);
    for g in 0..2 {
        let leakage = s.leakage[g];
        assert!(leakage > 0.0);
        for side in [0, 2] {
            assert!(
                s.currents[side][g].abs() * 2.52 < 1e-6 * leakage,
                "{:?}",
                s.currents
            );
        }
        for side in [1, 3] {
            assert!(s.currents[side][g] > 0.0, "{:?}", s.currents);
        }
        let out = 2.52 * (s.currents[1][g] + s.currents[3][g]);
        assert!((out - leakage).abs() < 1e-6 * leakage, "{out} vs {leakage}");
    }
}

#[test]
fn homogeneous_cell_matches_infinite_medium() {
    let cell = PinCell::new(1.26, vec![0.3, 0.45]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 8, 0.1).unwrap();
    let m = fuel();
    let solution = solve(&tracking, &[&m, &m, &m], &Options::default()).unwrap();
    let expected = infinite_medium_k(&m);
    assert!(
        (solution.k - expected).abs() < 1e-6,
        "{} vs {expected}",
        solution.k
    );
}

fn heterogeneous(n_azimuthal: usize, spacing: f64) -> Solution {
    let cell = PinCell::new(1.26, vec![0.25, 0.4096, 0.55]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), n_azimuthal, spacing).unwrap();
    let (f, m) = (fuel(), moderator());
    solve(&tracking, &[&f, &f, &m, &m], &Options::default()).unwrap()
}

#[test]
fn heterogeneous_cell_balances_neutrons() {
    let s = heterogeneous(16, 0.05);
    let materials = [fuel(), fuel(), moderator(), moderator()];
    let (mut absorption, mut production) = (0.0, 0.0);
    for (i, m) in materials.iter().enumerate() {
        for g in 0..2 {
            absorption += s.volumes[i] * m.absorption()[g] * s.flux[i][g];
            production += s.volumes[i] * m.nu_fission[g] * s.flux[i][g];
        }
    }
    assert!((absorption - production / s.k).abs() < 1e-6 * absorption);
    // 熱中性子束は減速材で高く、燃料の中心で低い
    assert!(s.flux[3][1] > s.flux[0][1]);
}

/// 燃料・ギャップ・減速材のピンセルで、ギャップを `gap` にした解
fn with_gap(gap: &Material) -> Solution {
    let cell = PinCell::new(1.26, vec![0.3, 0.4]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 16, 0.02).unwrap();
    let (f, m) = (fuel(), moderator());
    solve(&tracking, &[&f, gap, &m], &Options::default()).unwrap()
}

#[test]
fn void_gap_matches_a_nearly_transparent_gap() {
    let void = with_gap(&Material::void(2));
    // 吸収も散乱もほとんどない物質（Σt だけがごく小さい）
    let thin = with_gap(&two_group(
        [1.0e-9, 1.0e-9],
        [[0.0, 0.0], [0.0, 0.0]],
        [0.0, 0.0],
    ));
    assert!((void.k - thin.k).abs() < 1e-7, "{} vs {}", void.k, thin.k);
    for g in 0..2 {
        let (a, b) = (void.flux[1][g], thin.flux[1][g]);
        assert!((a - b).abs() < 1e-6 * b, "group {g}: {a} vs {b}");
    }
}

#[test]
fn rejects_partly_void_materials() {
    let cell = PinCell::new(1.26, vec![0.3, 0.4]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.05).unwrap();
    let (f, m) = (fuel(), moderator());
    let partly = two_group([1.0, 0.0], [[0.5, 0.0], [0.0, 0.0]], [0.0, 0.0]);
    let result = solve(&tracking, &[&f, &partly, &m], &Options::default());
    assert!(
        matches!(result, Err(MocError::NonPositiveTotal)),
        "{result:?}"
    );
}

#[test]
fn rejects_cells_that_are_all_void() {
    // ピンと、全体が真空のセル
    let pin = PinCell::new(1.26, vec![0.4]).unwrap();
    let empty = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::new(vec![pin, empty], vec![vec![0, 1], vec![1, 0]]).unwrap();
    let tracking = Tracking::new(&lattice, 8, 0.05).unwrap();
    let (f, m, v) = (fuel(), moderator(), Material::void(2));
    let result = solve(&tracking, &[&f, &m, &v, &v, &f, &m], &Options::default());
    assert!(matches!(result, Err(MocError::VoidCell(1))), "{result:?}");
}

#[test]
fn scattering_without_total_is_not_void() {
    let cell = PinCell::new(1.26, vec![0.3, 0.4]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.05).unwrap();
    let (f, m) = (fuel(), moderator());
    let broken = two_group([0.0, 0.0], [[0.1, 0.0], [0.0, 0.0]], [0.0, 0.0]);
    let result = solve(&tracking, &[&f, &broken, &m], &Options::default());
    assert!(
        matches!(result, Err(MocError::NonPositiveTotal)),
        "{result:?}"
    );
}

#[test]
fn reflective_cell_has_no_leakage() {
    assert_eq!(heterogeneous(8, 0.1).leakage, vec![0.0, 0.0]);
}

#[test]
fn vacuum_lattice_balances_production_with_absorption_and_leakage() {
    // 2×2 の格子（燃料のピン2本と水のセル2つ）の外周を真空にすると、漏れの分だけ k が下がる
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap();
    let (f, m) = (fuel(), moderator());
    let materials = [&f, &m, &m, &m, &f, &m];
    let solve_with = |boundaries| {
        let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, boundaries).unwrap();
        solve(&tracking, &materials, &Options::default()).unwrap()
    };
    let all_vacuum = Boundaries {
        left: Boundary::Vacuum,
        right: Boundary::Vacuum,
        top: Boundary::Vacuum,
        bottom: Boundary::Vacuum,
    };
    let s = solve_with(all_vacuum);
    let (mut absorption, mut production) = (0.0, 0.0);
    for (i, m) in materials.iter().enumerate() {
        for g in 0..2 {
            absorption += s.volumes[i] * m.absorption()[g] * s.flux[i][g];
            production += s.volumes[i] * m.nu_fission[g] * s.flux[i][g];
        }
    }
    let leakage: f64 = s.leakage.iter().sum();
    assert!(leakage > 0.0);
    // 反復は相対変化で打ち切るので、収支の誤差は全体の量（生成/k）に比例する
    assert!(
        (absorption + leakage - production / s.k).abs() < 1e-6 * production / s.k,
        "{absorption} + {leakage} vs {}",
        production / s.k
    );
    assert!(s.k < solve_with(Boundaries::REFLECTIVE).k);
}

#[test]
fn albedo_sides_lie_between_vacuum_and_reflection() {
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap();
    let (f, m) = (fuel(), moderator());
    let materials = [&f, &m, &m, &m, &f, &m];
    let solve_with = |right: Boundary| {
        let boundaries = Boundaries {
            right,
            ..Boundaries::REFLECTIVE
        };
        let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, boundaries).unwrap();
        solve(&tracking, &materials, &Options::default()).unwrap()
    };
    let reflective = solve_with(Boundary::Reflective);
    let vacuum = solve_with(Boundary::Vacuum);
    // α = 1 は鏡面反射、α = 0 は真空と同じ
    assert!((solve_with(Boundary::Albedo(1.0)).k - reflective.k).abs() < 1e-9);
    assert!((solve_with(Boundary::Albedo(0.0)).k - vacuum.k).abs() < 1e-9);
    // 間のアルベドでは k も間で、返らなかった分は漏れとして収支が合う
    let half = solve_with(Boundary::Albedo(0.5));
    assert!(vacuum.k < half.k && half.k < reflective.k, "{}", half.k);
    let (mut absorption, mut production) = (0.0, 0.0);
    for (i, m) in materials.iter().enumerate() {
        for g in 0..2 {
            absorption += half.volumes[i] * m.absorption()[g] * half.flux[i][g];
            production += half.volumes[i] * m.nu_fission[g] * half.flux[i][g];
        }
    }
    let leakage: f64 = half.leakage.iter().sum();
    assert!(
        (absorption + leakage - production / half.k).abs() < 1e-6 * production / half.k,
        "{absorption} + {leakage} vs {}",
        production / half.k
    );
}

#[test]
fn result_does_not_depend_on_the_thread_count() {
    // スイープの和を足す順がスレッドの分け方で変わると、結果が最下位のビットで揺れる
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap();
    let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
    let (f, m) = (fuel(), moderator());
    let materials = [&f, &m, &m, &m, &f, &m];
    let solve_on = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| solve(&tracking, &materials, &Options::default()).unwrap())
    };
    let serial = solve_on(1);
    for threads in [2, 5, 8] {
        assert_eq!(solve_on(threads), serial, "{threads} threads");
    }
}

/// P1〜P(L) の散乱行列を、P0 の行列に `factors` の係数を掛けて付けた物質
fn with_moments(m: &Material, factors: &[f64]) -> Material {
    let mut m = m.clone();
    m.anisotropic = factors
        .iter()
        .map(|f| {
            m.scatter
                .iter()
                .map(|row| row.iter().map(|s| f * s).collect())
                .collect()
        })
        .collect();
    m
}

fn order(scattering_order: usize) -> Options {
    Options {
        scattering_order,
        ..Options::default()
    }
}

/// 燃料のピンと水のセルを市松に並べた n×n の格子（外周は真空）
fn vacuum_checkerboard(n: usize) -> Tracking {
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let layout = (0..n)
        .map(|i| (0..n).map(|j| (i + j) % 2).collect())
        .collect();
    let lattice = Lattice::new(vec![pin, water], layout).unwrap();
    let vacuum = Boundaries {
        left: Boundary::Vacuum,
        right: Boundary::Vacuum,
        top: Boundary::Vacuum,
        bottom: Boundary::Vacuum,
    };
    Tracking::with_boundaries(&lattice, 16, 0.05, vacuum).unwrap()
}

/// 市松の格子の領域ごとの物質（ピンは燃料と水、水のセルは水）
fn checkerboard_materials<'a>(n: usize, f: &'a Material, m: &'a Material) -> Vec<&'a Material> {
    (0..n * n)
        .flat_map(|c| {
            if (c / n + c % n).is_multiple_of(2) {
                vec![f, m]
            } else {
                vec![m]
            }
        })
        .collect()
}

#[test]
fn isotropic_moments_reproduce_p0() {
    // 高次の散乱行列が 0 なら、高次のモーメントは線源に入らず、P0 と同じ解（丸め誤差まで）
    let tracking = vacuum_checkerboard(2);
    let (f, m) = (
        with_moments(&fuel(), &[0.0, 0.0]),
        with_moments(&moderator(), &[0.0, 0.0]),
    );
    let materials = checkerboard_materials(2, &f, &m);
    let p0 = solve(&tracking, &materials, &order(0)).unwrap();
    let p2 = solve(&tracking, &materials, &order(2)).unwrap();
    assert!((p2.k - p0.k).abs() < 1e-12 * p0.k, "{} vs {}", p2.k, p0.k);
    for (a, b) in p2.flux.iter().flatten().zip(p0.flux.iter().flatten()) {
        assert!((a - b).abs() < 1e-10 * b, "{a} vs {b}");
    }
}

#[test]
fn anisotropy_does_not_change_an_infinite_medium() {
    // 角度中性子束が等方なので、高次の散乱は効かない
    let cell = PinCell::new(1.26, vec![0.3, 0.45]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.1).unwrap();
    let m = with_moments(&fuel(), &[0.4, 0.1, 0.02]);
    let solution = solve(&tracking, &[&m, &m, &m], &order(3)).unwrap();
    let expected = infinite_medium_k(&m);
    assert!(
        (solution.k - expected).abs() < 1e-6,
        "{} vs {expected}",
        solution.k
    );
}

#[test]
fn forward_scattering_acts_like_the_transport_correction() {
    // 前方に偏った散乱では、外向きに流れる中性子が散乱しても外向きに進み続けるので、等方散乱より
    // 漏れが増えて k が下がる。その効果を拡散の近似で表すのが輸送補正（Σt と自群の散乱から
    // Σ_g' σ1(g→g') を引いた P0。拡散係数が大きくなる）なので、P1 の k は P0 より輸送補正の k に近い
    let n = 6;
    let tracking = vacuum_checkerboard(n);
    let (f, m) = (
        with_moments(&fuel(), &[0.05]),
        with_moments(&moderator(), &[0.5]),
    );
    let corrected = |m: &Material| {
        let mut c = m.clone();
        for g in 0..2 {
            let p1: f64 = m.anisotropic[0][g].iter().sum();
            c.total[g] -= p1;
            c.scatter[g][g] -= p1;
        }
        c.anisotropic.clear();
        c
    };
    let (fc, mc) = (corrected(&f), corrected(&m));
    let k = |f, m, options| {
        solve(&tracking, &checkerboard_materials(n, f, m), &options)
            .unwrap()
            .k
    };
    let (p0, p1, tc) = (
        k(&f, &m, order(0)),
        k(&f, &m, order(1)),
        k(&fc, &mc, order(0)),
    );
    assert!(p1 < p0, "P1 {p1} vs P0 {p0}");
    assert!(
        (p1 - tc).abs() < (p0 - tc).abs(),
        "P0 {p0}, P1 {p1}, transport-corrected {tc}"
    );
}

#[test]
fn anisotropic_result_does_not_depend_on_the_thread_count() {
    let tracking = vacuum_checkerboard(2);
    let (f, m) = (
        with_moments(&fuel(), &[0.05, 0.01]),
        with_moments(&moderator(), &[0.5, 0.2]),
    );
    let materials = checkerboard_materials(2, &f, &m);
    let solve_on = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| solve(&tracking, &materials, &order(2)).unwrap())
    };
    let serial = solve_on(1);
    for threads in [2, 5] {
        assert_eq!(solve_on(threads), serial, "{threads} threads");
    }
}

#[test]
fn available_order_is_the_lowest_among_materials() {
    let (f, m, void) = (
        with_moments(&fuel(), &[0.05, 0.01]),
        with_moments(&moderator(), &[0.5]),
        Material::void(2),
    );
    assert_eq!(available_order(&[&f, &void]), 2);
    assert_eq!(available_order(&[&f, &m, &void]), 1);
    assert_eq!(available_order(&[&f, &fuel()]), 0);
    assert_eq!(available_order(&[&void]), 0);
}

#[test]
fn scattering_order_needs_the_moments() {
    let tracking = vacuum_checkerboard(2);
    let (f, m) = (with_moments(&fuel(), &[0.05]), moderator());
    let materials = checkerboard_materials(2, &f, &m);
    assert!(matches!(
        solve(&tracking, &materials, &order(1)),
        Err(MocError::ScatteringOrder {
            region: 1,
            available: 0,
            requested: 1
        })
    ));
    // 足りなければ使わない次数は無視する
    let m = with_moments(&moderator(), &[0.5]);
    let materials = checkerboard_materials(2, &f, &m);
    assert!(solve(&tracking, &materials, &order(1)).is_ok());
}

#[test]
fn sweep_currents_balance_each_coarse_cell() {
    // 1回のスイープで、セルの中の領域の減衰の和 Σ Δ（= 入る - 出る）と、セルの面から外向きの
    // 正味の中性子流の和は打ち消し合う。線源と入ってくる角度中性子束は任意の値でよい
    let pin = PinCell::new(1.26, vec![0.3, 0.4096]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap();
    let boundaries = Boundaries {
        right: Boundary::Vacuum,
        bottom: Boundary::Vacuum,
        ..Boundaries::REFLECTIVE
    };
    let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, boundaries).unwrap();
    let mesh = CoarseMesh::new(&tracking);
    // 2×2 の格子の中心を通るトラックがあり、セルの角を斜めに抜ける場合も確かめる
    let corner = |track: &crate::tracking::Track| {
        track.segments.windows(2).any(|w| {
            let (a, b) = (mesh.cell(w[0].region), mesh.cell(w[1].region));
            a / 2 != b / 2 && a % 2 != b % 2
        })
    };
    assert!(tracking.tracks.iter().any(corner));
    let (f, m) = (fuel(), moderator());
    let materials = [&f, &f, &m, &m, &m, &f, &f, &m];
    let ng = 2;
    let sweeper = Sweeper::new(
        &tracking,
        &materials,
        ng,
        &mesh,
        &[(0, 0)],
        Representatives::all(&tracking, &mesh, 1),
    );
    let equilibrium: Vec<Vec<f64>> = (0..materials.len())
        .map(|i| {
            (0..ng)
                .map(|g| 0.1 + ((7 * i + 3 * g) % 11) as f64 * 0.01)
                .collect()
        })
        .collect();
    let incoming: Vec<f64> = (0..sweeper.outgoing_len())
        .map(|j| 0.05 + ((37 * j) % 101) as f64 * 1.0e-3)
        .collect();
    let mut outgoing = vec![0.0; sweeper.outgoing_len()];
    let mut tallies = sweeper.tallies(materials.len());
    sweeper
        .sweep(&equilibrium, &incoming, &mut outgoing, &mut tallies)
        .unwrap();
    let sums = sum_tallies(&sweeper, &tallies);
    for cell in 0..mesh.cells() {
        for g in 0..ng {
            let inside: f64 = lattice
                .cell_regions(cell)
                .map(|i| sums.delta[i * ng + g])
                .sum();
            let outward: f64 = mesh
                .cell_faces(cell)
                .iter()
                .map(|&(face, sign)| sign * sums.current[face * ng + g])
                .sum();
            assert!(
                (inside + outward).abs() < 1e-12 * inside.abs().max(outward.abs()),
                "cell {cell} group {g}: {inside} + {outward}"
            );
        }
    }
}

/// CMFD あり・なしで解き、同じ解に、より少ない反復で収束することを確かめる。
fn assert_cmfd_accelerates(tracking: &Tracking, materials: &[&Material]) {
    let plain = Options {
        acceleration: Acceleration::None,
        ..Options::default()
    };
    let slow = solve(tracking, materials, &plain).unwrap();
    let fast = solve(tracking, materials, &Options::default()).unwrap();
    assert!((fast.k - slow.k).abs() < 1e-6, "{} vs {}", fast.k, slow.k);
    for (a, b) in fast.flux.iter().flatten().zip(slow.flux.iter().flatten()) {
        assert!(relative_change(*a, *b) < 1e-5, "{a} vs {b}");
    }
    assert!(
        5 * fast.iterations <= slow.iterations,
        "{} vs {} iterations",
        fast.iterations,
        slow.iterations
    );
}

#[test]
fn cmfd_accelerates_a_pin_cell() {
    let cell = PinCell::new(1.26, vec![0.25, 0.4096, 0.55]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 16, 0.05).unwrap();
    let (f, m) = (fuel(), moderator());
    assert_cmfd_accelerates(&tracking, &[&f, &f, &m, &m]);
}

#[test]
fn cmfd_accelerates_a_lattice_with_vacuum_sides() {
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let layout = vec![vec![0, 1, 1], vec![1, 0, 1], vec![1, 1, 1]];
    let lattice = Lattice::new(vec![pin, water], layout).unwrap();
    let boundaries = Boundaries {
        right: Boundary::Vacuum,
        bottom: Boundary::Vacuum,
        ..Boundaries::REFLECTIVE
    };
    let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, boundaries).unwrap();
    let (f, m) = (fuel(), moderator());
    // セルの順: ピン, 水, 水 / 水, ピン, 水 / 水, 水, 水（ピンは燃料と水の2領域）
    let materials = [&f, &m, &m, &m, &m, &f, &m, &m, &m, &m, &m];
    assert_cmfd_accelerates(&tracking, &materials);
}

#[test]
fn cmfd_accelerates_a_one_row_slab_with_vacuum_ends() {
    // 軸方向の平板（ADR-0128）と同じく、均質なセルの 1 行の両端が真空で上下が鏡面反射
    let lattice =
        Lattice::rectilinear(vec![1.0; 8], vec![1.0], vec![Cell::BLOCK], vec![vec![0; 8]]).unwrap();
    let boundaries = Boundaries {
        left: Boundary::Vacuum,
        right: Boundary::Vacuum,
        ..Boundaries::REFLECTIVE
    };
    let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, boundaries).unwrap();
    let (f, m) = (fuel(), moderator());
    let materials = [&m, &m, &f, &f, &f, &f, &m, &m];
    assert_cmfd_accelerates(&tracking, &materials);
}

#[test]
fn a_mirrored_slab_matches_the_full_slab() {
    // 左右対称な 1 行の平板を鏡映で解くと、粗メッシュの代表のチェーンの中央で、奇数の列では中央のセルの両隣が同じ代表に、
    // 偶数の列では中央の面の隣が自分の像になる（ADR-0158）。1 列に立てた平板も同じ
    let (f, m) = (fuel(), moderator());
    for (count, row) in [(7, true), (8, true), (7, false), (8, false)] {
        let (widths, heights, layout) = if row {
            (vec![1.0; count], vec![1.0], vec![vec![0; count]])
        } else {
            (vec![1.0], vec![1.0; count], vec![vec![0]; count])
        };
        let lattice = Lattice::rectilinear(widths, heights, vec![Cell::BLOCK], layout).unwrap();
        let boundaries = if row {
            Boundaries {
                left: Boundary::Vacuum,
                right: Boundary::Vacuum,
                ..Boundaries::REFLECTIVE
            }
        } else {
            Boundaries {
                top: Boundary::Vacuum,
                bottom: Boundary::Vacuum,
                ..Boundaries::REFLECTIVE
            }
        };
        let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, boundaries).unwrap();
        let materials: Vec<&Material> = (0..count)
            .map(|c| if c == 0 || c + 1 == count { &m } else { &f })
            .collect();
        let options = Options::default();
        let full = solve_as(&tracking, &materials, &options, Symmetry::NONE).unwrap();
        let mirrored = solve_as(&tracking, &materials, &options, Symmetry::QUARTER).unwrap();
        assert_same_solution(&full, &mirrored, &options);
    }
}

#[test]
fn cmfd_accelerates_a_lattice_with_a_margin() {
    // 集合体の間の水のギャップのように、格子の周りに余白がある。縁のセルの粗メッシュは余白の分だけ広い
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let layout = vec![vec![0, 1, 0], vec![1, 0, 1], vec![0, 1, 0]];
    let lattice = Lattice::with_margin(vec![pin, water], layout, 0.3).unwrap();
    let boundaries = Boundaries {
        right: Boundary::Vacuum,
        ..Boundaries::REFLECTIVE
    };
    let tracking = Tracking::with_boundaries(&lattice, 8, 0.1, boundaries).unwrap();
    let (f, m) = (fuel(), moderator());
    // セルの順にピン（燃料と水）と水、その後に余白の 4 × 3 + 4 領域は水
    let mut materials = Vec::new();
    for row in [
        [true, false, true],
        [false, true, false],
        [true, false, true],
    ] {
        for pin in row {
            if pin {
                materials.extend([&f, &m]);
            } else {
                materials.push(&m);
            }
        }
    }
    materials.extend(std::iter::repeat_n(&m, 16));
    assert_eq!(materials.len(), lattice.regions());
    assert_cmfd_accelerates(&tracking, &materials);
}

#[test]
fn cmfd_converges_with_an_optically_thick_reflector() {
    // C5G7 の反射体と同じく、熱群で光学的に厚い（Σt h ≈ 4）水のセルが真空境界まで続く。
    // 従来の CMFD はこのような粗メッシュで振動して収束しない
    let pin = PinCell::new(1.26, vec![0.54]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let mut layout = vec![vec![1; 6]; 6];
    for row in layout.iter_mut().take(3) {
        row[..3].fill(0);
    }
    let lattice = Lattice::new(vec![pin, water], layout).unwrap();
    let boundaries = Boundaries {
        right: Boundary::Vacuum,
        bottom: Boundary::Vacuum,
        ..Boundaries::REFLECTIVE
    };
    let tracking = Tracking::with_boundaries(&lattice, 8, 0.1, boundaries).unwrap();
    let f = fuel();
    let reflector = two_group([0.6, 3.2], [[0.5, 0.098], [0.002, 3.18]], [0.0, 0.0]);
    let mut materials = Vec::new();
    for cell in 0..36 {
        let (row, column) = (cell / 6, cell % 6);
        if row < 3 && column < 3 {
            materials.extend([&f, &reflector]);
        } else {
            materials.push(&reflector);
        }
    }
    assert_cmfd_accelerates(&tracking, &materials);
}

#[test]
fn refinement_converges() {
    let coarse = heterogeneous(16, 0.05);
    let fine = heterogeneous(32, 0.02);
    assert!(
        (coarse.k - fine.k).abs() < 1e-3,
        "{} vs {}",
        coarse.k,
        fine.k
    );
}

#[test]
fn no_fission_is_error() {
    let cell = PinCell::new(1.26, vec![0.4]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 8, 0.1).unwrap();
    let m = moderator();
    assert!(solve(&tracking, &[&m, &m], &Options::default()).is_err());
}

#[test]
fn group_without_neutrons_is_error() {
    // 3群目は χ = 0 で、どの群からも散乱で入らないので、中性子束が常に 0 になる
    let three_group = |nu_fission: [f64; 3]| {
        let chi = vec![1.0, 0.0, 0.0];
        Material {
            total: vec![0.5, 1.2, 0.8],
            scatter: vec![
                vec![0.45, 0.02, 0.0],
                vec![0.0, 0.9, 0.0],
                vec![0.0, 0.0, 0.7],
            ],
            anisotropic: Vec::new(),
            fission: nu_fission.iter().map(|v| v / 2.5).collect(),
            nu_fission: nu_fission.to_vec(),
            fission_terms: if nu_fission.iter().any(|&v| v > 0.0) {
                vec![FissionTerm {
                    chi,
                    nu_fission: nu_fission.to_vec(),
                }]
            } else {
                vec![]
            },
        }
    };
    let (f, m) = (three_group([0.008, 0.25, 0.0]), three_group([0.0; 3]));
    let cell = PinCell::new(1.26, vec![0.4]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.1).unwrap();
    for acceleration in [Acceleration::None, Options::default().acceleration] {
        let options = Options {
            acceleration,
            ..Options::default()
        };
        let err = solve(&tracking, &[&f, &m], &options).unwrap_err();
        assert!(matches!(err, MocError::GroupWithoutSource(2)), "{err}");
    }
}

#[test]
fn region_without_tracks_is_error() {
    // 中心の極小の円には、中心からずれたトラックが1本も通らない
    let cell = PinCell::new(1.26, vec![1.0e-6, 0.4]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 4, 0.3).unwrap();
    let m = fuel();
    assert!(tracking.volumes(3)[0] == 0.0);
    let err = solve(&tracking, &[&m, &m, &m], &Options::default())
        .unwrap_err()
        .to_string();
    assert!(err.contains("region 0"), "{err}");
}

#[test]
fn worst_volume_error_finds_the_largest_relative_mismatch() {
    let exact = [1.0, 0.1, 2.0];
    let tracked = [1.001, 0.098, 2.002];
    let (region, error) = worst_volume_error(&exact, &tracked);
    assert_eq!(region, 1);
    assert!((error + 0.02).abs() < 1e-12, "{error}");
}

#[test]
fn cmfd_tolerance_follows_the_moc_change() {
    let cmfd = CmfdOptions::default();
    let moc = Options::default();
    // MOC の許容誤差以下の変化では設定のまま
    assert_eq!(
        inexact_cmfd(&cmfd, &moc, 0.5 * moc.k_tolerance, moc.flux_tolerance),
        cmfd
    );
    // 許容誤差の 100 倍と 30 倍の変化では、k と中性子束の精度をそれぞれ同じ倍率で緩める
    let tuned = inexact_cmfd(
        &cmfd,
        &moc,
        100.0 * moc.k_tolerance,
        30.0 * moc.flux_tolerance,
    );
    assert!((tuned.k_tolerance / cmfd.k_tolerance - 100.0).abs() < 1e-12);
    assert!((tuned.flux_tolerance / cmfd.flux_tolerance - 30.0).abs() < 1e-12);
    assert_eq!(tuned.max_iterations, cmfd.max_iterations);
}

#[test]
fn rejects_tolerances_that_do_not_bound_cmfd() {
    let cell = PinCell::new(1.26, vec![0.4]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.1).unwrap();
    let (f, m) = (fuel(), moderator());
    let cmfd = CmfdOptions::default();
    let with = |k: f64, flux: f64, cmfd: CmfdOptions| Options {
        k_tolerance: k,
        flux_tolerance: flux,
        acceleration: Acceleration::Cmfd(cmfd),
        ..Options::default()
    };
    // MOC の許容誤差が 0、負、NaN、無限大
    let bad = [
        with(0.0, 1e-6, cmfd),
        with(1e-7, -1e-6, cmfd),
        with(f64::NAN, 1e-6, cmfd),
        with(1e-7, f64::INFINITY, cmfd),
        // CMFD の許容誤差が 0 か、MOC より大きい（k と中性子束で別に）
        with(
            1e-7,
            1e-6,
            CmfdOptions {
                k_tolerance: 0.0,
                ..cmfd
            },
        ),
        with(
            1e-7,
            1e-6,
            CmfdOptions {
                k_tolerance: 1e-6,
                ..cmfd
            },
        ),
        with(
            1e-7,
            1e-6,
            CmfdOptions {
                flux_tolerance: 1e-5,
                ..cmfd
            },
        ),
    ];
    for options in bad {
        let result = solve(&tracking, &[&f, &m], &options);
        assert!(
            matches!(result, Err(MocError::Tolerances { .. })),
            "{options:?}: {result:?}"
        );
    }
}

#[test]
fn relative_change_handles_zero_flux() {
    assert_eq!(relative_change(0.0, 0.0), 0.0);
    assert_eq!(relative_change(1.0, 0.0), 1.0);
    assert_eq!(relative_change(0.0, 2.0), 1.0);
    assert!((relative_change(1.1, 1.0) - 0.1 / 1.1).abs() < 1e-15);
}

#[test]
fn remaining_error_estimate() {
    let mut r = Remaining::default();
    assert_eq!(r.update(1.0), f64::INFINITY);
    // 縮小率 0.5 なら、残りは最新の変化と同じだけ
    assert!((r.update(0.5) - 1.0).abs() < 1e-15);
    // 変化が 0 になれば収束
    let mut r = Remaining::default();
    r.update(1.0e-16);
    assert_eq!(r.update(0.0), 0.0);
    let mut r = Remaining::default();
    r.update(0.0);
    assert_eq!(r.update(0.0), 0.0);
}

#[test]
fn material_errors_keep_the_cause_in_the_chain_only() {
    let cell = PinCell::new(1.26, vec![0.4]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 8, 0.1).unwrap();
    let mut m = fuel();
    m.nu_fission[1] *= 2.0;
    let err = solve(&tracking, &[&m, &m], &Options::default()).unwrap_err();
    let cause = std::error::Error::source(&err).expect("the material error is the cause");
    assert!(!err.to_string().contains(&cause.to_string()), "{err}");
}

#[test]
fn inconsistent_material_is_error() {
    let cell = PinCell::new(1.26, vec![0.4]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 8, 0.1).unwrap();
    let mut m = fuel();
    m.fission_terms.clear();
    let err = solve(&tracking, &[&m, &m], &Options::default()).unwrap_err();
    assert!(
        matches!(
            err,
            MocError::Material {
                region: 0,
                source: MaterialError::NuFission
            }
        ),
        "{err}"
    );
}

#[test]
fn region_count_must_match() {
    let cell = PinCell::new(1.26, vec![0.4]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell.clone()), 8, 0.1).unwrap();
    let m = fuel();
    assert!(solve(&tracking, &[&m], &Options::default()).is_err());
}

/// 燃料の均質なブロック（1 cm 角、4 × 4 の領域）を x に `columns` 個並べた帯の、[領域] の x 座標（中心が原点）
fn strip(columns: usize, n: usize) -> (Lattice, Vec<f64>) {
    let lattice = Lattice::rectilinear(
        vec![1.0; columns],
        vec![1.0],
        vec![Cell::Block { nx: n, ny: n }],
        vec![vec![0; columns]],
    )
    .unwrap();
    let x = (0..columns)
        .flat_map(|c| {
            (0..n * n)
                .map(move |r| -0.5 * columns as f64 + c as f64 + ((r % n) as f64 + 0.5) / n as f64)
        })
        .collect();
    (lattice, x)
}

/// P1 の散乱を P0 の 0.3 倍で持つ物質
fn with_p1(mut m: Material) -> Material {
    m.anisotropic = vec![
        m.scatter
            .iter()
            .map(|row| row.iter().map(|s| 0.3 * s).collect())
            .collect(),
    ];
    m
}

#[test]
fn region_currents_vanish_in_an_infinite_medium() {
    let (lattice, _) = strip(2, 4);
    let tracking = Tracking::new(&lattice, 16, 0.05).unwrap();
    let m = fuel();
    let s = solve(&tracking, &[&m; 32], &Options::default()).unwrap();
    assert_eq!(s.region_currents.len(), 32);
    for (r, current) in s.region_currents.iter().enumerate() {
        for (c, component) in current.iter().enumerate() {
            for (g, j) in component.iter().enumerate() {
                assert!(
                    j.abs() < 1e-8 * s.flux[r][g],
                    "region {r} {c} group {g}: {j}"
                );
            }
        }
    }
}

/// 燃料の帯（x ∈ [−2, 2]、右が真空、ほかは鏡面反射。1 cm のブロックを n × n の領域に分け、1つを真空にする）の、
/// x の1次のモーメントの式
///   Σ_r V_r J_x,r = 2 × 1 × J_右 − Σ_r x_r V_r (q_r − Σt,r φ_r)
/// （部分積分と領域ごとのバランス ∫∇·J = V (q − Σt φ)）の [群] の相対差。x_r は領域の中心なので、領域の中での
/// 変化による誤差（領域の幅の2乗に比例）が残る
fn first_moment_mismatch(n: usize, order: usize) -> Vec<f64> {
    let (lattice, x) = strip(4, n);
    let boundaries = Boundaries {
        right: Boundary::Vacuum,
        ..Boundaries::REFLECTIVE
    };
    let tracking = Tracking::with_boundaries(&lattice, 32, 0.02, boundaries).unwrap();
    let (f, void) = (with_p1(fuel()), Material::void(2));
    let materials: Vec<&Material> = (0..4 * n * n)
        .map(|r| if r == n + 1 { &void } else { &f })
        .collect();
    let options = Options {
        scattering_order: order,
        ..Options::default()
    };
    let s = solve(&tracking, &materials, &options).unwrap();
    (0..2)
        .map(|g| {
            let mut moment = 0.0;
            let mut divergence = 0.0;
            for (r, m) in materials.iter().enumerate() {
                let v = s.volumes[r];
                moment += v * s.region_currents[r][0][g];
                let scattered: f64 = (0..2)
                    .map(|from| m.scatter[from][g] * s.flux[r][from])
                    .sum();
                let fission: f64 = m
                    .fission_terms
                    .iter()
                    .map(|t| {
                        t.chi[g]
                            * t.nu_fission
                                .iter()
                                .zip(&s.flux[r])
                                .map(|(n, phi)| n * phi)
                                .sum::<f64>()
                    })
                    .sum();
                let q = scattered + fission / s.k;
                divergence += x[r] * v * (q - m.total[g] * s.flux[r][g]);
            }
            // 真空の領域の中性子流も右向き
            assert!(s.region_currents[n + 1][0][g] > 0.0);
            let expected = 2.0 * s.currents[1][g] - divergence;
            (moment - expected).abs() / expected.abs()
        })
        .collect()
}

#[test]
fn region_currents_satisfy_the_first_moment_balance() {
    // 等方散乱（スイープで数えた中性子流）と P1 の散乱（1次のモーメント）の両方で、領域を細かくすると差は
    // 幅の2乗で減る（0.25 cm で熱群 2%、0.125 cm で 0.6% 程度）
    for order in [0, 1] {
        let (coarse, fine) = (
            first_moment_mismatch(4, order),
            first_moment_mismatch(8, order),
        );
        for g in 0..2 {
            assert!(fine[g] < 1e-2, "order {order} group {g}: {fine:?}");
            assert!(
                coarse[g] > 3.0 * fine[g],
                "order {order} group {g}: {coarse:?} {fine:?}"
            );
        }
    }
}

#[test]
fn region_currents_are_antisymmetric_in_a_mirror_lattice() {
    // 両側が真空の対称な帯では、x の中性子流は左右で符号が逆、y の中性子流は 0（領域ごとには推定の誤差がある）
    let (lattice, _) = strip(4, 4);
    let boundaries = Boundaries {
        left: Boundary::Vacuum,
        right: Boundary::Vacuum,
        ..Boundaries::REFLECTIVE
    };
    let tracking = Tracking::with_boundaries(&lattice, 32, 0.02, boundaries).unwrap();
    let m = fuel();
    let s = solve(&tracking, &[&m; 64], &Options::default()).unwrap();
    let mirror = |r: usize| {
        let (cell, inner) = (r / 16, r % 16);
        (3 - cell) * 16 + (inner / 4) * 4 + 3 - inner % 4
    };
    let scale = s
        .region_currents
        .iter()
        .map(|c| c[0][0].abs())
        .fold(0.0, f64::max);
    for r in 0..64 {
        for g in 0..2 {
            let (a, b) = (
                s.region_currents[r][0][g],
                s.region_currents[mirror(r)][0][g],
            );
            assert!(
                (a + b).abs() < 1e-6 * scale,
                "region {r} group {g}: {a} vs {b}"
            );
            assert!(s.region_currents[r][1][g].abs() < 1e-2 * scale);
        }
    }
    // 1D の問題なので、セルの中の y の中性子流は（領域の推定の誤差で ±0.2% 程度あるが）セルで足すと打ち消す
    for cell in 0..4 {
        for g in 0..2 {
            let sum: f64 = (cell * 16..(cell + 1) * 16)
                .map(|r| s.volumes[r] * s.region_currents[r][1][g])
                .sum();
            assert!(sum.abs() < 1e-8 * scale, "cell {cell} group {g}: {sum}");
        }
    }
}

#[test]
fn reflective_faces_carry_no_net_current_where_track_ends_meet_a_cell_line() {
    // 左の列が燃料、右の列が水の 2 × 2 の格子（鏡面反射）。トラックの端が上と下の辺の上で列の境（x = 0）にちょうど
    // 来ると、出ていく流れと、反射してつながるトラックの入る流れを、境の両側の別の面に数えうる
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![0, 1]]).unwrap();
    let (f, w) = (fuel(), moderator());
    let materials: Vec<&Material> = (0..lattice.cells())
        .flat_map(|cell| {
            if cell % 2 == 0 {
                vec![&f, &w]
            } else {
                vec![&w]
            }
        })
        .collect();
    let tracking = Tracking::new(&lattice, 16, 0.05).unwrap();
    // 前提: 上か下の辺の上で、端が x = 0 にあるトラックがある
    let half = 0.5 * lattice.height();
    let on_line = |p: (f64, f64)| (p.1.abs() - half).abs() < 1e-9 && p.0.abs() < 1e-9;
    assert!(
        tracking
            .tracks
            .iter()
            .any(|t| on_line(t.start) || on_line(t.end)),
        "no track ends on the column line"
    );
    let options = Options::default();
    let s = solve(&tracking, &materials, &options).unwrap();
    let (columns, rows) = (lattice.columns(), lattice.rows());
    for line in [0, rows] {
        for column in 0..columns {
            let face = crate::cmfd::horizontal_face(rows, columns, line, column);
            for g in 0..2 {
                let (current, flux) = (s.face_currents[face][g], s.face_fluxes[face][g]);
                assert!(
                    current.abs() < options.flux_tolerance * flux,
                    "face {face}, group {g}: net current {current} with face flux {flux}"
                );
            }
        }
    }
}

#[test]
fn faces_of_a_mirror_symmetric_lattice_are_symmetric() {
    // 格子の中心（4 つのセルの角）や外周のセルの角をちょうど通るトラックがある、上下と左右で対称な格子。角を通る移動が
    // 丸めで片側の隣のセルを経由したことになると、鏡像の面の量が揃わない（#0051）
    let pin = PinCell::with_sectors(1.26, vec![0.35, 0.4096], 8).unwrap();
    let water = PinCell::grid(1.26, 2).unwrap();
    let layout = vec![
        vec![0, 1, 1, 0],
        vec![0, 0, 0, 0],
        vec![0, 0, 0, 0],
        vec![0, 1, 1, 0],
    ];
    let lattice = Lattice::with_margin(vec![pin, water], layout, 0.04).unwrap();
    let (f, w) = (fuel(), moderator());
    let materials: Vec<&Material> = (0..lattice.cells())
        .flat_map(|cell| match lattice.cell_type(cell) {
            Cell::Pin(pin) if !pin.is_grid() => (0..pin.regions())
                .map(|r| if r < 8 { &f } else { &w })
                .collect::<Vec<_>>(),
            cell => vec![&w; cell.regions()],
        })
        .collect();
    let vacuum = Boundaries {
        left: Boundary::Vacuum,
        right: Boundary::Vacuum,
        top: Boundary::Vacuum,
        bottom: Boundary::Vacuum,
    };
    let tracking = Tracking::with_boundaries(&lattice, 8, 0.05, vacuum).unwrap();
    // 前提: 内側のセルの角をちょうど通るトラックがある（角から直線までの距離が丸めの程度）
    let edges = |sizes: &[f64]| -> Vec<f64> {
        let mut at = -0.5 * sizes.iter().sum::<f64>();
        sizes
            .iter()
            .map(|s| {
                at += s;
                at
            })
            .collect()
    };
    // 内側の線だけ（最後の要素は外周の辺）
    let (xs, ys) = (edges(lattice.widths()), edges(lattice.heights()));
    let (xs, ys) = (&xs[..xs.len() - 1], &ys[..ys.len() - 1]);
    let through_corner = tracking.tracks.iter().any(|t| {
        let d = (t.end.0 - t.start.0, t.end.1 - t.start.1);
        let norm = d.0.hypot(d.1);
        xs.iter().any(|&x| {
            ys.iter()
                .any(|&y| ((x - t.start.0) * d.1 - (y - t.start.1) * d.0).abs() / norm < 1e-12)
        })
    });
    assert!(through_corner, "no track passes through a corner");
    let options = Options::default();
    let s = solve(&tracking, &materials, &options).unwrap();
    let (columns, rows) = (lattice.columns(), lattice.rows());
    // 面の量を、同じ面の同じ群の中性子束の大きさで比べる（正味の中性子流は 0 に近いことがある）
    let check = |face: usize, image: usize, sign: f64| {
        for g in 0..2 {
            let scale = s.face_fluxes[face][g];
            for (name, a, b) in [
                ("flux", s.face_fluxes[face][g], s.face_fluxes[image][g]),
                (
                    "current",
                    s.face_currents[face][g],
                    sign * s.face_currents[image][g],
                ),
            ] {
                assert!(
                    (a - b).abs() < options.flux_tolerance * scale,
                    "{name} of face {face} and {image}, group {g}: {a} vs {b}"
                );
            }
        }
    };
    // 上下の鏡映: 縦の面は行が反転して流れ（x 向き）はそのまま、横の面は線が反転して流れ（y 向き）の符号が変わる
    for row in 0..rows {
        for line in 0..=columns {
            let face = crate::cmfd::vertical_face(columns, row, line);
            check(
                face,
                crate::cmfd::vertical_face(columns, rows - 1 - row, line),
                1.0,
            );
        }
    }
    for line in 0..=rows {
        for column in 0..columns {
            let face = crate::cmfd::horizontal_face(rows, columns, line, column);
            let image = crate::cmfd::horizontal_face(rows, columns, rows - line, column);
            check(face, image, -1.0);
        }
    }
    // 左右の鏡映: 縦の面は線が反転して流れ（x 向き）の符号が変わり、横の面は列が反転して流れはそのまま
    for row in 0..rows {
        for line in 0..=columns {
            let face = crate::cmfd::vertical_face(columns, row, line);
            check(
                face,
                crate::cmfd::vertical_face(columns, row, columns - line),
                -1.0,
            );
        }
    }
    for line in 0..=rows {
        for column in 0..columns {
            let face = crate::cmfd::horizontal_face(rows, columns, line, column);
            let image = crate::cmfd::horizontal_face(rows, columns, line, columns - 1 - column);
            check(face, image, 1.0);
        }
    }
}

/// 1/4 対称だが 1/8 対称でない格子と、その領域の物質（燃料 0、水 1、別の燃料 2、真空 3 の番号）。ピンは扇形 8
/// （燃料、真空のギャップ、水）。`odd` なら 3 × 3 で、真ん中の列は幅 1.0 の 2 × 2 の水のブロック（左右の鏡映の軸が
/// セルの中を通る長方形）、真ん中の行のピンは別の燃料。そうでなければ 4 × 4 で、ピンと 2 × 2 の水を左右と上下で
/// 対称に（対角線では非対称に）置き、周りに水の余白を付けたもの（軸がセルの境を通る）
fn quarter_lattice(odd: bool) -> (Lattice, Vec<usize>) {
    let pin = PinCell::with_sectors(1.26, vec![0.35, 0.4096], 8).unwrap();
    let lattice = if odd {
        Lattice::rectilinear(
            vec![1.26, 1.0, 1.26],
            vec![1.26; 3],
            vec![Cell::Pin(pin), Cell::Block { nx: 2, ny: 2 }],
            vec![vec![0, 1, 0]; 3],
        )
        .unwrap()
    } else {
        let water = PinCell::grid(1.26, 2).unwrap();
        let layout = vec![
            vec![0, 1, 1, 0],
            vec![0, 0, 0, 0],
            vec![0, 0, 0, 0],
            vec![0, 1, 1, 0],
        ];
        Lattice::with_margin(vec![pin, water], layout, 0.04).unwrap()
    };
    let columns = lattice.columns();
    let materials = (0..lattice.cells())
        .flat_map(|cell| match lattice.cell_type(cell) {
            Cell::Pin(pin) if !pin.is_grid() => {
                let fuel = if odd && cell / columns == 1 { 2 } else { 0 };
                (0..pin.regions())
                    .map(|r| [fuel, 3, 1][r / 8])
                    .collect::<Vec<_>>()
            }
            cell => vec![1; cell.regions()],
        })
        .collect();
    (lattice, materials)
}

fn order_options(scattering_order: usize) -> Options {
    Options {
        scattering_order,
        ..Options::default()
    }
}

/// 対称性 `symmetry` を確かめて（形状とトラック）、その対称性で解く。恒等だけなら全体を解く
fn solve_as(
    tracking: &Tracking,
    materials: &[&Material],
    options: &Options,
    symmetry: Symmetry,
) -> Result<Solution, MocError> {
    let verified = require(tracking, symmetry)?;
    solve_symmetric(tracking, materials, options, &verified)
}

/// 全体を解いた解 `full` と、対称性を使って解いた解 `reduced` が、反復の判定の程度で一致すること
fn assert_same_solution(full: &Solution, reduced: &Solution, options: &Options) {
    assert!(
        (full.k - reduced.k).abs() < options.k_tolerance,
        "k {} vs {}",
        full.k,
        reduced.k
    );
    // 中性子流は、同じ場所・同じ群の中性子束の大きさに対して比べる（鏡面反射の辺や対称線の上の中性子流は丸めの
    // 程度の 0 なので、自分の大きさに対しては比べられない）。`scale` は [位置][群] の尺度
    let close = |a: &[Vec<f64>], b: &[Vec<f64>], scale: &[Vec<f64>], what: &str| {
        for (i, ((x, y), s)) in a.iter().zip(b).zip(scale).enumerate() {
            for (g, ((x, y), s)) in x.iter().zip(y).zip(s).enumerate() {
                assert!(
                    (x - y).abs() <= options.flux_tolerance * s.abs(),
                    "{what} {i} group {g}: {x} vs {y} (scale {s})"
                );
            }
        }
    };
    for (r, (a, b)) in full.flux.iter().zip(&reduced.flux).enumerate() {
        for (g, (a, b)) in a.iter().zip(b).enumerate() {
            assert!(
                relative_change(*a, *b) < options.flux_tolerance,
                "region {r} group {g}: {a} vs {b}"
            );
        }
    }
    let faces = &full.face_fluxes;
    close(faces, &reduced.face_fluxes, faces, "face flux");
    close(
        &full.face_currents,
        &reduced.face_currents,
        faces,
        "face current",
    );
    close(&full.sides, &reduced.sides, &full.sides, "side flux");
    close(
        &full.currents,
        &reduced.currents,
        &full.sides,
        "side current",
    );
    // 漏れは外周にわたる和なので、辺の中性子束 × 辺の長さの和に対して（辺は左、右、上、下）
    let [width, height] = full.size;
    let boundary: Vec<f64> = (0..full.leakage.len())
        .map(|g| {
            [height, height, width, width]
                .iter()
                .zip(&full.sides)
                .map(|(l, side)| l * side[g])
                .sum()
        })
        .collect();
    close(
        std::slice::from_ref(&full.leakage),
        std::slice::from_ref(&reduced.leakage),
        std::slice::from_ref(&boundary),
        "leakage",
    );
    // 領域の中性子流は、領域の中性子束に対して（[領域][x, y の順に群]）
    let currents =
        |s: &Solution| -> Vec<Vec<f64>> { s.region_currents.iter().map(|c| c.concat()).collect() };
    let fluxes: Vec<Vec<f64>> = full
        .flux
        .iter()
        .map(|f| [&f[..], &f[..]].concat())
        .collect();
    close(
        &currents(full),
        &currents(reduced),
        &fluxes,
        "region current",
    );
    assert_eq!(full.volumes, reduced.volumes);
}

#[test]
fn quarter_symmetry_reproduces_the_full_solution() {
    let (f, w, g, void) = (
        with_p1(fuel()),
        with_p1(moderator()),
        with_p1(two_group(
            [0.55, 1.4],
            [[0.5, 0.03], [0.0, 1.1]],
            [0.01, 0.3],
        )),
        Material::void(2),
    );
    let library = [&f, &w, &g, &void];
    let vacuum = Boundaries {
        left: Boundary::Vacuum,
        right: Boundary::Vacuum,
        top: Boundary::Vacuum,
        bottom: Boundary::Vacuum,
    };
    // 左右だけ真空（1/8 対称ではない）
    let sides = Boundaries {
        left: Boundary::Vacuum,
        right: Boundary::Vacuum,
        ..Boundaries::REFLECTIVE
    };
    // アルベドの辺（左右と上下で違う値）
    let albedo = Boundaries {
        left: Boundary::Albedo(0.3),
        right: Boundary::Albedo(0.3),
        top: Boundary::Albedo(0.6),
        bottom: Boundary::Albedo(0.6),
    };
    // 左右がアルベドで上下が鏡面反射
    let mixed = Boundaries {
        left: Boundary::Albedo(0.3),
        right: Boundary::Albedo(0.3),
        ..Boundaries::REFLECTIVE
    };
    // 軸がセルの中を通る長方形と、セルの境を通る正方形。鏡面反射、真空、アルベド、等方散乱と P1
    for odd in [true, false] {
        let (lattice, index) = quarter_lattice(odd);
        let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
        for boundaries in [Boundaries::REFLECTIVE, vacuum, sides, albedo, mixed] {
            let tracking = Tracking::with_boundaries(&lattice, 8, 0.07, boundaries).unwrap();
            for order in [0, 1] {
                let full =
                    solve_as(&tracking, &materials, &order_options(order), Symmetry::NONE).unwrap();
                let options = order_options(order);
                let quarter = solve_as(&tracking, &materials, &options, Symmetry::QUARTER).unwrap();
                assert_same_solution(&full, &quarter, &options);
            }
        }
    }
}

/// 正方形の 10 通りの対称性（ADR-0112）
fn all_symmetries() -> Vec<Symmetry> {
    [
        &[][..],
        &[2],
        &[4],
        &[5],
        &[6],
        &[7],
        &[2, 4, 5],
        &[2, 6, 7],
        &[1, 2, 3],
        &[1, 2, 3, 4, 5, 6, 7],
    ]
    .iter()
    .map(|ops| Symmetry::new(ops).unwrap())
    .collect()
}

/// 1/8 対称な格子と、その領域の物質（燃料 0、水 1、別の燃料 2、真空 3 の番号）。ピンは扇形 8（燃料、真空のギャップ、
/// 水）。`odd` なら 3 × 3 で中央のピンは別の燃料（軸と対角線がセルの中を通る）。そうでなければ 4 × 4 で、辺の中央の
/// 8 つが 2 × 2 の水、周りに水の余白（軸がセルの境を通る）
fn octant_lattice(odd: bool) -> (Lattice, Vec<usize>) {
    let pin = PinCell::with_sectors(1.26, vec![0.35, 0.4096], 8).unwrap();
    let lattice = if odd {
        Lattice::new(vec![pin], vec![vec![0; 3]; 3]).unwrap()
    } else {
        let water = PinCell::grid(1.26, 2).unwrap();
        let layout = vec![
            vec![0, 1, 1, 0],
            vec![1, 0, 0, 1],
            vec![1, 0, 0, 1],
            vec![0, 1, 1, 0],
        ];
        Lattice::with_margin(vec![pin, water], layout, 0.04).unwrap()
    };
    let materials = (0..lattice.cells())
        .flat_map(|cell| match lattice.cell_type(cell) {
            Cell::Pin(pin) if !pin.is_grid() => {
                let fuel = if odd && cell == 4 { 2 } else { 0 };
                (0..pin.regions())
                    .map(|r| [fuel, 3, 1][r / 8])
                    .collect::<Vec<_>>()
            }
            cell => vec![1; cell.regions()],
        })
        .collect();
    (lattice, materials)
}

#[test]
fn every_symmetry_reproduces_the_full_solution() {
    let (f, w, g, void) = (
        with_p1(fuel()),
        with_p1(moderator()),
        with_p1(two_group(
            [0.55, 1.4],
            [[0.5, 0.03], [0.0, 1.1]],
            [0.01, 0.3],
        )),
        Material::void(2),
    );
    let library = [&f, &w, &g, &void];
    let vacuum = Boundaries {
        left: Boundary::Vacuum,
        right: Boundary::Vacuum,
        top: Boundary::Vacuum,
        bottom: Boundary::Vacuum,
    };
    let albedo = Boundaries {
        left: Boundary::Albedo(0.3),
        right: Boundary::Albedo(0.3),
        top: Boundary::Albedo(0.3),
        bottom: Boundary::Albedo(0.3),
    };
    // 1/8 対称な格子では 10 通りの対称性がどれも成り立つ。方位角 8（45° なし）は等方散乱で、10（45° あり）は P1 で。間隔
    // 0.05 cm では、方位角 10 の 45° のトラックが余白の角のセル（0.04 cm 角）の面を通らないので 0.04 cm にする（0.04 cm で通ることを確かめた）
    for odd in [true, false] {
        let (lattice, index) = octant_lattice(odd);
        let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
        for boundaries in [Boundaries::REFLECTIVE, vacuum, albedo] {
            for (azimuths, order, spacing) in [(8, 0, 0.05), (10, 1, 0.04)] {
                let tracking =
                    Tracking::with_boundaries(&lattice, azimuths, spacing, boundaries).unwrap();
                let full =
                    solve_as(&tracking, &materials, &order_options(order), Symmetry::NONE).unwrap();
                for symmetry in all_symmetries() {
                    let options = order_options(order);
                    let reduced = solve_as(&tracking, &materials, &options, symmetry)
                        .unwrap_or_else(|e| panic!("{}: {e}", symmetry.name()));
                    assert_same_solution(&full, &reduced, &options);
                }
            }
        }
    }
}

/// 3 × 3 の、扇形 8 のピン（燃料、真空のギャップ、水）の中央に、扇形ごとの物質が `pattern`（`library` の番号）の
/// 扇形 8 のセルを置いた格子の領域の物質（`octant_lattice(true)` の中央を置き換えたもの）
fn patterned<'a>(pattern: [usize; 8], library: &[&'a Material]) -> (Lattice, Vec<&'a Material>) {
    let pin = PinCell::with_sectors(1.26, vec![0.35, 0.4096], 8).unwrap();
    let center = PinCell::with_sectors(1.26, vec![], 8).unwrap();
    let layout = vec![vec![0, 0, 0], vec![0, 1, 0], vec![0, 0, 0]];
    let lattice = Lattice::new(vec![pin, center], layout).unwrap();
    let materials = (0..9)
        .flat_map(|cell| {
            if cell == 4 {
                pattern.iter().map(|&i| library[i]).collect::<Vec<_>>()
            } else {
                (0..24).map(|r| library[[0, 3, 1][r / 8]]).collect()
            }
        })
        .collect();
    (lattice, materials)
}

#[test]
fn each_symmetry_alone_reproduces_the_full_solution() {
    // それぞれの対称性だけが成り立つ格子（中央のセルの扇形の物質の並びと境界条件で、ほかの操作を外す）を、物質の参照で
    // 判定した対称性で解き、全体を解いた解と比べる（P1、方位角 10）
    let (f, w, void) = (with_p1(fuel()), with_p1(moderator()), Material::void(2));
    let variant =
        |total: [f64; 2]| with_p1(two_group(total, [[0.5, 0.098], [0.002, 1.97]], [0.0, 0.0]));
    let (a, b, c) = (
        variant([0.65, 2.1]),
        variant([0.7, 2.2]),
        variant([0.75, 2.3]),
    );
    // 番号: 燃料 0、水 1、水の変種 2、4、5、真空 3（ピンのギャップ）
    let library = [&f, &w, &a, &void, &b, &c];
    let (left, right, bottom) = (
        Boundaries {
            left: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        },
        Boundaries {
            right: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        },
        Boundaries {
            bottom: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        },
    );
    let corner = Boundaries {
        right: Boundary::Vacuum,
        ..bottom
    };
    let opposite = Boundaries {
        left: Boundary::Vacuum,
        ..bottom
    };
    let reflective = Boundaries::REFLECTIVE;
    let symmetry = |ops: &[usize]| Symmetry::new(ops).unwrap();
    let cases = [
        ([1; 8], right, symmetry(&[5])),
        ([1; 8], bottom, symmetry(&[4])),
        ([1; 8], corner, symmetry(&[6])),
        ([1; 8], opposite, symmetry(&[7])),
        ([1, 2, 4, 5, 1, 2, 4, 5], reflective, symmetry(&[2])),
        ([1, 2, 2, 1, 1, 2, 2, 1], reflective, Symmetry::QUARTER),
        ([1, 1, 2, 2, 1, 1, 2, 2], reflective, symmetry(&[2, 6, 7])),
        ([1, 2, 1, 2, 1, 2, 1, 2], reflective, symmetry(&[1, 2, 3])),
        ([1, 2, 1, 2, 1, 2, 1, 2], left, Symmetry::NONE),
    ];
    let options = order_options(1);
    for (pattern, boundaries, expected) in cases {
        let (lattice, materials) = patterned(pattern, &library);
        let tracking = Tracking::with_boundaries(&lattice, 10, 0.05, boundaries).unwrap();
        let identity: Vec<usize> = materials
            .iter()
            .map(|&m| std::ptr::from_ref(m) as usize)
            .collect();
        let detected = detect(&tracking, &identity).unwrap();
        assert_eq!(
            detected.symmetry(),
            expected,
            "{pattern:?}: {} instead of {}",
            detected.symmetry().name(),
            expected.name()
        );
        let full = solve(&tracking, &materials, &options).unwrap();
        let reduced = solve_symmetric(&tracking, &materials, &options, &detected).unwrap();
        assert_same_solution(&full, &reduced, &options);
    }
}

#[test]
fn axis_swapping_symmetries_keep_high_moments() {
    // P2 と P3 では、軸を入れ替える操作でモーメント (l, m) と (l, -m) が入れ替わり、符号が付く。1/8 対称な格子を、軸を
    // 入れ替える操作を含む対称性で解き、全体を解いた解と比べる
    let (f, w, g, void) = (
        with_moments(&fuel(), &[0.3, 0.15, 0.05]),
        with_moments(&moderator(), &[0.4, 0.2, 0.1]),
        with_moments(&fuel(), &[0.2, 0.1, 0.05]),
        Material::void(2),
    );
    let library = [&f, &w, &g, &void];
    let (lattice, index) = octant_lattice(true);
    let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
    let tracking = Tracking::new(&lattice, 10, 0.05).unwrap();
    let swapping = [
        Symmetry::new(&[6]).unwrap(),
        Symmetry::new(&[7]).unwrap(),
        Symmetry::new(&[2, 6, 7]).unwrap(),
        Symmetry::new(&[1, 2, 3]).unwrap(),
        Symmetry::OCTANT,
    ];
    for order in [2, 3] {
        let options = order_options(order);
        let full = solve(&tracking, &materials, &options).unwrap();
        for symmetry in swapping {
            let reduced = solve_as(&tracking, &materials, &options, symmetry)
                .unwrap_or_else(|e| panic!("{}: {e}", symmetry.name()));
            assert_same_solution(&full, &reduced, &options);
        }
    }
}

#[test]
fn axis_swapping_symmetries_need_a_square_lattice() {
    // 1/4 対称な長方形（真ん中の列が幅 1.0）では、軸を入れ替えない対称性は解けて全体の解と一致し、軸を入れ替える
    // 操作を含む対称性はエラー
    let (lattice, index) = quarter_lattice(true);
    let (f, w, void) = (fuel(), moderator(), Material::void(2));
    let library = [&f, &w, &f, &void];
    let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
    let tracking = Tracking::new(&lattice, 8, 0.05).unwrap();
    let full = solve_as(&tracking, &materials, &order_options(0), Symmetry::NONE).unwrap();
    for symmetry in all_symmetries() {
        let options = order_options(0);
        let result = solve_as(&tracking, &materials, &options, symmetry);
        let swaps = symmetry
            .operations()
            .iter()
            .any(|&i| [1, 3, 6, 7].contains(&i));
        if swaps {
            assert!(
                matches!(result, Err(MocError::Symmetry(_))),
                "{}: {result:?}",
                symmetry.name()
            );
        } else {
            assert_same_solution(&full, &result.unwrap(), &options);
        }
    }
}

#[test]
fn octant_result_does_not_depend_on_the_thread_count() {
    let (lattice, index) = octant_lattice(true);
    let (f, w, void) = (with_p1(fuel()), with_p1(moderator()), Material::void(2));
    let library = [&f, &w, &f, &void];
    let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
    let tracking = Tracking::new(&lattice, 10, 0.05).unwrap();
    let options = order_options(1);
    let solve_on = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| solve_as(&tracking, &materials, &options, Symmetry::OCTANT).unwrap())
    };
    let serial = solve_on(1);
    for threads in [2, 5] {
        assert_eq!(solve_on(threads), serial, "{threads} threads");
    }
}

#[test]
fn quarter_symmetry_solves_only_representatives() {
    // 軸がセルの境を通る格子では、どの領域とセルも 4 つずつ移り合う
    let (lattice, index) = quarter_lattice(false);
    let (f, w, void) = (fuel(), moderator(), Material::void(2));
    let library = [&f, &w, &f, &void];
    let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
    let tracking = Tracking::new(&lattice, 8, 0.05).unwrap();
    let mesh = CoarseMesh::new(&tracking);
    let reduced = Representatives::by_images(
        &tracking,
        &materials,
        &mesh,
        &[(0, 0)],
        require(&tracking, Symmetry::QUARTER).unwrap().images(),
    )
    .unwrap();
    assert_eq!(4 * reduced.regions.representatives.len(), lattice.regions());
    assert_eq!(4 * reduced.cells.representatives.len(), lattice.cells());
}

#[test]
fn quarter_result_does_not_depend_on_the_thread_count() {
    // 代表のトラックの組と、代表の領域に集める和の順は、スレッドの分け方によらない
    let (lattice, index) = quarter_lattice(true);
    let (f, w, void) = (with_p1(fuel()), with_p1(moderator()), Material::void(2));
    let library = [&f, &w, &f, &void];
    let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
    let tracking = Tracking::new(&lattice, 8, 0.05).unwrap();
    let options = order_options(1);
    let solve_on = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| solve_as(&tracking, &materials, &options, Symmetry::QUARTER).unwrap())
    };
    let serial = solve_on(1);
    for threads in [2, 5] {
        assert_eq!(solve_on(threads), serial, "{threads} threads");
    }
}

#[test]
fn quarter_symmetry_rejects_asymmetric_materials() {
    // 左上のピンの1つの扇形だけ物質が違う
    let (lattice, mut index) = quarter_lattice(false);
    let sector = index.iter().position(|&i| i == 0).unwrap();
    index[sector] = 2;
    let (f, w, void) = (fuel(), moderator(), Material::void(2));
    let g = Material {
        total: vec![0.55, 1.4],
        ..fuel()
    };
    let library = [&f, &w, &g, &void];
    let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
    let tracking = Tracking::new(&lattice, 8, 0.05).unwrap();
    assert!(solve_as(&tracking, &materials, &order_options(0), Symmetry::NONE).is_ok());
    let result = solve_as(&tracking, &materials, &order_options(0), Symmetry::QUARTER);
    assert!(
        matches!(result, Err(MocError::AsymmetricMaterial { .. })),
        "{result:?}"
    );
    // 操作は番号ではなく名前で示す
    let message = result.unwrap_err().to_string();
    assert!(
        message.contains("under the rotation by 180 degrees")
            || message.contains("under the mirror across"),
        "{message}"
    );
}

#[test]
fn quarter_symmetry_needs_one_material_for_images() {
    // 値が同じでも別の物質（自己遮蔽で別々に解いたピンにあたる）は、移り合う領域に使えない。値の比較で決めると、
    // 別々に解いたピンが丸めの違いだけで対称かどうかが変わる
    let (lattice, index) = quarter_lattice(false);
    let (f, w, void) = (fuel(), moderator(), Material::void(2));
    let copy = fuel();
    let library = [&f, &w, &f, &void];
    let mut materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
    let last = index.iter().rposition(|&i| i == 0).unwrap();
    materials[last] = &copy;
    let tracking = Tracking::new(&lattice, 8, 0.05).unwrap();
    assert!(solve_as(&tracking, &materials, &order_options(0), Symmetry::NONE).is_ok());
    let result = solve_as(&tracking, &materials, &order_options(0), Symmetry::QUARTER);
    assert!(
        matches!(result, Err(MocError::AsymmetricMaterial { .. })),
        "{result:?}"
    );
    // 事実（別の物質を指している）を示し、自己遮蔽で別々に解かれた可能性はヒントにとどめる
    let message = result.unwrap_err().to_string();
    assert!(
        message.contains("refer to different materials") && message.contains("ADR-0112"),
        "{message}"
    );
}

#[test]
fn quarter_symmetry_rejects_an_asymmetric_boundary() {
    // 右だけ真空: 上下の鏡映しか残らない
    let (lattice, index) = quarter_lattice(false);
    let (f, w, void) = (fuel(), moderator(), Material::void(2));
    let library = [&f, &w, &f, &void];
    let materials: Vec<&Material> = index.iter().map(|&i| library[i]).collect();
    let boundaries = Boundaries {
        right: Boundary::Vacuum,
        ..Boundaries::REFLECTIVE
    };
    let tracking = Tracking::with_boundaries(&lattice, 8, 0.05, boundaries).unwrap();
    let result = solve_as(&tracking, &materials, &order_options(0), Symmetry::QUARTER);
    assert!(matches!(result, Err(MocError::Symmetry(_))), "{result:?}");
}
