use super::*;
use std::f64::consts::PI;

/// どの位置にもノードを置く物質の並び
fn placed(materials: impl IntoIterator<Item = usize>) -> Vec<Option<usize>> {
    materials.into_iter().map(Some).collect()
}

/// 2群の均質な物質（χ = (1, 0)、上方散乱なし）
fn material() -> FewGroup {
    let (d1, d2) = (1.4, 0.4);
    let (a1, a2) = (0.010, 0.080);
    let s12 = 0.017;
    let (nf1, nf2) = (0.006, 0.11);
    let total = [1.0 / (3.0 * d1), 1.0 / (3.0 * d2)];
    FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        // NEM は使わないが、少数群定数として正しい値にする
        k_inf: 1.0,
        total: total.to_vec(),
        absorption: vec![a1, a2],
        nu_fission: vec![nf1, nf2],
        fission: vec![nf1 / 2.4, nf2 / 2.4],
        chi: vec![1.0, 0.0],
        diffusion: vec![d1, d2],
        scatter: vec![vec![total[0] - a1 - s12, s12], vec![0.0, total[1] - a2]],
        adf: vec![vec![1.0, 1.0]; 6],
        vacuum_ratio: vec![None; 6],
        normalized_adf: vec![false; 6],
        height: None,
        pins: None,
    }
}

/// 2群・χ = (1, 0)・上方散乱なしの均質な体系の k（バックリング b2）
fn analytic_k(m: &FewGroup, b2: f64) -> f64 {
    let s12 = m.scatter[0][1];
    let r1 = m.absorption[0] + s12 + m.diffusion[0] * b2;
    let r2 = m.absorption[1] + m.diffusion[1] * b2;
    (m.nu_fission[0] + m.nu_fission[1] * s12 / r2) / r1
}

fn cuboid(size: [f64; 3], nodes: [usize; 3], boundaries: [Boundary; 6]) -> Problem {
    let split = |l: f64, n: usize| vec![l / n as f64; n];
    Problem {
        mesh: Mesh {
            dx: split(size[0], nodes[0]),
            dy: split(size[1], nodes[1]),
            dz: split(size[2], nodes[2]),
        },
        materials: vec![material()],
        reflectors: Vec::new(),
        node_material: placed(vec![0; nodes.iter().product()]),
        boundaries,
        outline: None,
    }
}

#[test]
fn reflective_box_gives_k_infinity() {
    let p = cuboid([30.0, 30.0, 40.0], [2, 3, 2], [Boundary::Reflective; 6]);
    let s = solve(&p, &Options::default()).unwrap();
    let expected = analytic_k(&material(), 0.0);
    assert!((s.k - expected).abs() < 1e-8, "{} vs {expected}", s.k);
}

#[test]
fn bare_cuboid_matches_analytic_buckling() {
    // 誤差は主に横方向の漏れの2次近似によるもので、ノードを細かくすると速く減る
    let size = [60.0, 80.0, 100.0];
    let b2: f64 = size.iter().map(|l| (PI / l).powi(2)).sum();
    let expected = analytic_k(&material(), b2);
    let error = |n: usize| {
        let p = cuboid(size, [n; 3], [Boundary::ZeroFlux; 6]);
        (solve(&p, &Options::default()).unwrap().k - expected).abs() / expected
    };
    let (coarse, fine) = (error(4), error(8));
    assert!(fine < 2e-4, "relative error {fine}");
    assert!(fine < coarse / 8.0, "{coarse} -> {fine}");
}

#[test]
fn symmetric_half_core_matches_full_core() {
    let full = cuboid([60.0, 60.0, 60.0], [6, 4, 4], [Boundary::ZeroFlux; 6]);
    let mut half_boundaries = [Boundary::ZeroFlux; 6];
    half_boundaries[0] = Boundary::Reflective;
    let half = cuboid([30.0, 60.0, 60.0], [3, 4, 4], half_boundaries);
    let (f, h) = (
        solve(&full, &Options::default()).unwrap(),
        solve(&half, &Options::default()).unwrap(),
    );
    assert!((f.k - h.k).abs() < 1e-7, "{} vs {}", f.k, h.k);
}

#[test]
fn neutron_balance_holds() {
    let p = cuboid([60.0, 80.0, 100.0], [3, 4, 5], [Boundary::ZeroFlux; 6]);
    let s = solve(&p, &Options::default()).unwrap();
    let m = material();
    let volume = 60.0 * 80.0 * 100.0 / 60.0;
    let (mut absorption, mut production) = (0.0, 0.0);
    for phi in &s.flux {
        for ((f, a), nf) in phi.iter().zip(&m.absorption).zip(&m.nu_fission) {
            absorption += volume * a * f;
            production += volume * nf * f;
        }
    }
    let leakage: f64 = s.leakage.iter().sum();
    assert!(leakage > 0.0);
    assert!(
        (production / s.k - absorption - leakage).abs() < 1e-6 * production,
        "{production} / {} vs {absorption} + {leakage}",
        s.k
    );
}

/// 解の面の正味の中性子流が、各ノードのバランスを（反復の収束の程度で）満たし、隣のノードと面で一致し、
/// 反射の面で 0、外周の面の和が漏れに等しく、Marshak の面で J = f φ/4·2 = f (o + i) であること
fn assert_face_currents(p: &Problem) {
    let (core, s, state) = solve_state(p);
    let ng = core.groups;
    let scale = s
        .currents
        .iter()
        .flatten()
        .flatten()
        .flatten()
        .fold(0.0, |a: f64, j| a.max(j.abs()));
    assert!(scale > 0.0);
    for n in 0..core.nodes() {
        let m = core.material(n);
        let phi = &s.flux[n];
        let fission: f64 = (0..ng).map(|g| m.nu_fission[g] * phi[g]).sum::<f64>() / s.k;
        for g in 0..ng {
            let leakage: f64 = (0..3)
                .map(|d| {
                    let [minus, plus] = s.currents[d][n][g];
                    (plus - minus) / core.width(n, d)
                })
                .sum();
            let scattering_in: f64 = (0..ng)
                .filter(|&h| h != g)
                .map(|h| m.scatter[h][g] * phi[h])
                .sum();
            let removal = (m.total[g] - m.scatter[g][g]) * phi[g];
            let residual = leakage + removal - scattering_in - m.chi[g] * fission;
            assert!(
                residual.abs() < 1e-6 * removal,
                "node {n} group {g}: {residual}"
            );
        }
    }
    let mut outer = vec![0.0; ng];
    for n in 0..core.nodes() {
        for d in 0..3 {
            let area = core.volume(n) / core.width(n, d);
            for (g, out) in outer.iter_mut().enumerate() {
                let [minus, plus] = s.currents[d][n][g];
                for (step, index, current) in [(-1, 0, -minus), (1, 1, plus)] {
                    if let Some(next) = core.neighbor(n, d, step) {
                        let other = s.currents[d][next][g][1 - index];
                        let own = if step > 0 { plus } else { minus };
                        assert!(
                            (own - other).abs() < 1e-12 * scale,
                            "node {n} direction {d}: {own} vs {other}"
                        );
                        continue;
                    }
                    *out += area * current;
                    match p.boundaries[2 * d + index] {
                        Boundary::Reflective => {
                            assert!(current.abs() < 1e-12 * scale, "reflective face of {n}");
                        }
                        Boundary::Marshak => {
                            let (o, i) = if step > 0 {
                                let (_, i) = core.incoming(&state, n, d, g);
                                (state.out_plus[d][n][g], i)
                            } else {
                                let (i, _) = core.incoming(&state, n, d, g);
                                (state.out_minus[d][n][g], i)
                            };
                            let f = core.adf(n, d, step > 0, g);
                            assert!(
                                (current - f * (o + i)).abs() < 1e-12 * scale,
                                "Marshak face of {n}"
                            );
                        }
                        _ => assert!(current > 0.0, "outer face of {n} must leak"),
                    }
                }
            }
        }
    }
    for (g, (out, leak)) in outer.iter().zip(&s.leakage).enumerate() {
        assert!(
            (out - leak).abs() < 1e-12 * leak.abs().max(1e-300),
            "group {g}: {out} vs {leak}"
        );
    }
}

#[test]
fn face_currents_balance_each_node() {
    // ピン出力の再構成（ADR-0096）は面の正味の中性子流を条件にする
    let mut boundaries = [Boundary::ZeroFlux; 6];
    boundaries[0] = Boundary::Reflective;
    boundaries[3] = Boundary::Reflective;
    assert_face_currents(&cuboid([60.0, 80.0, 100.0], [3, 4, 5], boundaries));
}

#[test]
fn face_currents_are_continuous_across_adf_jumps() {
    // 燃料と反射体（核分裂がないので ANM。ADR-0102）の 2D で、ADF が面で違い、+x は Marshak。
    // 正味の中性子流は ADF の跳びによらず面で連続
    let mut fuel = material();
    fuel.adf = vec![
        vec![1.1, 0.9],
        vec![0.95, 1.2],
        vec![1.05, 0.85],
        vec![1.0, 1.1],
        vec![1.0, 1.0],
        vec![1.0, 1.0],
    ];
    let mut water = reflector();
    water.adf = vec![
        vec![0.8, 1.3],
        vec![1.2, 0.7],
        vec![0.9, 1.4],
        vec![1.1, 0.75],
        vec![1.0, 1.0],
        vec![1.0, 1.0],
    ];
    let mut boundaries = [Boundary::Reflective; 6];
    boundaries[1] = Boundary::Marshak;
    boundaries[2] = Boundary::ZeroFlux;
    let p = Problem {
        mesh: Mesh {
            dx: vec![10.0, 10.0, 10.0, 10.0],
            dy: vec![10.0, 10.0],
            dz: vec![20.0],
        },
        materials: vec![fuel, water],
        reflectors: Vec::new(),
        node_material: placed(vec![0, 0, 1, 1, 0, 1, 1, 1]),
        boundaries,
        outline: None,
    };
    assert_face_currents(&p);
}

#[test]
fn assembly_sized_nodes_converge() {
    // 集合体程度のノード（16.7 cm）の裸の立方体。外周の境界条件の遅れで発散していた
    let size = [200.0; 3];
    let b2: f64 = size.iter().map(|l| (PI / l).powi(2)).sum();
    let expected = analytic_k(&material(), b2);
    let p = cuboid(size, [12; 3], [Boundary::ZeroFlux; 6]);
    let s = solve(&p, &Options::default()).unwrap();
    assert!(
        (s.k - expected).abs() < 1e-4 * expected,
        "{} vs {expected}",
        s.k
    );
}

#[test]
fn single_reflective_node_direction_stays_bounded() {
    // z に 300 cm のノード1つで両側が反射。その方向の解は一様で、展開係数は 0 のまま
    // でなければならない。前の反復の a1, a2 から a3, a4 を求めると、丸め誤差が反復ごとに
    // 約 -5 倍と -2.3 倍になり、400 回ほどで NaN になっていた（#0025）
    let mut boundaries = [Boundary::Reflective; 6];
    boundaries[0] = Boundary::ZeroFlux;
    boundaries[1] = Boundary::ZeroFlux;
    let p = cuboid([60.0, 20.0, 300.0], [6, 1, 1], boundaries);
    let tight = Options {
        k_tolerance: 1e-14,
        flux_tolerance: 1e-13,
        max_iterations: 10_000,
    };
    let converged = solve(&p, &tight).unwrap();
    let mut core = Core::new(&p).unwrap();
    // 許容差 0 では収束しないので、上限の回数まで反復を続ける
    let long = Options {
        k_tolerance: 0.0,
        flux_tolerance: 0.0,
        max_iterations: 2000,
    };
    match iterate(&mut core, &long, None) {
        Err(NemError::NotConverged { k, .. }) => assert!(
            (k - converged.k).abs() < 1e-12 * converged.k,
            "{k} vs {}",
            converged.k
        ),
        Err(e) => panic!("expected NotConverged, got {e}"),
        Ok((s, _)) => panic!("converged with zero tolerances (k = {})", s.k),
    }
}

#[test]
fn non_fundamental_mode_is_error() {
    // 断面 10 cm の裸の細い柱を1ノードで表すと、横方向の漏れが除去の約5倍になり、
    // 符号が交互に変わる解に収束する。基本モードではないので誤りとして返す
    let p = cuboid([10.0, 10.0, 100.0], [1, 1, 10], [Boundary::ZeroFlux; 6]);
    let err = solve(&p, &Options::default()).unwrap_err().to_string();
    assert!(err.contains("positive"), "{err}");
}

#[test]
fn changes_at_rounding_level_count_as_converged() {
    let mut r = Remaining::default();
    r.update(1.0e-15);
    assert_eq!(r.update(2.0e-15), 0.0);
    let mut r = Remaining::default();
    assert_eq!(r.update(1.0), f64::INFINITY);
    // 縮小率 0.5 なら、最新の変化を含めて残りは 2 倍
    assert!((r.update(0.5) - 1.0).abs() < 1e-15);
}

#[test]
fn inconsistent_problem_is_error() {
    let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
    p.node_material.pop();
    assert!(solve(&p, &Options::default()).is_err());
    let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
    p.node_material[0] = Some(1);
    assert!(solve(&p, &Options::default()).is_err());
    let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
    p.mesh.dx[0] = 0.0;
    assert!(solve(&p, &Options::default()).is_err());
    // 有限でない断面積は、反復に入る前にエラーにする
    let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
    p.materials[0].nu_fission[1] = f64::NAN;
    let err = solve(&p, &Options::default()).unwrap_err().to_string();
    assert!(err.contains("finite"), "{err}");
    let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
    p.materials[0].diffusion[0] = f64::INFINITY;
    assert!(
        solve(&p, &Options::default())
            .unwrap_err()
            .to_string()
            .contains("finite")
    );
}

/// 1群の均質な物質（拡散係数 1 cm）
fn one_group() -> FewGroup {
    let (d, absorption, nu_fission) = (1.0, 0.07, 0.08);
    let total = 1.0 / (3.0 * d);
    FewGroup {
        bounds: vec![1.0e7, 1.0e-5],
        // NEM は使わないが、少数群定数として正しい値にする
        k_inf: 1.0,
        total: vec![total],
        absorption: vec![absorption],
        nu_fission: vec![nu_fission],
        fission: vec![nu_fission / 2.4],
        chi: vec![1.0],
        diffusion: vec![d],
        scatter: vec![vec![total - absorption]],
        adf: vec![vec![1.0]; 6],
        vacuum_ratio: vec![None; 6],
        normalized_adf: vec![false; 6],
        height: None,
        pins: None,
    }
}

/// 1群の均質な体系の k（バックリング b2）
fn one_group_k(m: &FewGroup, b2: f64) -> f64 {
    m.nu_fission[0] / (m.absorption[0] + m.diffusion[0] * b2)
}

/// (0, π/x_max) にある f の最小の根（符号の変わる区間を探して二分法）
fn smallest_root(f: impl Fn(f64) -> f64, upper: f64) -> f64 {
    let steps = 10_000;
    let mut a = upper * 1e-6;
    for i in 1..=steps {
        let b = upper * i as f64 / steps as f64;
        if f(a).signum() != f(b).signum() {
            let (mut lo, mut hi) = (a, b);
            for _ in 0..200 {
                let mid = 0.5 * (lo + hi);
                if f(lo).signum() == f(mid).signum() {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            return 0.5 * (lo + hi);
        }
        a = b;
    }
    panic!("no root below {upper}");
}

/// x 方向の平板（y, z は反射）
fn slab(dx: Vec<f64>, materials: Vec<FewGroup>, node_material: Vec<usize>, x: Boundary) -> Problem {
    let mut boundaries = [Boundary::Reflective; 6];
    boundaries[0] = x;
    boundaries[1] = x;
    Problem {
        mesh: Mesh {
            dx,
            dy: vec![20.0],
            dz: vec![20.0],
        },
        materials,
        reflectors: Vec::new(),
        node_material: placed(node_material),
        boundaries,
        outline: None,
    }
}

/// x 方向に1列の問題（y と z は1ノード）を z 方向に並べ直した問題。x と z の幅、境界条件、少数群定数の
/// 辺（左、右）と（−z、+z）を入れ替える。物質の配置の並びはそのまま（ノードの番号は x でも z でも列の順）
fn to_z(p: &Problem) -> Problem {
    assert_eq!((p.mesh.dy.len(), p.mesh.dz.len()), (1, 1));
    fn swap<T>(sides: &mut [T]) {
        sides.swap(0, 4);
        sides.swap(1, 5);
    }
    let mut q = p.clone();
    q.mesh.dx = p.mesh.dz.clone();
    q.mesh.dz = p.mesh.dx.clone();
    swap(&mut q.boundaries);
    for m in &mut q.materials {
        swap(&mut m.adf);
        swap(&mut m.vacuum_ratio);
        swap(&mut m.normalized_adf);
    }
    q
}

/// x 方向に 6 ノードの平板。燃料と、断面積の違う燃料が交互に並ぶ
fn two_material_slab(adf_a: Vec<Vec<f64>>, adf_b: Vec<Vec<f64>>) -> Problem {
    let mut a = material();
    a.adf = adf_a;
    let mut b = material();
    b.diffusion = vec![1.2, 0.5];
    b.absorption = vec![0.012, 0.10];
    b.nu_fission = vec![0.004, 0.09];
    b.adf = adf_b;
    slab(
        vec![10.0; 6],
        vec![a, b],
        vec![0, 1, 0, 1, 0, 1],
        Boundary::ZeroFlux,
    )
}

/// 収束した状態も返す
fn solve_state(p: &Problem) -> (Core<'_>, Solution, State) {
    let mut core = Core::new(p).unwrap();
    let (s, state) = iterate(&mut core, &Options::default(), None).unwrap();
    (core, s, state)
}

/// 内側のすべての面で、不均質な中性子束 f φ（φ = 2(J+ + J-)）と正味の中性子流が連続であること。
/// 辺の番号は（左、右、上、下）で、-x は左、+x は右、-y は下、+y は上、z の面の ADF は 1
fn assert_interface_conditions(p: &Problem) {
    let (core, _, state) = solve_state(p);
    let adf = |n: usize, d: usize, plus: bool, g: usize| {
        let side = match (d, plus) {
            (0, false) => 0,
            (0, true) => 1,
            (1, true) => 2,
            (1, false) => 3,
            _ => return 1.0,
        };
        p.materials[p.node_material[n].unwrap()].adf[side][g]
    };
    let mut checked = 0;
    for n in 0..core.nodes() {
        for d in 0..3 {
            let Some(m) = core.neighbor(n, d, 1) else {
                continue;
            };
            for g in 0..core.groups {
                // n の + 面と m の - 面
                let (_, in_n) = core.incoming(&state, n, d, g);
                let (in_m, _) = core.incoming(&state, m, d, g);
                let (out_n, out_m) = (state.out_plus[d][n][g], state.out_minus[d][m][g]);
                let (phi_n, phi_m) = (2.0 * (out_n + in_n), 2.0 * (out_m + in_m));
                let (f_n, f_m) = (adf(n, d, true, g), adf(m, d, false, g));
                assert!(
                    (phi_n / phi_m - f_m / f_n).abs() < 1e-6 * f_m / f_n,
                    "node {n} dir {d} group {g}: {} vs {}",
                    phi_n / phi_m,
                    f_m / f_n
                );
                let (j_n, j_m) = (out_n - in_n, in_m - out_m);
                assert!(
                    (j_n - j_m).abs() < 1e-6 * phi_n,
                    "node {n} dir {d} group {g}: {j_n} vs {j_m}"
                );
                checked += 1;
            }
        }
    }
    assert!(checked > 0);
}

#[test]
fn unit_adfs_leave_results_unchanged() {
    // ADF を入れる前（e3b627c）と同じ反復（a3, a4 を前の反復の a1, a2 から求める。#0025 の修正の前で、
    // 既定の許容差では e3b627c の結果に 1e-12 で一致していた）を、許容差を詰めて収束させた値。
    // 既定の許容差の値は反復の道筋で 1e-9 程度変わるので、収束させた値で比べる
    let p = two_material_slab(vec![vec![1.0, 1.0]; 6], vec![vec![1.0, 1.0]; 6]);
    let tight = Options {
        k_tolerance: 1e-14,
        flux_tolerance: 1e-13,
        max_iterations: 10_000,
    };
    let s = solve(&p, &tight).unwrap();
    let k = 8.510692798259939e-1;
    assert!((s.k - k).abs() < 1e-12, "{} vs {k}", s.k);
    let flux = [
        [0.0005911585762519943, 0.00012519241937033706],
        [0.001661405206417996, 0.0003510047208700676],
    ];
    for (phi, expected) in s.flux.iter().zip(flux) {
        for (a, b) in phi.iter().zip(expected) {
            assert!((a - b).abs() < 1e-10 * b, "{a} vs {b}");
        }
    }
    assert_interface_conditions(&p);
}

#[test]
fn adfs_make_homogeneous_flux_jump_at_faces() {
    // 群と辺で違う値にして、辺の対応も確かめる
    let a = vec![
        vec![1.05, 1.3],
        vec![1.1, 1.25],
        vec![0.9, 1.4],
        vec![0.95, 1.2],
        vec![1.0, 1.0],
        vec![1.0, 1.0],
    ];
    let b = vec![
        vec![0.97, 0.8],
        vec![0.93, 0.85],
        vec![1.02, 0.75],
        vec![0.99, 0.9],
        vec![1.0, 1.0],
        vec![1.0, 1.0],
    ];
    let p = two_material_slab(a.clone(), b.clone());
    assert_interface_conditions(&p);
    let unit = two_material_slab(vec![vec![1.0, 1.0]; 6], vec![vec![1.0, 1.0]; 6]);
    let (k, k_unit) = (
        solve(&p, &Options::default()).unwrap().k,
        solve(&unit, &Options::default()).unwrap().k,
    );
    assert!((k - k_unit).abs() > 1e-4, "{k} vs {k_unit}");

    // x と y に 3 × 2 の市松模様、z に 2 層
    let mut q = cuboid([30.0, 20.0, 40.0], [3, 2, 2], [Boundary::ZeroFlux; 6]);
    let mut m = material();
    m.adf = b;
    m.absorption[1] = 0.09;
    q.materials[0].adf = a;
    q.materials.push(m);
    q.node_material = placed((0..12).map(|n| (n + n / 3) % 2));
    assert_interface_conditions(&q);
}

#[test]
fn two_region_slab_matches_analytic_jump() {
    // 同じ断面積の領域 A (0, a) と B (a, a + b) の1群の平板、両端は中性子束 0。
    // φ_A = sin(Bx)、φ_B = c sin(B(a + b - x)) を f_A φ_A = f_B φ_B と中性子流の連続に入れると
    // f_A sin(Ba) cos(Bb) + f_B cos(Ba) sin(Bb) = 0
    let (a, b) = (30.0, 50.0);
    let (f_a, f_b) = (1.0, 1.5);
    let buckling = smallest_root(
        |x| f_a * (x * a).sin() * (x * b).cos() + f_b * (x * a).cos() * (x * b).sin(),
        PI / a,
    );
    let mut m_a = one_group();
    m_a.adf = vec![vec![f_a]; 6];
    let mut m_b = one_group();
    m_b.adf = vec![vec![f_b]; 6];
    let expected = one_group_k(&m_a, buckling * buckling);
    let p = slab(
        vec![5.0; 16],
        vec![m_a, m_b],
        (0..16).map(|i| usize::from(i >= 6)).collect(),
        Boundary::ZeroFlux,
    );
    let s = solve(&p, &Options::default()).unwrap();
    assert!(
        (s.k - expected).abs() < 1e-8 * expected,
        "{} vs {expected}",
        s.k
    );
    // ADF がなければ k は B = π/(a + b) のもの
    let plain = one_group_k(&one_group(), (PI / (a + b)).powi(2));
    assert!((expected - plain).abs() > 1e-3 * plain);
}

#[test]
fn uniform_adf_scaling_keeps_k_without_vacuum() {
    // 外周が反射と中性子束 0 だけなら、面の条件は ADF の比だけで決まる（真空の面では変わる。
    // outer_adf_at_vacuum_matches_analytic）
    let a = vec![
        vec![1.05, 1.3],
        vec![1.1, 1.25],
        vec![0.9, 1.4],
        vec![0.95, 1.2],
        vec![1.0, 1.0],
        vec![1.0, 1.0],
    ];
    let b = vec![
        vec![0.97, 0.8],
        vec![0.93, 0.85],
        vec![1.02, 0.75],
        vec![0.99, 0.9],
        vec![1.0, 1.0],
        vec![1.0, 1.0],
    ];
    // 群ごとに違う倍率で、すべての物質とすべての辺の ADF を変える
    let scale = |adf: &[Vec<f64>]| -> Vec<Vec<f64>> {
        adf.iter().map(|s| vec![2.0 * s[0], 0.5 * s[1]]).collect()
    };
    let base = solve(
        &two_material_slab(a.clone(), b.clone()),
        &Options::default(),
    )
    .unwrap();
    let scaled = solve(
        &two_material_slab(scale(&a), scale(&b)),
        &Options::default(),
    )
    .unwrap();
    assert!(
        (base.k - scaled.k).abs() < 1e-9,
        "{} vs {}",
        base.k,
        scaled.k
    );
    for (x, y) in base.flux.iter().flatten().zip(scaled.flux.iter().flatten()) {
        assert!((x - y).abs() < 1e-6 * x, "{x} vs {y}");
    }
}

#[test]
fn invalid_adfs_are_errors() {
    let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
    p.materials[0].adf.pop();
    assert!(solve(&p, &Options::default()).is_err());
    let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
    p.materials[0].adf[1].pop();
    assert!(solve(&p, &Options::default()).is_err());
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], [Boundary::Reflective; 6]);
        p.materials[0].adf[2][1] = bad;
        let err = solve(&p, &Options::default()).unwrap_err().to_string();
        assert!(err.contains("discontinuity"), "{err}");
    }
}

#[test]
fn vacuum_face_has_no_incoming_heterogeneous_current() {
    // 不均質な中性子束 f φ = 2f(o + i) と正味の中性子流 o - i から、入射の部分中性子流は
    // f φ/4 - (o - i)/2 = (f(o + i) - (o - i))/2 で、真空の面ではこれが 0
    let mut m = material();
    m.adf[0] = vec![1.3, 0.8];
    m.adf[1] = vec![0.7, 1.6];
    let p = slab(vec![5.0; 6], vec![m.clone()], vec![0; 6], Boundary::Marshak);
    let (core, s, state) = solve_state(&p);
    for g in 0..2 {
        let (in_minus, _) = core.incoming(&state, 0, 0, g);
        let (_, in_plus) = core.incoming(&state, 5, 0, g);
        let (out_minus, out_plus) = (state.out_minus[0][0][g], state.out_plus[0][5][g]);
        for (f, o, i) in [
            (m.adf[0][g], out_minus, in_minus),
            (m.adf[1][g], out_plus, in_plus),
        ] {
            assert!(o > 0.0);
            let incoming = (f * (o + i) - (o - i)) / 2.0;
            assert!(incoming.abs() < 1e-14 * o, "{incoming} for {o}");
        }
        // 外へ出る正味の中性子数は、外周の面の正味の中性子流の和
        let area = 20.0 * 20.0;
        let out = area * (out_minus - in_minus + out_plus - in_plus);
        assert!((s.leakage[g] - out).abs() < 1e-12 * out);
    }
}

#[test]
fn vacuum_k_lies_between_zero_flux_and_reflective() {
    let k = |x: Boundary| {
        let p = slab(vec![5.0; 6], vec![material()], vec![0; 6], x);
        solve(&p, &Options::default()).unwrap().k
    };
    let (zero, vacuum, reflective) = (
        k(Boundary::ZeroFlux),
        k(Boundary::Marshak),
        k(Boundary::Reflective),
    );
    assert!(
        zero < vacuum && vacuum < reflective,
        "{zero} {vacuum} {reflective}"
    );
}

/// 1群の裸の平板（幅 l、外周の面の ADF f）に、不均質な中性子束 f φ で Marshak の条件
/// J- = f φ/4 + (D/2) dφ/dx = 0 を課した解 φ = cos(Bx) は tan(Bl/2) = f/(2DB) を満たす
/// （外挿距離は 2D/f）
fn marshak_buckling(m: &FewGroup, l: f64, f: f64) -> f64 {
    let d = m.diffusion[0];
    smallest_root(
        |b| (b * l / 2.0).sin() * 2.0 * d * b - f * (b * l / 2.0).cos(),
        PI / l,
    )
}

#[test]
fn outer_adf_at_vacuum_matches_analytic() {
    // 1種類の物質なので内側の面には跳びがなく、外周の真空の面の ADF だけが効く。
    // ADF をすべて同じ倍率で変えると、k は解析解のとおり変わる
    let l = 40.0;
    // f が大きいほど外挿距離 2D/f が短く、漏れが増えて k は下がる
    let mut previous = f64::INFINITY;
    for f in [0.7, 1.0, 1.5, 2.0] {
        let mut m = one_group();
        m.adf = vec![vec![f]; 6];
        let b = marshak_buckling(&m, l, f);
        let expected = one_group_k(&m, b * b);
        let p = slab(vec![l / 8.0; 8], vec![m], vec![0; 8], Boundary::Marshak);
        let s = solve(&p, &Options::default()).unwrap();
        assert!(
            (s.k - expected).abs() < 1e-8 * expected,
            "f = {f}: {} vs {expected}",
            s.k
        );
        assert!(expected < previous, "f = {f}");
        previous = expected;
    }
}

#[test]
fn bare_slab_matches_marshak_analytic() {
    let m = one_group();
    let l = 40.0;
    let b = marshak_buckling(&m, l, 1.0);
    let expected = one_group_k(&m, b * b);
    let p = slab(vec![l / 8.0; 8], vec![m], vec![0; 8], Boundary::Marshak);
    let s = solve(&p, &Options::default()).unwrap();
    assert!(
        (s.k - expected).abs() < 1e-8 * expected,
        "{} vs {expected}",
        s.k
    );
    // 外挿距離 2D を無視した中性子束 0 の解とは区別できる
    let zero = one_group_k(&one_group(), (PI / l).powi(2));
    assert!((expected - zero).abs() > 1e-3 * zero);
}

#[test]
fn bare_cuboid_matches_marshak_analytic() {
    // 各方向の解は分離できるので B² は各方向の Marshak の解の和。横方向の漏れの近似を通る
    let m = one_group();
    let size = [40.0, 50.0, 60.0];
    let b2: f64 = size
        .iter()
        .map(|&l| marshak_buckling(&m, l, 1.0).powi(2))
        .sum();
    let expected = one_group_k(&m, b2);
    let error = |n: usize| {
        let mut p = cuboid(size, [n; 3], [Boundary::Marshak; 6]);
        p.materials = vec![one_group()];
        (solve(&p, &Options::default()).unwrap().k - expected).abs() / expected
    };
    let (coarse, fine) = (error(4), error(8));
    assert!(fine < 1e-4, "relative error {fine}");
    assert!(fine < coarse / 8.0, "{coarse} -> {fine}");
}

/// 外周の x の2つの面（辺の左と右）に、真空の面の比 r を群ごとに持たせる
fn with_x_ratios(mut m: FewGroup, r: &[f64]) -> FewGroup {
    m.vacuum_ratio[0] = Some(r.to_vec());
    m.vacuum_ratio[1] = Some(r.to_vec());
    m
}

#[test]
fn half_ratio_vacuum_equals_marshak() {
    // r = 1/2 の真空の面の条件 J = r f φ は、Marshak の条件 f φ/4 = J/2 と同じ（ADR-0128）
    let mut m = material();
    m.adf[0] = vec![1.3, 0.8];
    m.adf[1] = vec![0.7, 1.6];
    let marshak = solve(
        &slab(vec![5.0; 6], vec![m.clone()], vec![0; 6], Boundary::Marshak),
        &Options::default(),
    )
    .unwrap();
    let p = slab(
        vec![5.0; 6],
        vec![with_x_ratios(m, &[0.5, 0.5])],
        vec![0; 6],
        Boundary::Vacuum,
    );
    let vacuum = solve(&p, &Options::default()).unwrap();
    assert!(
        (vacuum.k - marshak.k).abs() < 1e-14,
        "{} vs {}",
        vacuum.k,
        marshak.k
    );
    for (a, b) in vacuum
        .flux
        .iter()
        .flatten()
        .zip(marshak.flux.iter().flatten())
    {
        assert!((a - b).abs() < 1e-12 * b, "{a} vs {b}");
    }
}

#[test]
fn vacuum_face_imposes_the_ratio() {
    // 真空の面では、不均質な中性子束 f φ = 2f(o + i) と正味の中性子流 o - i が J = r f φ を満たす
    let mut m = material();
    m.adf[0] = vec![1.3, 0.8];
    m.adf[1] = vec![0.7, 1.6];
    let r = [0.62, 0.55];
    let p = slab(
        vec![5.0; 6],
        vec![with_x_ratios(m.clone(), &r)],
        vec![0; 6],
        Boundary::Vacuum,
    );
    let (core, _, state) = solve_state(&p);
    for (g, ratio) in r.into_iter().enumerate() {
        let (in_minus, _) = core.incoming(&state, 0, 0, g);
        let (_, in_plus) = core.incoming(&state, 5, 0, g);
        let (out_minus, out_plus) = (state.out_minus[0][0][g], state.out_plus[0][5][g]);
        for (f, o, i) in [
            (m.adf[0][g], out_minus, in_minus),
            (m.adf[1][g], out_plus, in_plus),
        ] {
            let (current, heterogeneous) = (o - i, 2.0 * f * (o + i));
            assert!(current > 0.0);
            assert!(
                (current - ratio * heterogeneous).abs() < 1e-14 * current,
                "{current} vs {}",
                ratio * heterogeneous
            );
        }
    }
}

#[test]
fn vacuum_ratio_matches_analytic() {
    // 1群の裸の平板に J = r f φ を課した解 φ = cos(Bx) は D B tan(Bl/2) = r f を満たす。
    // r が大きいほど漏れが増えて k は下がる
    let l = 40.0;
    for f in [0.8, 1.3] {
        let mut previous = f64::INFINITY;
        for r in [0.5, 0.55, 0.6, 0.9] {
            let mut m = one_group();
            m.adf = vec![vec![f]; 6];
            let d = m.diffusion[0];
            let b = smallest_root(
                |b| (b * l / 2.0).sin() * d * b - r * f * (b * l / 2.0).cos(),
                PI / l,
            );
            let expected = one_group_k(&m, b * b);
            let p = slab(
                vec![l / 8.0; 8],
                vec![with_x_ratios(m, &[r])],
                vec![0; 8],
                Boundary::Vacuum,
            );
            // z の面も同じ規則（ADR-0128）
            for p in [p.clone(), to_z(&p)] {
                let s = solve(&p, &Options::default()).unwrap();
                assert!(
                    (s.k - expected).abs() < 1e-8 * expected,
                    "f = {f}, r = {r}: {} vs {expected}",
                    s.k
                );
            }
            assert!(expected < previous, "f = {f}, r = {r}");
            previous = expected;
        }
    }
}

#[test]
fn vacuum_needs_ratios() {
    // 物質がその面の辺の比を持たなければエラー。黙って Marshak の条件にしない
    let p = slab(vec![5.0; 6], vec![material()], vec![0; 6], Boundary::Vacuum);
    let err = solve(&p, &Options::default()).unwrap_err().to_string();
    assert!(err.contains("ratio") && err.contains("marshak"), "{err}");
    // 片側だけ持つ
    let mut m = material();
    m.vacuum_ratio[0] = Some(vec![0.5, 0.5]);
    let p = slab(vec![5.0; 6], vec![m], vec![0; 6], Boundary::Vacuum);
    assert!(solve(&p, &Options::default()).is_err());
    // z の面も、その面の辺（+z）の比を使う（ADR-0128）。径方向の辺がすべて比を持っていても、+z の辺が
    // 持たなければエラー
    let mut m = material();
    m.vacuum_ratio = (0..6).map(|s| (s < 4).then(|| vec![0.5, 0.5])).collect();
    let mut boundaries = [Boundary::Reflective; 6];
    boundaries[5] = Boundary::Vacuum;
    let mut p = cuboid([30.0, 30.0, 30.0], [2, 2, 2], boundaries);
    p.materials = vec![m];
    let err = solve(&p, &Options::default()).unwrap_err().to_string();
    assert!(err.contains("side 5") && err.contains("marshak"), "{err}");
    p.materials[0].vacuum_ratio[5] = Some(vec![0.5, 0.5]);
    assert!(solve(&p, &Options::default()).is_ok());
}

#[test]
fn side_mapping_matches_two_region_analytic() {
    // 2領域の平板を x、y、z の各方向に置き、ADF の跳びを解析解（two_region_slab_matches_analytic_jump）
    // と比べる。跳びを境の2つのノードだけに置く: 境の A のノードは - 面で 1、+ 面で f_A、境の B の
    // ノードは - 面で f_B、+ 面で 1、ほかのノードと辺は 1 か関係のない値（その方向に隣がなく、
    // 反射の面では ADF は効かない）。- 面と + 面の辺を取り違えれば、跳びが境からずれて解析解から外れる
    let (a, b) = (30.0, 50.0);
    let (f_a, f_b) = (1.0, 1.5);
    let analytic = |f_a: f64, f_b: f64| {
        let buckling = smallest_root(
            |x| f_a * (x * a).sin() * (x * b).cos() + f_b * (x * a).cos() * (x * b).sin(),
            PI / a,
        );
        one_group_k(&one_group(), buckling * buckling)
    };
    // 方向ごとの（- 面, + 面）の辺（左、右、上、下、−z、+z の番号）。z の面も少数群定数の辺の ADF を使う
    // （ADR-0128）
    let expected = analytic(f_a, f_b);
    for (d, (minus, plus)) in [(0, 1), (3, 2), (4, 5)].into_iter().enumerate() {
        let bulk = one_group();
        let mut edge_a = one_group();
        let mut edge_b = one_group();
        edge_a.adf = vec![vec![3.0]; 6];
        edge_a.adf[minus] = vec![1.0];
        edge_a.adf[plus] = vec![f_a];
        edge_b.adf = vec![vec![0.5]; 6];
        edge_b.adf[minus] = vec![f_b];
        edge_b.adf[plus] = vec![1.0];
        let mut widths = [vec![20.0], vec![20.0], vec![20.0]];
        widths[d] = vec![5.0; 16];
        let mut boundaries = [Boundary::Reflective; 6];
        boundaries[2 * d] = Boundary::ZeroFlux;
        boundaries[2 * d + 1] = Boundary::ZeroFlux;
        let [dx, dy, dz] = widths;
        let p = Problem {
            mesh: Mesh { dx, dy, dz },
            materials: vec![bulk, edge_a, edge_b],
            reflectors: Vec::new(),
            node_material: placed((0..16).map(|i| match i {
                5 => 1,
                6 => 2,
                _ => 0,
            })),
            boundaries,
            outline: None,
        };
        let k = solve(&p, &Options::default()).unwrap().k;
        assert!(
            (k - expected).abs() < 1e-8 * expected,
            "direction {d}: {k} vs {expected}"
        );
    }
}

#[test]
fn vacuum_corner_converges_with_refinement() {
    // 40 cm 角の 2D（-x, -y は反射、+x, +y は真空）。20 cm 角の4区画に2種類の物質を市松に置き、
    // 真空の角にも物質の境がある分離できない問題で、細かいノードの解に近づくこと。
    // 横方向の漏れの真空の面での扱いと、物質の境の ADF を通る
    let mut a = material();
    a.adf = vec![vec![1.1, 1.3]; 6];
    let mut b = material();
    b.diffusion = vec![1.2, 0.5];
    b.absorption = vec![0.012, 0.10];
    b.nu_fission = vec![0.004, 0.09];
    b.adf = vec![vec![0.95, 0.8]; 6];
    let k = |n: usize| {
        let mut boundaries = [Boundary::Reflective; 6];
        boundaries[1] = Boundary::Marshak;
        boundaries[3] = Boundary::Marshak;
        let p = Problem {
            mesh: Mesh {
                dx: vec![40.0 / n as f64; n],
                dy: vec![40.0 / n as f64; n],
                dz: vec![20.0],
            },
            materials: vec![a.clone(), b.clone()],
            reflectors: Vec::new(),
            node_material: placed(
                (0..n * n).map(|i| usize::from((i % n < n / 2) != (i / n < n / 2))),
            ),
            boundaries,
            outline: None,
        };
        solve(&p, &Options::default()).unwrap().k
    };
    let reference = k(32);
    let error = |n: usize| (k(n) - reference).abs() / reference;
    let (coarse, fine) = (error(4), error(8));
    // 測った値は 9.6e-5 -> 1.1e-5
    assert!(coarse < 2e-4, "relative error {coarse}");
    assert!(fine < coarse / 6.0, "{coarse} -> {fine}");
}

#[test]
fn no_flux_next_to_a_vacuum_face_is_error() {
    // 真空の面の横方向の漏れは面の中性子束 / ノード平均を使うので、ノード平均が正でなければ止める
    let p = slab(
        vec![5.0; 6],
        vec![material()],
        vec![0; 6],
        Boundary::Marshak,
    );
    let core = Core::new(&p).unwrap();
    for bad in [0.0, -1.0, f64::NAN] {
        let mut state = State::new(6, 2);
        state.flux[0][1] = bad;
        let err = core
            .sweep_group(&mut state, 1, 1.0)
            .unwrap_err()
            .to_string();
        assert!(err.contains("vacuum"), "{err}");
    }
    // 真空の面がなければ問題にしない
    let p = slab(
        vec![5.0; 6],
        vec![material()],
        vec![0; 6],
        Boundary::ZeroFlux,
    );
    let core = Core::new(&p).unwrap();
    let mut state = State::new(6, 2);
    state.flux[0][1] = 0.0;
    assert!(core.sweep_group(&mut state, 1, 1.0).is_ok());
}

/// 細かいノードの解を参照（不均質な解）とし、粗いノード（x, y に `blocks` 個ずつの細かいノードの
/// まとまり）ごとに反応率を保存して定数を均質化し、GET の1ノードの問題（ADR-0103）で ADF を作った
/// 粗いノードの問題と、細かいノードの解を返す。z は1層。
///
/// 1ノードの問題の外周の扱い（`Edge`）は、粗いノードの面が外周なら細かいノードの問題の境界、内側の面
/// なら隣の粗いノード（幅と参照の横方向の漏れ。ADR-0095）とする
fn coarsen_with_get(fine: &Problem, blocks: [usize; 2]) -> (Problem, Solution) {
    use xlc_numeric::nem::{Edge, OneNode, surface_fluxes};
    let (core, solution, state) = solve_state(fine);
    let core = &core;
    let [nx, ny, nz] = core.problem.mesh.counts();
    assert_eq!(nz, 1);
    let [bx, by] = blocks;
    let (cx, cy) = (nx / bx, ny / by);
    assert_eq!((cx * bx, cy * by), (nx, ny));
    let ng = core.groups;
    // [粗いノード] の定数、ノード平均、面（-x, +x, -y, +y）の不均質な中性子束と外向きの正味の中性子流、幅
    let mut nodes = Vec::new();
    for (cj, ci) in (0..cy).flat_map(|j| (0..cx).map(move |i| (j, i))) {
        // どの位置にもノードがあるので、ノードの番号は位置の番号
        let members: Vec<usize> = (cj * by..(cj + 1) * by)
            .flat_map(|j| (ci * bx..(ci + 1) * bx).map(move |i| i + nx * j))
            .collect();
        // 反応率を保存する均質化（D は中性子束の重みの平均。GET の ADF が差を吸収する）
        let mut volume = 0.0;
        let mut flux = vec![0.0; ng];
        let mut m = FewGroup {
            bounds: fine.materials[0].bounds.clone(),
            k_inf: solution.k,
            total: vec![0.0; ng],
            absorption: vec![0.0; ng],
            nu_fission: vec![0.0; ng],
            fission: vec![0.0; ng],
            chi: vec![0.0; ng],
            diffusion: vec![0.0; ng],
            scatter: vec![vec![0.0; ng]; ng],
            adf: vec![vec![1.0; ng]; 6],
            vacuum_ratio: vec![None; 6],
            normalized_adf: vec![false; 6],
            height: None,
            pins: None,
        };
        for &n in &members {
            let (v, f, x) = (core.volume(n), &solution.flux[n], core.material(n));
            volume += v;
            let production: f64 = (0..ng).map(|g| x.nu_fission[g] * f[g]).sum();
            for g in 0..ng {
                let w = v * f[g];
                flux[g] += w;
                m.total[g] += w * x.total[g];
                m.nu_fission[g] += w * x.nu_fission[g];
                m.fission[g] += w * x.fission[g];
                m.diffusion[g] += w * x.diffusion[g];
                m.chi[g] += v * x.chi[g] * production;
                for to in 0..ng {
                    m.scatter[g][to] += w * x.scatter[g][to];
                }
            }
        }
        for (g, f) in flux.iter().enumerate() {
            for value in [
                &mut m.total[g],
                &mut m.nu_fission[g],
                &mut m.fission[g],
                &mut m.diffusion[g],
            ] {
                *value /= f;
            }
            for value in &mut m.scatter[g] {
                *value /= f;
            }
            m.absorption[g] = m.total[g] - m.scatter[g].iter().sum::<f64>();
        }
        let produced: f64 = m.chi.iter().sum();
        if produced > 0.0 {
            for c in &mut m.chi {
                *c /= produced;
            }
        }
        let average: Vec<f64> = flux.iter().map(|f| f / volume).collect();

        // 面（-x, +x, -y, +y）の不均質な中性子束と外向きの正味の中性子流の、面にわたる平均
        let mut het = [vec![0.0; ng], vec![0.0; ng], vec![0.0; ng], vec![0.0; ng]];
        let mut current = het.clone();
        for (face, d, plus) in [(0, 0, false), (1, 0, true), (2, 1, false), (3, 1, true)] {
            let on_face: Vec<usize> = members
                .iter()
                .copied()
                .filter(|&n| {
                    let p = core.position(n);
                    let lowest = if d == 0 { ci * bx } else { cj * by };
                    let count = if d == 0 { bx } else { by };
                    p[d] == if plus { lowest + count - 1 } else { lowest }
                })
                .collect();
            let length: f64 = on_face
                .iter()
                .map(|&n| core.volume(n) / core.width(n, d))
                .sum();
            for &n in &on_face {
                let w = core.volume(n) / core.width(n, d) / length;
                for g in 0..ng {
                    let (in_minus, in_plus) = core.incoming(&state, n, d, g);
                    let (o, i) = if plus {
                        (state.out_plus[d][n][g], in_plus)
                    } else {
                        (state.out_minus[d][n][g], in_minus)
                    };
                    het[face][g] += w * 2.0 * (o + i);
                    current[face][g] += w * (o - i);
                }
            }
        }
        let widths: [f64; 2] = [
            core.widths[0][ci * bx..(ci + 1) * bx].iter().sum(),
            core.widths[1][cj * by..(cj + 1) * by].iter().sum(),
        ];
        nodes.push((m, average, het, current, widths));
    }
    // [粗いノード][方向] の横方向の漏れ（もう一方の方向の正味の漏れ）
    let transverse: Vec<[Vec<f64>; 2]> = nodes
        .iter()
        .map(|(_, _, _, current, widths)| {
            let net = |d: usize| -> Vec<f64> {
                (0..ng)
                    .map(|g| (current[2 * d][g] + current[2 * d + 1][g]) / widths[d])
                    .collect()
            };
            [net(1), net(0)]
        })
        .collect();
    let mut materials = Vec::new();
    for (c, (mut m, average, het, current, widths)) in nodes.clone().into_iter().enumerate() {
        let (ci, cj) = (c % cx, c / cx);
        let edge = |d: usize, plus: bool| {
            let neighbor = match (d, plus) {
                (0, false) => (ci > 0).then(|| c - 1),
                (0, true) => (ci + 1 < cx).then_some(c + 1),
                (_, false) => (cj > 0).then(|| c - cx),
                (_, true) => (cj + 1 < cy).then_some(c + cx),
            };
            match (neighbor, fine.boundaries[2 * d + usize::from(plus)]) {
                (Some(other), _) => Edge::Neighbor {
                    width: nodes[other].4[d],
                    leakage: &transverse[other][d],
                },
                (None, Boundary::Marshak) => Edge::Vacuum,
                (None, Boundary::Reflective) => Edge::Reflective,
                (None, Boundary::ZeroFlux) => panic!("no zero-flux edge in GET"),
                (None, Boundary::Vacuum) => panic!("the fine reference uses marshak"),
            }
        };
        let surface = surface_fluxes(&OneNode {
            diffusion: &m.diffusion,
            total: &m.total,
            scatter: &m.scatter,
            nu_fission: &m.nu_fission,
            chi: &m.chi,
            k: solution.k,
            flux: &average,
            widths,
            currents: [&current[0], &current[1], &current[2], &current[3]],
            edges: [edge(0, false), edge(0, true), edge(1, false), edge(1, true)],
        })
        .unwrap();
        // 少数群定数の辺は（左、右、上、下）で、-x が左、+x が右、-y が下、+y が上
        for (side, face) in [(0, 0), (1, 1), (2, 3), (3, 2)] {
            m.adf[side] = (0..ng).map(|g| het[face][g] / surface[face][g]).collect();
            // 参照の真空の面の比 J/φ_het（ADR-0128。参照は Marshak の条件なので 1/2）
            let (d, plus) = (face / 2, face % 2 == 1);
            if edge(d, plus) == Edge::Vacuum {
                m.vacuum_ratio[side] =
                    Some((0..ng).map(|g| current[face][g] / het[face][g]).collect());
            }
        }
        materials.push(m);
    }
    let merge = |widths: &[f64], b: usize| -> Vec<f64> {
        widths.chunks(b).map(|c| c.iter().sum()).collect()
    };
    let coarse = Problem {
        mesh: Mesh {
            dx: merge(&fine.mesh.dx, bx),
            dy: merge(&fine.mesh.dy, by),
            dz: fine.mesh.dz.clone(),
        },
        reflectors: Vec::new(),
        node_material: placed(0..materials.len()),
        materials,
        // 粗いノードの真空の面は、参照から作った比 J/φ_het を使う
        boundaries: fine.boundaries.map(|b| match b {
            Boundary::Marshak => Boundary::Vacuum,
            other => other,
        }),
        outline: None,
    };
    (coarse, solution)
}

/// 粗いノードの問題の k とノード平均の中性子束が、細かいノードの解（同じ規格化: 全体の生成率 1）を
/// 再現することを確かめ、ADF をすべて 1 にしたときの k の相対差を返す
fn assert_get_equivalence(coarse: &Problem, fine: &Solution, blocks: [usize; 2], nx: usize) {
    let s = solve(coarse, &Options::default()).unwrap();
    assert!(
        (s.k - fine.k).abs() < 1e-7 * fine.k,
        "{} vs {}",
        s.k,
        fine.k
    );
    let [bx, by] = blocks;
    let cx = nx / bx;
    for (c, phi) in s.flux.iter().enumerate() {
        let (ci, cj) = (c % cx, c / cx);
        // どの位置にもノードがあるので、ノードの番号は位置の番号
        let members: Vec<usize> = (cj * by..(cj + 1) * by)
            .flat_map(|j| (ci * bx..(ci + 1) * bx).map(move |i| i + nx * j))
            .collect();
        for (g, value) in phi.iter().enumerate() {
            // 細かいノードはどれも同じ大きさ
            let expected =
                members.iter().map(|&n| fine.flux[n][g]).sum::<f64>() / members.len() as f64;
            assert!(
                (value - expected).abs() < 1e-6 * expected,
                "node {c} group {g}: {value} vs {expected}"
            );
        }
    }
    let mut unit = coarse.clone();
    for m in &mut unit.materials {
        m.adf = vec![vec![1.0; m.groups()]; 6];
    }
    let k_unit = solve(&unit, &Options::default()).unwrap().k;
    // ADF がなければ、粗いノードの NEM の誤差の分だけ k が外れる
    assert!(
        (k_unit - fine.k).abs() > 1e-4 * fine.k,
        "{k_unit} vs {}",
        fine.k
    );
}

/// 細かいノードの物質: 燃料、断面積の違う燃料、核分裂の少ない燃料
fn three_fuels() -> Vec<FewGroup> {
    let a = material();
    let mut b = material();
    b.diffusion = vec![1.2, 0.5];
    b.absorption = vec![0.012, 0.10];
    b.nu_fission = vec![0.004, 0.09];
    let mut c = material();
    c.absorption = vec![0.009, 0.06];
    c.nu_fission = vec![0.002, 0.04];
    c.fission = vec![0.002 / 2.4, 0.04 / 2.4];
    vec![a, b, c]
}

#[test]
fn get_adfs_reproduce_a_two_node_slab() {
    // x 方向に 20 cm の粗いノード2つ（細かいノード 2.5 cm × 8 ずつ、中に2種類の物質）。-x は反射、
    // +x は真空。粗いノードの間の面で、ADF をかけた不均質な中性子束と正味の中性子流がつながる
    let mut boundaries = [Boundary::Reflective; 6];
    boundaries[1] = Boundary::Marshak;
    let fine = Problem {
        mesh: Mesh {
            dx: vec![2.5; 16],
            dy: vec![20.0],
            dz: vec![20.0],
        },
        materials: three_fuels(),
        reflectors: Vec::new(),
        node_material: placed(vec![0, 1, 0, 1, 1, 0, 1, 0, 2, 2, 1, 1, 2, 2, 0, 0]),
        boundaries,
        outline: None,
    };
    let (coarse, solution) = coarsen_with_get(&fine, [8, 1]);
    assert_get_equivalence(&coarse, &solution, [8, 1], 16);
}

#[test]
fn get_adfs_reproduce_a_node_with_vacuum_sides() {
    // 20 cm 角の粗いノード1つ（細かいノード 5 cm 角 × 16、3種類の物質）。-x と -y は反射、+x と +y は
    // 真空。真空の面での横方向の漏れの扱いと、2方向の横方向の漏れを通る
    let mut boundaries = [Boundary::Reflective; 6];
    boundaries[1] = Boundary::Marshak;
    boundaries[3] = Boundary::Marshak;
    let fine = Problem {
        mesh: Mesh {
            dx: vec![5.0; 4],
            dy: vec![5.0; 4],
            dz: vec![20.0],
        },
        materials: three_fuels(),
        reflectors: Vec::new(),
        node_material: placed(vec![0, 1, 0, 2, 1, 0, 2, 0, 0, 2, 1, 1, 2, 0, 1, 0]),
        boundaries,
        outline: None,
    };
    let (coarse, solution) = coarsen_with_get(&fine, [4, 4]);
    assert_get_equivalence(&coarse, &solution, [4, 4], 4);
}

/// 核分裂のない2群の物質（水に近い。上方散乱あり）
fn reflector() -> FewGroup {
    let mut m = material();
    m.diffusion = vec![1.25, 0.30];
    m.total = vec![0.60, 1.90];
    m.scatter = vec![vec![0.565, 0.032], vec![0.0015, 1.885]];
    m.absorption = vec![0.003, 0.0135];
    m.nu_fission = vec![0.0; 2];
    m.fission = vec![0.0; 2];
    m.chi = vec![0.0; 2];
    m
}

#[test]
fn get_adfs_reproduce_fuel_and_a_wide_reflector() {
    // 燃料（20 cm）と反射体（21.5 cm。細かいノード 2.5 cm × 8 と 2.6875 cm × 8）の 1D、+x は真空。
    // 粗い反射体のノードは core と1ノードの問題の両方で ANM で解き（ADR-0102）、GET の定数で細かいノードの
    // 解を再現する
    let mut boundaries = [Boundary::Reflective; 6];
    boundaries[1] = Boundary::Marshak;
    let mut materials = three_fuels();
    materials.push(reflector());
    let mut dx = vec![2.5; 8];
    dx.extend(vec![21.5 / 8.0; 8]);
    let fine = Problem {
        mesh: Mesh {
            dx,
            dy: vec![20.0],
            dz: vec![20.0],
        },
        materials,
        reflectors: Vec::new(),
        node_material: placed(vec![0, 1, 0, 1, 1, 0, 2, 0, 3, 3, 3, 3, 3, 3, 3, 3]),
        boundaries,
        outline: None,
    };
    let (coarse, solution) = coarsen_with_get(&fine, [8, 1]);
    assert!(coarse.materials[1].nu_fission.iter().all(|&v| v == 0.0));
    assert_get_equivalence(&coarse, &solution, [8, 1], 16);
}

#[test]
fn non_fissile_nodes_need_two_groups() {
    // ANM の核は2群（ADR-0102）
    let mut m = reflector();
    m.bounds = vec![1.0e7, 1.0e-5];
    m.total = vec![0.6];
    m.absorption = vec![0.01];
    m.diffusion = vec![1.0];
    m.scatter = vec![vec![0.59]];
    m.nu_fission = vec![0.0];
    m.fission = vec![0.0];
    m.chi = vec![0.0];
    m.adf = vec![vec![1.0]; 6];
    let mut fuel = m.clone();
    fuel.nu_fission = vec![0.02];
    fuel.fission = vec![0.008];
    fuel.chi = vec![1.0];
    let p = slab(
        vec![10.0; 2],
        vec![fuel, m],
        vec![0, 1],
        Boundary::Reflective,
    );
    let err = solve(&p, &Options::default());
    assert!(matches!(err, Err(NemError::InvalidProblem(_))), "{err:?}");
}

#[test]
fn negative_nu_fission_is_an_error() {
    // νΣf で核の種類（NEM か ANM）を選ぶので、負の νΣf は推測せずにエラー（ADR-0102）
    let mut m = material();
    m.nu_fission[0] = -1.0e-3;
    let p = slab(vec![10.0; 2], vec![m], vec![0, 0], Boundary::Reflective);
    let err = solve(&p, &Options::default());
    assert!(matches!(err, Err(NemError::InvalidProblem(_))), "{err:?}");
}

#[test]
fn get_adfs_reproduce_inner_faces_with_transverse_leakage() {
    // 2×2 の粗いノード（細かいノード 2.5 cm × 4 × 4 ずつ）。+x と +y が真空なので、内側の面にも横方向の
    // 漏れがあり、1ノードの問題の内側の辺は隣の粗いノードの幅と漏れで2次式を作る（ADR-0095）
    let mut boundaries = [Boundary::Reflective; 6];
    boundaries[1] = Boundary::Marshak;
    boundaries[3] = Boundary::Marshak;
    let pattern = [0, 1, 2, 0, 1, 1, 0, 2];
    let fine = Problem {
        mesh: Mesh {
            dx: vec![2.5; 8],
            dy: vec![2.5; 8],
            dz: vec![20.0],
        },
        materials: three_fuels(),
        reflectors: Vec::new(),
        node_material: placed((0..64).map(|n| pattern[(n * 5 + n / 8) % 8])),
        boundaries,
        outline: None,
    };
    let (coarse, solution) = coarsen_with_get(&fine, [4, 4]);
    assert_get_equivalence(&coarse, &solution, [4, 4], 8);
}

#[test]
fn normalized_adfs_take_the_neighbor_adf() {
    // 反射体の定数のように、辺の ADF を f* = f_自分 / f_隣 で持ち、印を付けた辺では、core はその面の
    // ADF を f* × 隣のノードの面の ADF にする（ADR-0095）。隣の ADF が作ったときと同じなら、元の ADF の
    // 解と同じ
    let mut boundaries = [Boundary::Reflective; 6];
    boundaries[1] = Boundary::Marshak;
    let fine = Problem {
        mesh: Mesh {
            dx: vec![2.5; 16],
            dy: vec![20.0],
            dz: vec![20.0],
        },
        materials: three_fuels(),
        reflectors: Vec::new(),
        node_material: placed(vec![0, 1, 0, 1, 1, 0, 1, 0, 2, 2, 1, 1, 2, 2, 0, 0]),
        boundaries,
        outline: None,
    };
    let (coarse, _) = coarsen_with_get(&fine, [8, 1]);
    let raw = solve(&coarse, &Options::default()).unwrap();
    let mut normalized = coarse.clone();
    let fuel = normalized.materials[0].adf[1].clone();
    let reflector = &mut normalized.materials[1];
    reflector.adf[0] = reflector.adf[0]
        .iter()
        .zip(&fuel)
        .map(|(f, g)| f / g)
        .collect();
    reflector.normalized_adf[0] = true;
    let s = solve(&normalized, &Options::default()).unwrap();
    assert!((s.k - raw.k).abs() < 1e-12 * raw.k, "{} vs {}", s.k, raw.k);
    for (a, b) in s.flux.iter().flatten().zip(raw.flux.iter().flatten()) {
        assert!((a - b).abs() < 1e-10 * b, "{a} vs {b}");
    }
    // 隣の ADF を変えると f* はそれに掛かるが、面の条件は両側の ADF の比だけで決まるので、解は変わらない
    // （f* の方法が比を反射体の問題から取る理由）
    let mut other = normalized.clone();
    other.materials[0].adf[1] = fuel.iter().map(|f| 1.1 * f).collect();
    let k = solve(&other, &Options::default()).unwrap().k;
    assert!((k - raw.k).abs() < 1e-12 * raw.k, "{k} vs {}", raw.k);
    // z の面の f* も同じ（ADR-0128）
    let s = solve(&to_z(&normalized), &Options::default()).unwrap();
    assert!((s.k - raw.k).abs() < 1e-12 * raw.k, "{} vs {}", s.k, raw.k);
    for (a, b) in s.flux.iter().flatten().zip(raw.flux.iter().flatten()) {
        assert!((a - b).abs() < 1e-10 * b, "{a} vs {b}");
    }
}

#[test]
fn normalized_adfs_need_a_fissile_neighbor() {
    // f* は燃料の ADF で割った値なので、隣が核分裂のない物質ならエラー（上下を取り違えた軸方向の反射体の
    // 定数を見つける。ADR-0128）。x でも z でも同じ
    let mut star = reflector();
    star.normalized_adf[0] = true;
    let p = slab(
        vec![10.0; 3],
        vec![material(), star, reflector()],
        vec![0, 1, 1],
        Boundary::ZeroFlux,
    );
    for p in [p.clone(), to_z(&p)] {
        let err = solve(&p, &Options::default()).unwrap_err().to_string();
        assert!(err.contains("fission"), "{err}");
    }
    // 燃料に向いていれば解ける
    let mut q = p.clone();
    q.node_material = placed(vec![0, 1, 2]);
    for q in [q.clone(), to_z(&q)] {
        assert!(solve(&q, &Options::default()).is_ok());
    }
}

#[test]
fn normalized_adfs_toward_assembly_layers_take_unit_adfs() {
    // 反射体の位置の f* の辺が、x と y の面で集合体の位置の核分裂のない層（端栓、プレナム）に向くときは、両側の
    // ADF を 1 にする（ADR-0132。I 型の定数は燃料と組で作ったもので、核分裂のない層と組で作った定数はない）。
    // 下の行は燃料と反射体、上の行は核分裂のない集合体の層と反射体
    let mut star = reflector();
    star.normalized_adf[0] = true;
    star.adf[0] = vec![1.3, 0.6];
    // 集合体の層は軸方向の定数（スラブの z の厚さを持つ。ADR-0128）
    let mut layer = reflector();
    layer.adf[1] = vec![0.8, 1.4];
    layer.height = Some(20.0);
    let problem = |star: FewGroup, layer: FewGroup, reflectors: Vec<usize>| {
        let mut boundaries = [Boundary::Reflective; 6];
        boundaries[1] = Boundary::Marshak;
        plane(
            vec![20.0, 21.5],
            vec![20.0, 20.0],
            vec![material(), star, layer],
            vec![Some(0), Some(1), Some(2), Some(1)],
            boundaries,
            None,
        )
        .with_reflectors(reflectors)
    };
    let p = problem(star.clone(), layer.clone(), vec![1]);
    // 上の行の面だけ両側を 1 にした定数（下の行の反射体は f* のまま）
    let mut unit_layer = layer.clone();
    unit_layer.adf[1] = vec![1.0, 1.0];
    let mut unit_star = star.clone();
    unit_star.normalized_adf[0] = false;
    unit_star.adf[0] = vec![1.0, 1.0];
    let mut q = problem(star.clone(), unit_layer, vec![1]);
    q.materials.push(unit_star);
    q.node_material[3] = Some(3);
    let (s, t) = (
        solve(&p, &Options::default()).unwrap(),
        solve(&q, &Options::default()).unwrap(),
    );
    assert!((s.k - t.k).abs() < 1e-12 * t.k, "{} vs {}", s.k, t.k);
    for (a, b) in s.flux.iter().flatten().zip(t.flux.iter().flatten()) {
        assert!((a - b).abs() < 1e-10 * b, "{a} vs {b}");
    }
    // 核分裂のない層が反射体の位置なら、回転を取り違えた定数としてエラー
    let err = solve(
        &problem(star.clone(), layer.clone(), vec![1, 2]),
        &Options::default(),
    )
    .unwrap_err()
    .to_string();
    assert!(err.contains("rotation"), "{err}");
    // f* の側が反射体の位置でなければ（燃料の位置に置いた反射体の定数）、エラー
    let err = solve(
        &problem(star.clone(), layer.clone(), vec![]),
        &Options::default(),
    )
    .unwrap_err()
    .to_string();
    assert!(err.contains("fission"), "{err}");
    // 集合体の位置の核分裂のない定数でも、軸方向の層（z の厚さを持つ）でなければ（集合体の位置に置いた径方向の
    // 反射体の定数）、エラー
    let mut radial = layer;
    radial.height = None;
    let err = solve(&problem(star, radial, vec![1]), &Options::default())
        .unwrap_err()
        .to_string();
    assert!(err.contains("fission"), "{err}");
}

#[test]
fn heights_must_match_the_node() {
    // z の厚さを持つ（軸方向の）定数は、z の幅が丸めの範囲で同じノードにだけ置ける（ADR-0128）
    let mut p = cuboid([30.0, 30.0, 30.0], [1, 1, 3], [Boundary::Reflective; 6]);
    p.materials[0].height = Some(10.0);
    assert!(solve(&p, &Options::default()).is_ok());
    p.materials[0].height = Some(10.0 * (1.0 + 1e-14));
    assert!(solve(&p, &Options::default()).is_ok());
    for h in [9.9, 20.0] {
        p.materials[0].height = Some(h);
        let err = solve(&p, &Options::default()).unwrap_err().to_string();
        assert!(err.contains("height"), "{err}");
    }
}

/// x 方向の問題を y 方向に並べ直す（x と y の幅、境界条件、少数群定数の辺（左、右）と（下、上）を入れ替える）
fn to_y(p: &Problem) -> Problem {
    let mut q = p.clone();
    q.mesh.dx = p.mesh.dy.clone();
    q.mesh.dy = p.mesh.dx.clone();
    q.boundaries.swap(0, 2);
    q.boundaries.swap(1, 3);
    // 辺は（左、右、上、下、−z、+z）: 新しい下は元の左、新しい上は元の右、新しい左右は元の下と上。z はそのまま
    fn turn<T: Clone>(sides: &[T]) -> Vec<T> {
        [3, 2, 1, 0, 4, 5]
            .iter()
            .map(|&s| sides[s].clone())
            .collect()
    }
    for m in &mut q.materials {
        m.adf = turn(&m.adf);
        m.vacuum_ratio = turn(&m.vacuum_ratio);
        m.normalized_adf = turn(&m.normalized_adf);
    }
    q
}

#[test]
fn faces_between_reflector_sites_keep_their_adfs() {
    // 炉心の形の反射体の位置に置いた物質（`reflectors`）どうしの x と y の面も、燃料どうしの面と同じく各側の
    // 少数群定数の ADF を使う（ADR-0132）
    let mut a = reflector();
    a.adf = vec![
        vec![1.3, 0.8],
        vec![0.7, 1.6],
        vec![1.2, 0.9],
        vec![0.85, 1.4],
        vec![1.0, 1.0],
        vec![1.0, 1.0],
    ];
    let mut b = reflector();
    b.adf = vec![
        vec![1.1, 0.6],
        vec![0.9, 1.2],
        vec![1.05, 0.75],
        vec![0.95, 1.3],
        vec![1.0, 1.0],
        vec![1.0, 1.0],
    ];
    let x = |reflectors: Vec<usize>, a: &FewGroup, b: &FewGroup| {
        let mut p = slab(
            vec![20.0, 21.5, 21.5],
            vec![material(), a.clone(), b.clone()],
            vec![0, 1, 2],
            Boundary::Reflective,
        );
        p.boundaries[1] = Boundary::Marshak;
        p.reflectors = reflectors;
        p
    };
    let unit = |m: &FewGroup, side: usize| {
        let mut m = m.clone();
        m.adf[side] = vec![1.0; 2];
        m
    };
    let same = |p: &Problem, q: &Problem| {
        let (s, t) = (
            solve(p, &Options::default()).unwrap(),
            solve(q, &Options::default()).unwrap(),
        );
        let close = (s.k - t.k).abs() < 1e-12 * t.k
            && s.flux
                .iter()
                .flatten()
                .zip(t.flux.iter().flatten())
                .all(|(u, v)| (u - v).abs() < 1e-10 * v);
        (close, s.k, t.k)
    };
    for turn in [false, true] {
        let along = |p: Problem| if turn { to_y(&p) } else { p };
        // a の右（+x）と b の左（−x）の面
        let marked = along(x(vec![1, 2], &a, &b));
        let plain = along(x(vec![], &a, &b));
        let (close, k, k_plain) = same(&marked, &plain);
        assert!(close, "turn {turn}: {k} vs {k_plain}");
        // ADF は効いている（1 にすれば解が変わる）
        let unit = along(x(vec![1, 2], &unit(&a, 1), &unit(&b, 0)));
        assert!(!same(&marked, &unit).0, "turn {turn}");
    }
    // 範囲の外の番号はエラー
    let err = solve(&x(vec![3], &a, &b), &Options::default());
    assert!(matches!(err, Err(NemError::InvalidProblem(_))), "{err:?}");
    // f* の辺は燃料の ADF で割った値なので、反射体の位置に向いていれば（回転を取り違えた定数）エラー
    let mut turned = a.clone();
    turned.normalized_adf[1] = true;
    let err = solve(&x(vec![1, 2], &turned, &b), &Options::default());
    assert!(matches!(err, Err(NemError::InvalidProblem(_))), "{err:?}");
}

#[test]
fn normalized_adfs_need_a_plain_neighbor() {
    // f* の辺は隣のノードの ADF を掛けるので、外周の面と、隣も f* の面はエラー
    let mut p = two_material_slab(vec![vec![1.0, 1.0]; 6], vec![vec![1.0, 1.0]; 6]);
    p.materials[0].normalized_adf[0] = true;
    assert!(matches!(
        solve(&p, &Options::default()),
        Err(NemError::InvalidProblem(_))
    ));
    let mut p = two_material_slab(vec![vec![1.0, 1.0]; 6], vec![vec![1.0, 1.0]; 6]);
    p.materials[0].normalized_adf[1] = true;
    p.materials[1].normalized_adf[0] = true;
    let err = solve(&p, &Options::default());
    assert!(matches!(err, Err(NemError::InvalidProblem(_))), "{err:?}");
    // 片側だけなら解ける
    p.materials[0].normalized_adf[1] = false;
    assert!(solve(&p, &Options::default()).is_ok());
}

/// z が1層（20 cm、上下は反射）の nx × ny の格子の問題。`sites` は y の下の行から、各行は x の小さい方から
fn plane(
    dx: Vec<f64>,
    dy: Vec<f64>,
    materials: Vec<FewGroup>,
    sites: Vec<Option<usize>>,
    boundaries: [Boundary; 6],
    outline: Option<Boundary>,
) -> Problem {
    Problem {
        mesh: Mesh {
            dx,
            dy,
            dz: vec![20.0],
        },
        materials,
        reflectors: Vec::new(),
        node_material: sites,
        boundaries,
        outline,
    }
}

/// どの辺にも真空の比と、辺ごとに違う ADF を持つ反射体
fn sided_reflector() -> FewGroup {
    let mut r = reflector();
    // 辺（左、右、上、下、−z、+z）
    r.adf = vec![
        vec![1.1, 0.8],
        vec![1.3, 0.7],
        vec![0.9, 1.2],
        vec![1.2, 0.6],
        vec![1.0, 1.0],
        vec![1.0, 1.0],
    ];
    r.vacuum_ratio = vec![
        Some(vec![0.45, 0.40]),
        Some(vec![0.50, 0.35]),
        Some(vec![0.42, 0.38]),
        Some(vec![0.47, 0.33]),
        None,
        None,
    ];
    r
}

#[test]
fn outline_faces_match_the_cuboid_faces() {
    // ノードのない位置に面する面は、同じノードを直方体の外周に置いて同じ条件にした面と同じ解になる（ADR-0132）。
    // ノードのない位置の側の直方体の外周の条件は効かない
    let materials = vec![material(), sided_reflector()];
    let h = 21.5;
    for condition in [
        Boundary::Vacuum,
        Boundary::Marshak,
        Boundary::ZeroFlux,
        Boundary::Reflective,
    ] {
        // （直方体だけの問題、ノードのない位置を足した問題、比べる位置の対応）
        let mut x_plus = [Boundary::Reflective; 6];
        x_plus[1] = condition;
        let mut y_plus = [Boundary::Reflective; 6];
        y_plus[3] = condition;
        let mut x_minus = [Boundary::Reflective; 6];
        x_minus[0] = condition;
        // ノードのない位置の外の、直方体の外周の面だけを変える
        let far = |face: usize| {
            let mut b = [Boundary::Reflective; 6];
            b[face] = Boundary::ZeroFlux;
            b
        };
        let cases = [
            (
                plane(
                    vec![h; 2],
                    vec![h],
                    materials.clone(),
                    placed([0, 1]),
                    x_plus,
                    None,
                ),
                plane(
                    vec![h; 3],
                    vec![h],
                    materials.clone(),
                    vec![Some(0), Some(1), None],
                    far(1),
                    Some(condition),
                ),
                vec![0, 1],
            ),
            (
                plane(
                    vec![h],
                    vec![h; 2],
                    materials.clone(),
                    placed([0, 1]),
                    y_plus,
                    None,
                ),
                plane(
                    vec![h],
                    vec![h; 3],
                    materials.clone(),
                    vec![Some(0), Some(1), None],
                    far(3),
                    Some(condition),
                ),
                vec![0, 1],
            ),
            (
                plane(
                    vec![h; 2],
                    vec![h],
                    materials.clone(),
                    placed([1, 0]),
                    x_minus,
                    None,
                ),
                plane(
                    vec![h; 3],
                    vec![h],
                    materials.clone(),
                    vec![None, Some(1), Some(0)],
                    far(0),
                    Some(condition),
                ),
                vec![1, 2],
            ),
        ];
        for (cuboid, stepped, sites) in cases {
            let a = solve(&cuboid, &Options::default()).unwrap();
            let b = solve(&stepped, &Options::default()).unwrap();
            assert!(
                (a.k - b.k).abs() < 1e-12 * a.k,
                "{condition:?}: {} vs {}",
                a.k,
                b.k
            );
            for (n, &site) in sites.iter().enumerate() {
                assert_eq!(b.flux[site].len(), 2, "site {site}");
                for (x, y) in a.flux[n].iter().zip(&b.flux[site]) {
                    assert!((x - y).abs() < 1e-10 * x, "{condition:?}: {x} vs {y}");
                }
            }
            let empty = (0..b.flux.len()).find(|s| !sites.contains(s)).unwrap();
            assert!(b.flux[empty].is_empty());
            assert!(b.currents.iter().all(|c| c[empty].is_empty()));
            for (x, y) in a.leakage.iter().zip(&b.leakage) {
                assert!((x - y).abs() < 1e-10 * x.abs().max(1e-12), "{x} vs {y}");
            }
        }
    }
}

/// 4 × 4 の集合体の位置の四隅にノードを置かない、燃料だけの階段の炉心。各位置を `split` × `split` に分ける
fn stepped_core(split: usize, outline: Boundary) -> Problem {
    let n = 4 * split;
    let corner = |i: usize| i / split == 0 || i / split == 3;
    let sites = (0..n * n)
        .map(|s| (!(corner(s % n) && corner(s / n))).then_some(0))
        .collect();
    let mut boundaries = [outline; 6];
    boundaries[4] = Boundary::Reflective;
    boundaries[5] = Boundary::Reflective;
    plane(
        vec![20.0 / split as f64; n],
        vec![20.0 / split as f64; n],
        vec![material()],
        sites,
        boundaries,
        Some(outline),
    )
}

#[test]
fn stepped_core_keeps_its_symmetry_and_balance() {
    // 4 回の回転の対称性を持つ階段の炉心の中性子束は同じ対称性を持ち、生成 = 吸収 + 外周と外形からの漏れ
    let p = stepped_core(1, Boundary::Marshak);
    let s = solve(&p, &Options::default()).unwrap();
    for (site, phi) in s.flux.iter().enumerate() {
        let (i, j) = (site % 4, site / 4);
        let turned = &s.flux[(3 - j) + 4 * i];
        assert_eq!(phi.len(), turned.len());
        for (a, b) in phi.iter().zip(turned) {
            assert!((a - b).abs() < 1e-6 * a, "site {site}: {a} vs {b}");
        }
    }
    let m = material();
    let volume = 20.0 * 20.0 * 20.0;
    let (mut absorption, mut production) = (0.0, 0.0);
    for phi in &s.flux {
        for ((f, a), nf) in phi.iter().zip(&m.absorption).zip(&m.nu_fission) {
            absorption += volume * a * f;
            production += volume * nf * f;
        }
    }
    let leakage: f64 = s.leakage.iter().sum();
    assert!(leakage > 0.0);
    assert!(
        (production / s.k - absorption - leakage).abs() < 1e-6 * production,
        "{production} / {} vs {absorption} + {leakage}",
        s.k
    );
}

#[test]
fn stepped_core_converges_with_refinement() {
    // 外形の中性子束ゼロの面は、ノードを細かくすると k が上から収束する（外形の凹んだ角を含む）。凹んだ角の特異性の
    // ため直方体より遅く、20 cm のノードの k は 2.5 cm のノードより 2% 以上大きい
    let k = |split: usize| {
        solve(
            &stepped_core(split, Boundary::ZeroFlux),
            &Options::default(),
        )
        .unwrap()
        .k
    };
    let ks = [k(1), k(2), k(4), k(8)];
    for w in ks.windows(3) {
        assert!((w[2] - w[1]).abs() < 0.7 * (w[1] - w[0]).abs(), "{ks:?}");
    }
    assert!(ks.windows(2).all(|w| w[1] < w[0]), "{ks:?}");
    assert!(ks[0] > 1.02 * ks[3], "{ks:?}");
    // 粗いノードの外れは、外形の面の扱いでなく凹んだ角の NEM の近似から来る: 四隅を強い吸収体（中性子束が
    // ほぼ 0 になる）のノードで埋めた直方体も、同じ 20 cm のノードで同じ程度に大きい（ADR-0132）
    let mut black = material();
    black.diffusion = vec![0.02, 0.02];
    black.total = vec![50.0, 50.0];
    black.absorption = vec![50.0, 50.0];
    black.scatter = vec![vec![0.0; 2]; 2];
    black.nu_fission = vec![0.0; 2];
    black.fission = vec![0.0; 2];
    black.chi = vec![0.0; 2];
    let mut filled = stepped_core(1, Boundary::ZeroFlux);
    filled.materials.push(black);
    for m in &mut filled.node_material {
        m.get_or_insert(1);
    }
    filled.outline = None;
    let k_filled = solve(&filled, &Options::default()).unwrap().k;
    assert!(
        (k_filled - ks[0]).abs() < 0.2 * (ks[0] - ks[3]),
        "{k_filled} vs {ks:?}"
    );
}

#[test]
fn outline_needs_a_condition_and_ratios() {
    let h = 21.5;
    // ノードのない位置があるのに外形の条件がない
    let p = plane(
        vec![h; 3],
        vec![h],
        vec![material(), sided_reflector()],
        vec![Some(0), Some(1), None],
        [Boundary::Reflective; 6],
        None,
    );
    let err = solve(&p, &Options::default()).unwrap_err().to_string();
    assert!(err.contains("outline"), "{err}");
    // ノードのない位置がないのに外形の条件がある
    let p = plane(
        vec![h; 2],
        vec![h],
        vec![material(), sided_reflector()],
        placed([0, 1]),
        [Boundary::Reflective; 6],
        Some(Boundary::Vacuum),
    );
    let err = solve(&p, &Options::default()).unwrap_err().to_string();
    assert!(err.contains("outline"), "{err}");
    // 外形の真空の面のノードが、その辺の真空の比を持たない
    let mut r = sided_reflector();
    r.vacuum_ratio[1] = None;
    let p = plane(
        vec![h; 3],
        vec![h],
        vec![material(), r],
        vec![Some(0), Some(1), None],
        [Boundary::Reflective; 6],
        Some(Boundary::Vacuum),
    );
    let err = solve(&p, &Options::default()).unwrap_err().to_string();
    assert!(err.contains("vacuum ratio"), "{err}");
    // ノードのない位置は径方向の図の位置なので、どの高さでも同じ（z の層ごとに違えばエラー）
    let mut p = plane(
        vec![h; 2],
        vec![h],
        vec![material(), sided_reflector()],
        vec![Some(0), Some(1), Some(0), None],
        [Boundary::Reflective; 6],
        Some(Boundary::Vacuum),
    );
    p.mesh.dz = vec![10.0, 10.0];
    let err = solve(&p, &Options::default()).unwrap_err().to_string();
    assert!(err.contains("height"), "{err}");
    // ノードが1つもない
    let p = plane(
        vec![h],
        vec![h],
        vec![material()],
        vec![None],
        [Boundary::Reflective; 6],
        Some(Boundary::Vacuum),
    );
    assert!(solve(&p, &Options::default()).is_err());
}

impl Problem {
    fn with_reflectors(mut self, reflectors: Vec<usize>) -> Self {
        self.reflectors = reflectors;
        self
    }
}
