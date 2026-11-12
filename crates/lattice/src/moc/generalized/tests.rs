use super::*;
use crate::geometry::{Lattice, PinCell};
use crate::moc::tests::{
    checkerboard_materials, moderator, tight, two_group, two_spectra_fuel, vacuum_checkerboard,
    with_moments,
};
use crate::moc::{Solution, response_importance, solve, solve_adjoint_from, solve_with};
use crate::symmetry::{Symmetry, detect, require};
use crate::tracking::{Boundaries, Boundary};

/// 前向きの解と状態、随伴の状態
fn eigenpair(
    tracking: &Tracking,
    materials: &[&Material],
    options: &Options,
    symmetry: Option<&Detected>,
) -> (Solution, Start, Start) {
    let (solution, forward) = solve_with(tracking, materials, options, symmetry, None).unwrap();
    let (_, adjoint) = solve_adjoint_from(tracking, materials, options, symmetry, None).unwrap();
    (solution, forward, adjoint)
}

/// ⟨h, φ⟩ = Σ V h φ
fn inner(solution: &Solution, h: &[Vec<f64>]) -> f64 {
    solution
        .volumes
        .iter()
        .zip(h)
        .zip(&solution.flux)
        .map(|((v, h), f)| v * h.iter().zip(f).map(|(a, b)| a * b).sum::<f64>())
        .sum()
}

/// 比の応答 R = ⟨h₁, φ⟩ / ⟨h₂, φ⟩ の一般化随伴の源 S = h₁ / ⟨h₁, φ⟩ − h₂ / ⟨h₂, φ⟩（⟨S, φ⟩ = 0）
fn ratio_source(solution: &Solution, h1: &[Vec<f64>], h2: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let (a, b) = (inner(solution, h1), inner(solution, h2));
    h1.iter()
        .zip(h2)
        .map(|(x, y)| x.iter().zip(y).map(|(p, q)| p / a - q / b).collect())
        .collect()
}

/// 領域ごとの [群] の値を、`regions` の領域の群 `group` だけ `value(region)` にし、ほかを 0 にしたもの
fn indicator(
    count: usize,
    ng: usize,
    group: usize,
    regions: impl Fn(usize) -> Option<f64>,
) -> Vec<Vec<f64>> {
    (0..count)
        .map(|r| {
            let mut h = vec![0.0; ng];
            if let Some(v) = regions(r) {
                h[group] = v;
            }
            h
        })
        .collect()
}

/// ⟨Γ, F φ⟩ = Σ V Σ_項 (χ · Γ)(νΣf · φ) と、その項の絶対値の和
fn fission_overlap(solution: &Solution, materials: &[&Material], gamma: &[Vec<f64>]) -> (f64, f64) {
    let mut sum = 0.0;
    let mut size = 0.0;
    for (r, m) in materials.iter().enumerate() {
        for t in &m.fission_terms {
            let emitted: f64 = t.chi.iter().zip(&gamma[r]).map(|(c, x)| c * x).sum();
            let rate: f64 = t
                .nu_fission
                .iter()
                .zip(&solution.flux[r])
                .map(|(n, x)| n * x)
                .sum();
            let absolute: f64 = t.chi.iter().zip(&gamma[r]).map(|(c, x)| c * x.abs()).sum();
            sum += solution.volumes[r] * emitted * rate;
            size += solution.volumes[r] * absolute * rate;
        }
    }
    (sum, size)
}

fn settings(tolerance: f64) -> GeneralizedOptions {
    GeneralizedOptions {
        method: GeneralizedMethod::SourceIteration,
        tolerance,
        max_iterations: 100_000,
        batch: false,
    }
}

/// 2 群の燃料（核分裂のスペクトル χ = (0.8, 0.2)、上方散乱あり）
fn spread_fuel() -> Material {
    let mut m = two_group([0.5, 1.2], [[0.40, 0.06], [0.01, 0.9]], [0.008, 0.25]);
    m.fission_terms[0].chi = vec![0.8, 0.2];
    m
}

#[test]
fn generalized_adjoint_of_an_infinite_medium_solves_the_scattering_balance() {
    // 均質な無限媒質では Γ は平らで、(diag Σt − Σs) Γ − νΣf (χ · Γ) / k = S（Σs は [入射][出射]）。⟨Γ, Fφ⟩ = 0 は χ · Γ = 0 なので、
    // 核分裂の項が消え、Γ = (diag Σt − Σs)⁻¹ S
    let cell = PinCell::new(1.26, vec![0.3, 0.45]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.1).unwrap();
    let m = spread_fuel();
    let materials = [&m, &m, &m];
    let options = tight(0);
    let (solution, forward, adjoint) = eigenpair(&tracking, &materials, &options, None);
    // ⟨S, φ⟩ = 0 の平らな源: S ∝ (φ₂, −φ₁)
    let phi = &solution.flux[0];
    let source = vec![vec![phi[1], -phi[0]]; 3];
    let gamma = solve_generalized(
        &tracking,
        &materials,
        &options,
        None,
        &forward,
        &adjoint,
        std::slice::from_ref(&source),
        &settings(1e-11),
    )
    .unwrap();
    let a = [
        [m.total[0] - m.scatter[0][0], -m.scatter[0][1]],
        [-m.scatter[1][0], m.total[1] - m.scatter[1][1]],
    ];
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    let s = &source[0];
    let expected = [
        (a[1][1] * s[0] - a[0][1] * s[1]) / det,
        (a[0][0] * s[1] - a[1][0] * s[0]) / det,
    ];
    for flux in &gamma[0].flux {
        for (x, e) in flux.iter().zip(&expected) {
            assert!((x - e).abs() < 1e-8 * e.abs(), "{flux:?} vs {expected:?}");
        }
    }
}

#[test]
fn generalized_adjoint_is_orthogonal_to_the_forward_fission_source() {
    // 真空の外周、上方散乱、2 つの核分裂のスペクトル。応答は燃料の領域の熱群の全反応率と、格子全体の熱群の中性子束の比
    let tracking = vacuum_checkerboard(2);
    let (f, m) = (two_spectra_fuel(), moderator());
    let materials = checkerboard_materials(2, &f, &m);
    let options = tight(0);
    let (solution, forward, adjoint) = eigenpair(&tracking, &materials, &options, None);
    let n = materials.len();
    let h1 = indicator(n, 2, 1, |r| (r == 0).then_some(f.total[1]));
    let h2 = indicator(n, 2, 1, |_| Some(1.0));
    let source = ratio_source(&solution, &h1, &h2);
    let gamma = solve_generalized(
        &tracking,
        &materials,
        &options,
        None,
        &forward,
        &adjoint,
        &[source],
        &settings(1e-10),
    )
    .unwrap();
    let (overlap, size) = fission_overlap(&solution, &materials, &gamma[0].flux);
    assert!(size > 0.0);
    assert!(overlap.abs() < 1e-9 * size, "{overlap} of {size}");
}

/// 応答 R = ⟨h₁, φ⟩ / ⟨h₂, φ⟩ の、領域 `region` の物質に h δM を足した問題の ln R の h による微分の中心差分（h と h/2 の
/// Richardson の外挿）
fn central_log_response(
    tracking: &Tracking,
    materials: &[&Material],
    region: usize,
    d: &Material,
    options: &Options,
    h: [&[Vec<f64>]; 2],
    step: f64,
) -> f64 {
    let at = |s: f64| {
        let perturbed = materials[region].added(d, s).unwrap();
        let mut m = materials.to_vec();
        m[region] = &perturbed;
        let solution = solve(tracking, &m, options).unwrap();
        (inner(&solution, h[0]) / inner(&solution, h[1])).ln()
    };
    let difference = |s: f64| (at(s) - at(-s)) / (2.0 * s);
    let (coarse, fine) = (difference(step), difference(step / 2.0));
    fine + (fine - coarse) / 3.0
}

/// 物質 `m` の核分裂の項だけを持つ物質（全断面積と散乱は 0。νΣf の摂動の向き）
fn fission_only(m: &Material) -> Material {
    Material {
        fission: m.fission.clone(),
        nu_fission: m.nu_fission.clone(),
        fission_terms: m.fission_terms.clone(),
        anisotropic: vec![vec![vec![0.0; m.groups()]; m.groups()]; m.anisotropic.len()],
        ..Material::void(m.groups())
    }
}

/// 核分裂の摂動 δF による比の応答の 1 次の変化 δR / R = ⟨Γ, δF φ⟩ / k（全断面積と散乱、応答の重みは変えない）が、直接の
/// 摂動の中心差分と一致すること
fn assert_fission_derivative(scattering_order: usize) {
    let tracking = vacuum_checkerboard(2);
    let (f, m) = if scattering_order == 0 {
        (two_spectra_fuel(), moderator())
    } else {
        (
            with_moments(&two_spectra_fuel(), &[0.3]),
            with_moments(&moderator(), &[0.4]),
        )
    };
    let materials = checkerboard_materials(2, &f, &m);
    let numerator = (0, f.total[1]);
    assert_fission_derivative_in(&tracking, &materials, numerator, &tight(scattering_order));
}

/// `assert_fission_derivative` の問題を `tracking` と `materials`（領域 0 が燃料）で解く。応答は `numerator` の領域の熱群の中性子束に
/// 重みを掛けた反応率と、すべての領域（真空を含む）の熱群の中性子束の比
fn assert_fission_derivative_in(
    tracking: &Tracking,
    materials: &[&Material],
    numerator: (usize, f64),
    options: &Options,
) {
    let f = materials[0];
    let (solution, forward, adjoint) = eigenpair(tracking, materials, options, None);
    let n = materials.len();
    let h1 = indicator(n, 2, 1, |r| (r == numerator.0).then_some(numerator.1));
    let h2 = indicator(n, 2, 1, |_| Some(1.0));
    let source = ratio_source(&solution, &h1, &h2);
    let gamma = solve_generalized(
        tracking,
        materials,
        options,
        None,
        &forward,
        &adjoint,
        &[source],
        &settings(1e-11),
    )
    .unwrap();
    // 燃料の 1 領域（領域 0）の νΣf を変える
    let d = fission_only(f);
    let g = &gamma[0].flux[0];
    let first_order: f64 = solution.volumes[0]
        * d.fission_terms
            .iter()
            .map(|t| {
                let emitted: f64 = t.chi.iter().zip(g).map(|(c, x)| c * x).sum();
                let rate: f64 = t
                    .nu_fission
                    .iter()
                    .zip(&solution.flux[0])
                    .map(|(a, b)| a * b)
                    .sum();
                emitted * rate
            })
            .sum::<f64>()
        / solution.k;
    let direct = central_log_response(tracking, materials, 0, &d, options, [&h1, &h2], 0.01);
    assert!(
        (first_order - direct).abs() < 2e-7 * direct.abs(),
        "P{}: {first_order} vs {direct}",
        options.scattering_order
    );
}

#[test]
fn generalized_adjoint_gives_the_response_derivative_of_a_fission_change() {
    assert_fission_derivative(0);
}

#[test]
fn generalized_adjoint_gives_the_response_derivative_with_anisotropic_scattering() {
    assert_fission_derivative(1);
}

#[test]
fn generalized_adjoint_takes_sources_in_void_regions() {
    // 燃料と被覆管の間の真空の隙間（領域 1）。応答は隙間の熱群の中性子束と、真空の領域を含む格子全体の熱群の中性子束の比
    let pin = PinCell::new(1.26, vec![0.4096, 0.418, 0.475]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap();
    let vacuum = Boundaries {
        left: Boundary::Vacuum,
        right: Boundary::Vacuum,
        top: Boundary::Vacuum,
        bottom: Boundary::Vacuum,
    };
    let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, vacuum).unwrap();
    let (f, m, gap) = (two_spectra_fuel(), moderator(), Material::void(2));
    let pin_materials = [&f, &gap, &m, &m];
    let materials: Vec<&Material> = [&pin_materials[..], &[&m], &[&m], &pin_materials[..]].concat();
    assert_fission_derivative_in(&tracking, &materials, (1, 1.0), &tight(0));
}

#[test]
fn generalized_adjoint_with_symmetry_matches_the_full_lattice() {
    let tracking = vacuum_checkerboard(3);
    let (f, m) = (two_spectra_fuel(), moderator());
    let materials = checkerboard_materials(3, &f, &m);
    let options = tight(0);
    let region_names: Vec<String> = materials
        .iter()
        .map(|&x| if std::ptr::eq(x, &f) { "f" } else { "m" }.to_string())
        .collect();
    let verified = detect(&tracking, &region_names).unwrap();
    assert_ne!(verified.symmetry(), Symmetry::NONE);
    let n = materials.len();
    let fuel = |r: usize| std::ptr::eq(materials[r], &f).then_some(f.total[0]);
    let h1 = indicator(n, 2, 0, fuel);
    let h2 = indicator(n, 2, 0, |_| Some(1.0));
    let solve_as = |symmetry: Option<&Detected>| {
        let (solution, forward, adjoint) = eigenpair(&tracking, &materials, &options, symmetry);
        let source = ratio_source(&solution, &h1, &h2);
        solve_generalized(
            &tracking,
            &materials,
            &options,
            symmetry,
            &forward,
            &adjoint,
            &[source],
            &settings(1e-10),
        )
        .unwrap()
        .remove(0)
    };
    let (full, reduced) = (solve_as(None), solve_as(Some(&verified)));
    let scale = full
        .flux
        .iter()
        .flatten()
        .fold(0.0_f64, |a, x| a.max(x.abs()));
    for (a, b) in full
        .flux
        .iter()
        .flatten()
        .zip(reduced.flux.iter().flatten())
    {
        assert!((a - b).abs() < 1e-7 * scale, "{a} vs {b}");
    }
}

#[test]
fn generalized_adjoint_rejects_a_source_that_is_not_orthogonal() {
    let tracking = vacuum_checkerboard(2);
    let (f, m) = (two_spectra_fuel(), moderator());
    let materials = checkerboard_materials(2, &f, &m);
    let options = tight(0);
    let (_, forward, adjoint) = eigenpair(&tracking, &materials, &options, None);
    // 応答の分子の重みだけ（⟨S, φ⟩ > 0）
    let h1 = indicator(materials.len(), 2, 1, |r| (r == 0).then_some(1.0));
    let result = solve_generalized(
        &tracking,
        &materials,
        &options,
        None,
        &forward,
        &adjoint,
        &[h1],
        &settings(1e-10),
    );
    assert!(
        matches!(result, Err(MocError::NonOrthogonalSource { index: 0, .. })),
        "{result:?}"
    );
}

#[test]
fn generalized_adjoint_rejects_invalid_sources() {
    let tracking = vacuum_checkerboard(2);
    let (f, m) = (two_spectra_fuel(), moderator());
    let materials = checkerboard_materials(2, &f, &m);
    let options = tight(0);
    let (solution, forward, adjoint) = eigenpair(&tracking, &materials, &options, None);
    let n = materials.len();
    let h1 = indicator(n, 2, 1, |r| (r == 0).then_some(1.0));
    let h2 = indicator(n, 2, 1, |_| Some(1.0));
    let good = ratio_source(&solution, &h1, &h2);
    let run = |sources: &[Vec<Vec<f64>>]| {
        solve_generalized(
            &tracking,
            &materials,
            &options,
            None,
            &forward,
            &adjoint,
            sources,
            &settings(1e-10),
        )
    };
    let short = good[..n - 1].to_vec();
    assert!(matches!(
        run(&[good.clone(), short]),
        Err(MocError::SourceShape { index: 1 })
    ));
    let mut nan = good.clone();
    nan[2][0] = f64::NAN;
    assert!(matches!(
        run(&[nan]),
        Err(MocError::InvalidSource {
            index: 0,
            region: 2
        })
    ));
    // 前向きの状態の代わりに、別のトラッキングの状態
    let other = vacuum_checkerboard(3);
    let (f3, m3) = (two_spectra_fuel(), moderator());
    let materials3 = checkerboard_materials(3, &f3, &m3);
    let (_, foreign, _) = eigenpair(&other, &materials3, &options, None);
    assert!(matches!(
        solve_generalized(
            &tracking,
            &materials,
            &options,
            None,
            &foreign,
            &adjoint,
            &[good],
            &settings(1e-10),
        ),
        Err(MocError::ForeignStart)
    ));
}

#[test]
fn accelerated_methods_reach_the_source_iteration_solution_in_fewer_sweeps() {
    // 鏡面反射の 3×3 の市松の格子（燃料のピンと水のセル。熱群の散乱の比が 1 に近く、加速のない反復が遅い）。粗メッシュの補正は
    // セルの間の分布を直すので、セルが 1 つのピンセルでは効かない（ADR-0247）。P1 の散乱と 1/4 の対称性も通す
    let pin = PinCell::new(1.26, vec![0.25, 0.4096, 0.55]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let layout = (0..3)
        .map(|i| (0..3).map(|j| (i + j) % 2).collect())
        .collect();
    let lattice = Lattice::new(vec![pin, water], layout).unwrap();
    let tracking = Tracking::new(&lattice, 16, 0.05).unwrap();
    let (f, m) = (
        with_moments(&two_spectra_fuel(), &[0.3]),
        with_moments(&moderator(), &[0.4]),
    );
    let materials: Vec<&Material> = (0..9)
        .flat_map(|c| {
            if c % 2 == 0 {
                vec![&f, &f, &m, &m]
            } else {
                vec![&m]
            }
        })
        .collect();
    let options = tight(1);
    let verified = require(&tracking, Symmetry::QUARTER).unwrap();
    for symmetry in [None, Some(&verified)] {
        let (solution, forward, adjoint) = eigenpair(&tracking, &materials, &options, symmetry);
        let n = materials.len();
        let h1 = indicator(n, 2, 1, |r| {
            std::ptr::eq(materials[r], &f).then_some(f.total[1])
        });
        let h2 = indicator(n, 2, 1, |_| Some(1.0));
        let source = ratio_source(&solution, &h1, &h2);
        let solve_by = |method: GeneralizedMethod| {
            let settings = GeneralizedOptions {
                method,
                ..settings(1e-11)
            };
            solve_generalized(
                &tracking,
                &materials,
                &options,
                symmetry,
                &forward,
                &adjoint,
                std::slice::from_ref(&source),
                &settings,
            )
            .unwrap()
            .remove(0)
        };
        let reference = solve_by(GeneralizedMethod::SourceIteration);
        let scale = reference
            .flux
            .iter()
            .flatten()
            .fold(0.0_f64, |a, x| a.max(x.abs()));
        for method in [
            GeneralizedMethod::Cmfd(CmfdOptions::default()),
            GeneralizedMethod::Gmres { restart: 20 },
            GeneralizedMethod::CmfdGmres {
                restart: 20,
                cmfd: CmfdOptions::default(),
            },
        ] {
            let accelerated = solve_by(method);
            for (a, b) in accelerated
                .flux
                .iter()
                .flatten()
                .zip(reference.flux.iter().flatten())
            {
                assert!((a - b).abs() < 1e-8 * scale, "{method:?}: {a} vs {b}");
            }
            assert!(
                accelerated.iterations < reference.iterations,
                "{method:?}: {} sweeps, source iteration {}",
                accelerated.iterations,
                reference.iterations
            );
        }
    }
}

#[test]
fn sources_swept_together_reach_the_same_solutions() {
    // 一緒にスイープしても、源ごとに解いた解と同じ（どの方法でも）。P1 の散乱と 1/4 の対称性で、3 つの源
    let pin = PinCell::new(1.26, vec![0.25, 0.4096, 0.55]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let layout = (0..3)
        .map(|i| (0..3).map(|j| (i + j) % 2).collect())
        .collect();
    let lattice = Lattice::new(vec![pin, water], layout).unwrap();
    let tracking = Tracking::new(&lattice, 16, 0.05).unwrap();
    let (f, m) = (
        with_moments(&two_spectra_fuel(), &[0.3]),
        with_moments(&moderator(), &[0.4]),
    );
    let materials: Vec<&Material> = (0..9)
        .flat_map(|c| {
            if c % 2 == 0 {
                vec![&f, &f, &m, &m]
            } else {
                vec![&m]
            }
        })
        .collect();
    let options = tight(1);
    let verified = require(&tracking, Symmetry::QUARTER).unwrap();
    let (solution, forward, adjoint) = eigenpair(&tracking, &materials, &options, Some(&verified));
    let n = materials.len();
    let fuel = |r: usize| std::ptr::eq(materials[r], &f);
    let all = |_: usize| Some(1.0);
    let sources = [
        ratio_source(
            &solution,
            &indicator(n, 2, 1, |r| fuel(r).then_some(f.total[1])),
            &indicator(n, 2, 1, all),
        ),
        ratio_source(
            &solution,
            &indicator(n, 2, 0, |r| fuel(r).then_some(f.nu_fission[0])),
            &indicator(n, 2, 0, all),
        ),
        ratio_source(
            &solution,
            &indicator(n, 2, 0, |r| Some(materials[r].scatter[0][1])),
            &indicator(n, 2, 0, all),
        ),
    ];
    for method in [
        GeneralizedMethod::SourceIteration,
        GeneralizedMethod::Cmfd(CmfdOptions::default()),
        GeneralizedMethod::Gmres { restart: 20 },
        GeneralizedMethod::CmfdGmres {
            restart: 20,
            cmfd: CmfdOptions::default(),
        },
    ] {
        let solve_with_batch = |batch: bool| {
            let settings = GeneralizedOptions {
                method,
                batch,
                ..settings(1e-11)
            };
            solve_generalized(
                &tracking,
                &materials,
                &options,
                Some(&verified),
                &forward,
                &adjoint,
                &sources,
                &settings,
            )
            .unwrap()
        };
        let (each, together) = (solve_with_batch(false), solve_with_batch(true));
        for (a, b) in each.iter().zip(&together) {
            let scale = a.flux.iter().flatten().fold(0.0_f64, |s, x| s.max(x.abs()));
            for (x, y) in a.flux.iter().flatten().zip(b.flux.iter().flatten()) {
                assert!((x - y).abs() < 1e-9 * scale, "{method:?}: {x} vs {y}");
            }
            assert_eq!(a.iterations, b.iterations, "{method:?}");
        }
    }
}

#[test]
fn gmres_stops_when_the_operator_maps_a_basis_vector_to_zero() {
    // 作用素が基底の値を 0 に移すと、Givens の回転と後退代入が 0 で割り、解が NaN になる
    let mut krylov = Krylov::new(vec![1.0, 0.0], 1e-12, 10);
    assert!(matches!(
        krylov.receive(vec![0.0, 0.0], 2),
        Err(MocError::GeneralizedDiverged { iterations: 2 })
    ));
}

#[test]
fn generalized_adjoint_checks_the_restart_and_takes_no_sources() {
    let tracking = vacuum_checkerboard(2);
    let (f, m) = (two_spectra_fuel(), moderator());
    let materials = checkerboard_materials(2, &f, &m);
    let options = tight(0);
    let (solution, forward, adjoint) = eigenpair(&tracking, &materials, &options, None);
    let n = materials.len();
    let source = ratio_source(
        &solution,
        &indicator(n, 2, 1, |r| (r == 0).then_some(1.0)),
        &indicator(n, 2, 1, |_| Some(1.0)),
    );
    let run = |method: GeneralizedMethod, batch: bool, sources: &[Vec<Vec<f64>>]| {
        let settings = GeneralizedOptions {
            method,
            batch,
            ..settings(1e-10)
        };
        solve_generalized(
            &tracking, &materials, &options, None, &forward, &adjoint, sources, &settings,
        )
    };
    // 再始動しない GMRES は基底を上限まで増やす
    for method in [
        GeneralizedMethod::Gmres { restart: 0 },
        GeneralizedMethod::CmfdGmres {
            restart: 0,
            cmfd: CmfdOptions::default(),
        },
    ] {
        assert!(matches!(
            run(method, false, std::slice::from_ref(&source)),
            Err(MocError::GeneralizedRestart)
        ));
    }
    for batch in [false, true] {
        for method in [
            GeneralizedMethod::SourceIteration,
            GeneralizedMethod::Gmres { restart: 10 },
        ] {
            assert!(run(method, batch, &[]).unwrap().is_empty());
        }
    }
}

/// 領域 `changes`（領域、δM）の物質を M + h δM にした問題の ln R（R = ⟨h₁, φ⟩ / ⟨h₂, φ⟩、h は変えない）の h による微分の中心差分
/// （h と h/2 の Richardson の外挿）
fn central_log_response_of(
    tracking: &Tracking,
    materials: &[&Material],
    changes: &[(usize, &Material)],
    options: &Options,
    h: [&[Vec<f64>]; 2],
    step: f64,
) -> f64 {
    let at = |s: f64| {
        let perturbed: Vec<Material> = changes
            .iter()
            .map(|&(region, d)| materials[region].added(d, s).unwrap())
            .collect();
        let mut m = materials.to_vec();
        for (&(region, _), p) in changes.iter().zip(&perturbed) {
            m[region] = p;
        }
        let solution = solve(tracking, &m, options).unwrap();
        (inner(&solution, h[0]) / inner(&solution, h[1])).ln()
    };
    let difference = |s: f64| (at(s) - at(-s)) / (2.0 * s);
    let (coarse, fine) = (difference(step), difference(step / 2.0));
    fine + (fine - coarse) / 3.0
}

/// 応答 R = ⟨h₁, φ⟩ / ⟨h₂, φ⟩ の一般化随伴 Γ から求めた、物質の変化 `changes` による δR / R の中性子束を通じた項
/// ⟨Γ, (δF / k − δL) φ⟩（`ResponseImportance::flux_term`）
fn response_flux_term(
    tracking: &Tracking,
    materials: &[&Material],
    changes: &[(usize, &Material)],
    options: &Options,
    symmetry: Option<&Detected>,
    h: [&[Vec<f64>]; 2],
) -> f64 {
    let (solution, forward, adjoint) = eigenpair(tracking, materials, options, symmetry);
    let source = ratio_source(&solution, h[0], h[1]);
    let gamma = solve_generalized(
        tracking,
        materials,
        options,
        symmetry,
        &forward,
        &adjoint,
        std::slice::from_ref(&source),
        &settings(1e-11),
    )
    .unwrap();
    let importance = response_importance(
        tracking, materials, options, symmetry, &forward, &gamma[0], &source,
    )
    .unwrap();
    importance.flux_term(changes).unwrap()
}

/// 吸収だけの 2 群の物質（全断面積だけを持つ）
fn absorber(total: [f64; 2]) -> Material {
    two_group(total, [[0.0, 0.0], [0.0, 0.0]], [0.0, 0.0])
}

#[test]
fn response_derivative_of_a_total_cross_section_change_uses_the_angular_fluxes() {
    // 吸収体を水の領域に足す。応答は燃料の領域の熱群の全反応率と、格子全体の熱群の中性子束の比。全断面積の摂動は、Γ と φ の
    // 角度中性子束の積を区間ごとに積分した値で効く（k の感度と同じ。ADR-0247 の「摂動の内積」）
    let tracking = vacuum_checkerboard(2);
    let (f, m) = (two_spectra_fuel(), moderator());
    let materials = checkerboard_materials(2, &f, &m);
    let options = tight(0);
    let n = materials.len();
    let h1 = indicator(n, 2, 1, |r| (r == 0).then_some(f.total[1]));
    let h2 = indicator(n, 2, 1, |_| Some(1.0));
    let d = absorber([0.3, 0.8]);
    let changes = [(1, &d)];
    let adjoint = response_flux_term(&tracking, &materials, &changes, &options, None, [&h1, &h2]);
    let direct =
        central_log_response_of(&tracking, &materials, &changes, &options, [&h1, &h2], 0.05);
    assert!(
        (adjoint - direct).abs() < 1e-6 * direct.abs(),
        "{adjoint} vs {direct}"
    );
}

#[test]
fn response_derivative_matches_central_differences_of_all_reactions() {
    // 燃料の 1 領域の全反応（全断面積、散乱、核分裂）を変える。真空の外周、上方散乱、2 つの核分裂のスペクトル
    let tracking = vacuum_checkerboard(2);
    let (f, m) = (two_spectra_fuel(), moderator());
    let materials = checkerboard_materials(2, &f, &m);
    let options = tight(0);
    let n = materials.len();
    let h1 = indicator(n, 2, 0, |r| (r == 3).then_some(m.scatter[0][1]));
    let h2 = indicator(n, 2, 0, |_| Some(1.0));
    let d = f.scaled(&[0.2, 0.2]).unwrap();
    let changes = [(0, &d)];
    let adjoint = response_flux_term(&tracking, &materials, &changes, &options, None, [&h1, &h2]);
    let direct =
        central_log_response_of(&tracking, &materials, &changes, &options, [&h1, &h2], 0.05);
    assert!(
        (adjoint - direct).abs() < 1e-6 * direct.abs(),
        "{adjoint} vs {direct}"
    );
}

#[test]
fn response_derivative_matches_central_differences_with_anisotropic_scattering() {
    // P1 の散乱の摂動は、Γ のモーメントの (−1)ˡ の符号で効く
    let tracking = vacuum_checkerboard(2);
    let (f, m) = (
        with_moments(&two_spectra_fuel(), &[0.3]),
        with_moments(&moderator(), &[0.4]),
    );
    let materials = checkerboard_materials(2, &f, &m);
    let options = tight(1);
    let n = materials.len();
    let h1 = indicator(n, 2, 1, |r| (r == 0).then_some(f.total[1]));
    let h2 = indicator(n, 2, 1, |_| Some(1.0));
    let d = with_moments(&m.scaled(&[0.2, 0.2]).unwrap(), &[0.9]);
    let changes = [(1, &d)];
    let adjoint = response_flux_term(&tracking, &materials, &changes, &options, None, [&h1, &h2]);
    let direct =
        central_log_response_of(&tracking, &materials, &changes, &options, [&h1, &h2], 0.05);
    assert!(
        (adjoint - direct).abs() < 1e-6 * direct.abs(),
        "{adjoint} vs {direct}"
    );
}

#[test]
fn response_derivative_carries_the_generalized_adjoint_across_void() {
    // 燃料と被覆管の間の真空の隙間（領域 1）に源を持つ応答（隙間の熱群の中性子束と、真空を含む格子全体の熱群の中性子束の比）。
    // Γ の角度中性子束は真空の区間で q / 4π の傾きで増えるので、隙間の外側の被覆管の全断面積の項はその分を含む
    let pin = PinCell::new(1.26, vec![0.4096, 0.418, 0.475]).unwrap();
    let water = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap();
    let vacuum = Boundaries {
        left: Boundary::Vacuum,
        right: Boundary::Vacuum,
        top: Boundary::Vacuum,
        bottom: Boundary::Vacuum,
    };
    let tracking = Tracking::with_boundaries(&lattice, 16, 0.05, vacuum).unwrap();
    let (f, m, gap) = (two_spectra_fuel(), moderator(), Material::void(2));
    let pin_materials = [&f, &gap, &m, &m];
    let materials: Vec<&Material> = [&pin_materials[..], &[&m], &[&m], &pin_materials[..]].concat();
    let options = tight(0);
    let n = materials.len();
    let h1 = indicator(n, 2, 1, |r| (r == 1).then_some(1.0));
    let h2 = indicator(n, 2, 1, |_| Some(1.0));
    let d = absorber([0.3, 0.8]);
    let changes = [(2, &d)];
    let adjoint = response_flux_term(&tracking, &materials, &changes, &options, None, [&h1, &h2]);
    let direct =
        central_log_response_of(&tracking, &materials, &changes, &options, [&h1, &h2], 0.05);
    assert!(
        (adjoint - direct).abs() < 1e-6 * direct.abs(),
        "{adjoint} vs {direct}"
    );
}

#[test]
fn response_derivative_with_symmetry_matches_the_full_lattice() {
    let tracking = vacuum_checkerboard(3);
    let (f, m) = (two_spectra_fuel(), moderator());
    let materials = checkerboard_materials(3, &f, &m);
    let options = tight(0);
    let region_names: Vec<String> = materials
        .iter()
        .map(|&x| if std::ptr::eq(x, &f) { "f" } else { "m" }.to_string())
        .collect();
    let verified = detect(&tracking, &region_names).unwrap();
    assert_ne!(verified.symmetry(), Symmetry::NONE);
    let n = materials.len();
    let fuel = |r: usize| std::ptr::eq(materials[r], &f).then_some(f.total[0]);
    let h1 = indicator(n, 2, 0, fuel);
    let h2 = indicator(n, 2, 0, |_| Some(1.0));
    let d = absorber([0.3, 0.8]);
    let changes: Vec<(usize, &Material)> = materials
        .iter()
        .enumerate()
        .filter(|&(_, &x)| std::ptr::eq(x, &m))
        .map(|(region, _)| (region, &d))
        .collect();
    let full = response_flux_term(&tracking, &materials, &changes, &options, None, [&h1, &h2]);
    let reduced = response_flux_term(
        &tracking,
        &materials,
        &changes,
        &options,
        Some(&verified),
        [&h1, &h2],
    );
    assert!(
        (full - reduced).abs() < 1e-7 * full.abs(),
        "{full} vs {reduced}"
    );
}

#[test]
fn response_importance_rejects_a_generalized_adjoint_of_another_problem() {
    let tracking = vacuum_checkerboard(2);
    let (f, m) = (two_spectra_fuel(), moderator());
    let materials = checkerboard_materials(2, &f, &m);
    let options = tight(0);
    let (solution, forward, adjoint) = eigenpair(&tracking, &materials, &options, None);
    let n = materials.len();
    let h1 = indicator(n, 2, 1, |r| (r == 0).then_some(f.total[1]));
    let h2 = indicator(n, 2, 1, |_| Some(1.0));
    let source = ratio_source(&solution, &h1, &h2);
    let gamma = solve_generalized(
        &tracking,
        &materials,
        &options,
        None,
        &forward,
        &adjoint,
        std::slice::from_ref(&source),
        &settings(1e-9),
    )
    .unwrap();
    // 別のトラッキングの前向きの状態
    let other = vacuum_checkerboard(3);
    let other_materials = checkerboard_materials(3, &f, &m);
    let (_, other_forward, _) = eigenpair(&other, &other_materials, &options, None);
    let other_source = vec![vec![0.0; 2]; other_materials.len()];
    assert!(matches!(
        response_importance(
            &other,
            &other_materials,
            &options,
            None,
            &other_forward,
            &gamma[0],
            &other_source
        ),
        Err(MocError::ForeignStart)
    ));
    // 源の形が問題と違う
    assert!(matches!(
        response_importance(
            &tracking,
            &materials,
            &options,
            None,
            &forward,
            &gamma[0],
            &source[1..]
        ),
        Err(MocError::SourceShape { index: 0 })
    ));
}

#[test]
fn response_importance_rejects_a_symmetry_of_another_tracking() {
    let tracking = vacuum_checkerboard(2);
    let (f, m) = (two_spectra_fuel(), moderator());
    let materials = checkerboard_materials(2, &f, &m);
    let options = tight(0);
    let (solution, forward, adjoint) = eigenpair(&tracking, &materials, &options, None);
    let n = materials.len();
    let h1 = indicator(n, 2, 1, |r| (r == 0).then_some(f.total[1]));
    let h2 = indicator(n, 2, 1, |_| Some(1.0));
    let source = ratio_source(&solution, &h1, &h2);
    let gamma = solve_generalized(
        &tracking,
        &materials,
        &options,
        None,
        &forward,
        &adjoint,
        std::slice::from_ref(&source),
        &settings(1e-9),
    )
    .unwrap();
    let other = vacuum_checkerboard(3);
    let other_names: Vec<String> = checkerboard_materials(3, &f, &m)
        .iter()
        .map(|&x| if std::ptr::eq(x, &f) { "f" } else { "m" }.to_string())
        .collect();
    let foreign = detect(&other, &other_names).unwrap();
    assert_ne!(foreign.symmetry(), Symmetry::NONE);
    assert!(matches!(
        response_importance(
            &tracking,
            &materials,
            &options,
            Some(&foreign),
            &forward,
            &gamma[0],
            &source
        ),
        Err(MocError::ForeignSymmetry)
    ));
}

/// 鏡面反射の 3×3 の格子で、燃料のピン 8 本が中央の水のセルを囲む（2A の案内管のまわりの形）。ピッチ 2.5 cm で、粗メッシュの
/// セルが光学的に厚い（熱群の水で Σt h = 5）
fn thick_cells_around_water() -> Tracking {
    let pin = PinCell::new(2.5, vec![0.25, 0.4096, 0.55]).unwrap();
    let water = PinCell::new(2.5, vec![]).unwrap();
    let layout = vec![vec![0, 0, 0], vec![0, 1, 0], vec![0, 0, 0]];
    Tracking::new(&Lattice::new(vec![pin, water], layout).unwrap(), 16, 0.05).unwrap()
}

/// `thick_cells_around_water` の領域の物質（ピンは燃料 2 つと水 2 つの領域、中央は水）
fn around_water<'a>(fuel: &'a Material, water: &'a Material) -> Vec<&'a Material> {
    (0..9)
        .flat_map(|c| {
            if c == 4 {
                vec![water]
            } else {
                vec![fuel, fuel, water, water]
            }
        })
        .collect()
}

#[test]
fn cmfd_reaches_the_gmres_solution_on_thick_cells_around_water() {
    // スイープは、φ† の成分を除いた状態からも φ† の成分を作る（除く重み Fφ は、スイープの左の固有ベクトル（前向きの散乱と
    // 核分裂の源）でない）。その成分は解の不定な分で、粗メッシュの補正の残差に入れると、細かい問題にない源を粗メッシュに与える。
    // この格子の高速群から熱群への散乱の応答では、そうした補正で反復が 1,000 スイープの中で収束しない
    let tracking = thick_cells_around_water();
    let (f, m) = (two_spectra_fuel(), moderator());
    let materials = around_water(&f, &m);
    let options = tight(0);
    let verified = require(&tracking, Symmetry::OCTANT).unwrap();
    let (solution, forward, adjoint) = eigenpair(&tracking, &materials, &options, Some(&verified));
    let n = materials.len();
    let source = ratio_source(
        &solution,
        &indicator(n, 2, 0, |r| Some(materials[r].scatter[0][1])),
        &indicator(n, 2, 0, |_| Some(1.0)),
    );
    let solve_by = |method: GeneralizedMethod| {
        let settings = GeneralizedOptions {
            method,
            max_iterations: 1000,
            ..settings(1e-10)
        };
        solve_generalized(
            &tracking,
            &materials,
            &options,
            Some(&verified),
            &forward,
            &adjoint,
            std::slice::from_ref(&source),
            &settings,
        )
        .map(|mut g| g.remove(0))
    };
    let reference = solve_by(GeneralizedMethod::Gmres { restart: 100 }).unwrap();
    let scale = reference
        .flux
        .iter()
        .flatten()
        .fold(0.0_f64, |a, x| a.max(x.abs()));
    let unaccelerated = solve_by(GeneralizedMethod::SourceIteration).unwrap();
    let cmfd = solve_by(GeneralizedMethod::Cmfd(CmfdOptions::default())).unwrap();
    for (a, b) in cmfd
        .flux
        .iter()
        .flatten()
        .zip(reference.flux.iter().flatten())
    {
        assert!((a - b).abs() < 1e-8 * scale, "{a} vs {b}");
    }
    assert!(
        3 * cmfd.iterations < unaccelerated.iterations,
        "{} sweeps, source iteration {}",
        cmfd.iterations,
        unaccelerated.iterations
    );
}

#[test]
fn gmres_of_the_cmfd_iteration_needs_fewer_sweeps_than_either_alone() {
    // CMFD の補正を足した反復の作用素の GMRES。補正の残る遅いモード（粗メッシュのセルの中の分布）を GMRES が除く
    let tracking = thick_cells_around_water();
    let (f, m) = (two_spectra_fuel(), moderator());
    let materials = around_water(&f, &m);
    let options = tight(0);
    let verified = require(&tracking, Symmetry::OCTANT).unwrap();
    let (solution, forward, adjoint) = eigenpair(&tracking, &materials, &options, Some(&verified));
    let n = materials.len();
    let source = ratio_source(
        &solution,
        &indicator(n, 2, 1, |r| {
            std::ptr::eq(materials[r], &f).then_some(f.total[1])
        }),
        &indicator(n, 2, 1, |_| Some(1.0)),
    );
    let solve_by = |method: GeneralizedMethod, tolerance: f64| {
        let settings = GeneralizedOptions {
            method,
            ..settings(tolerance)
        };
        solve_generalized(
            &tracking,
            &materials,
            &options,
            Some(&verified),
            &forward,
            &adjoint,
            std::slice::from_ref(&source),
            &settings,
        )
        .unwrap()
        .remove(0)
    };
    let reference = solve_by(GeneralizedMethod::Gmres { restart: 100 }, 1e-12);
    let scale = reference
        .flux
        .iter()
        .flatten()
        .fold(0.0_f64, |a, x| a.max(x.abs()));
    let gmres = solve_by(GeneralizedMethod::Gmres { restart: 100 }, 1e-10);
    let cmfd = solve_by(GeneralizedMethod::Cmfd(CmfdOptions::default()), 1e-10);
    let both = solve_by(
        GeneralizedMethod::CmfdGmres {
            restart: 100,
            cmfd: CmfdOptions::default(),
        },
        1e-10,
    );
    for (a, b) in both
        .flux
        .iter()
        .flatten()
        .zip(reference.flux.iter().flatten())
    {
        assert!((a - b).abs() < 1e-8 * scale, "{a} vs {b}");
    }
    assert!(
        both.iterations < gmres.iterations && both.iterations < cmfd.iterations,
        "{} sweeps, GMRES {}, CMFD {}",
        both.iterations,
        gmres.iterations,
        cmfd.iterations
    );
}
