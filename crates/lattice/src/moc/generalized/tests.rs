use super::*;
use crate::geometry::{Lattice, PinCell};
use crate::moc::tests::{
    checkerboard_materials, moderator, tight, two_group, two_spectra_fuel, vacuum_checkerboard,
    with_moments,
};
use crate::moc::{Solution, solve, solve_adjoint_from, solve_with};
use crate::symmetry::{Symmetry, detect};

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
        &[source.clone()],
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
    let options = tight(scattering_order);
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
        &settings(1e-11),
    )
    .unwrap();
    // 燃料の 1 領域（領域 0）の νΣf を変える
    let d = fission_only(&f);
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
    let direct = central_log_response(&tracking, &materials, 0, &d, &options, [&h1, &h2], 0.05);
    assert!(
        (first_order - direct).abs() < 1e-6 * direct.abs(),
        "P{scattering_order}: {first_order} vs {direct}"
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
