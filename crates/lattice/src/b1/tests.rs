use super::*;

fn close(a: f64, b: f64, tolerance: f64) -> bool {
    (a - b).abs() <= tolerance * b.abs()
}

fn two_group(scatter1: [[f64; 2]; 2], nu_fission: [f64; 2]) -> Medium {
    Medium {
        total: vec![0.55, 1.4],
        scatter: vec![vec![0.48, 0.03], vec![0.001, 1.25]],
        scatter1: scatter1.iter().map(|r| r.to_vec()).collect(),
        nu_fission: nu_fission.to_vec(),
        chi: vec![1.0, 0.0],
    }
}

#[test]
fn gamma_series_joins_the_closed_form() {
    for s in [SERIES_LIMIT, -SERIES_LIMIT] {
        let series = gamma(s * (1.0 - 1e-9));
        let x = s.abs().sqrt();
        let closed = if s > 0.0 {
            s * x.atan() / (3.0 * (x - x.atan()))
        } else {
            -s * x.atanh() / (3.0 * (x.atanh() - x))
        };
        assert!(close(series, closed, 1e-10), "{s}: {series} vs {closed}");
    }
    assert_eq!(gamma(0.0), 1.0);
}

#[test]
fn isotropic_scattering_gives_one_third_total_at_zero_buckling() {
    let medium = two_group([[0.0; 2]; 2], [0.01, 0.3]);
    let mode = medium.mode(0.0).unwrap();
    for (d, t) in mode.diffusion().iter().zip(&medium.total) {
        assert!(close(*d, 1.0 / (3.0 * t), 1e-13), "{d} vs {t}");
    }
}

#[test]
fn isotropic_scattering_follows_gamma_at_finite_buckling() {
    // Σs1 = 0 なら中性子流の式は群ごとに閉じ、D_g = 1/(3 γ(B²/Σt²) Σt)
    let medium = two_group([[0.0; 2]; 2], [0.01, 0.3]);
    for buckling in [0.05, -0.2] {
        let mode = medium.mode(buckling).unwrap();
        for (d, t) in mode.diffusion().iter().zip(&medium.total) {
            let expected = 1.0 / (3.0 * gamma(buckling / (t * t)) * t);
            assert!(close(*d, expected, 1e-13), "{buckling}: {d} vs {expected}");
        }
    }
}

#[test]
fn diagonal_p1_scattering_gives_the_transport_cross_section() {
    let medium = two_group([[0.2, 0.0], [0.0, 0.6]], [0.01, 0.3]);
    let mode = medium.mode(0.0).unwrap();
    let expected = [1.0 / (3.0 * (0.55 - 0.2)), 1.0 / (3.0 * (1.4 - 0.6))];
    for (d, e) in mode.diffusion().iter().zip(expected) {
        assert!(close(*d, e, 1e-13), "{d} vs {e}");
    }
}

/// 1群の等方散乱の平面波の輸送方程式の分散関係から、臨界の x = B/Σt に対する c = (Σs + νΣf)/Σt を作る。
/// 正の B²: c = x / arctan(x)。負の B²（x = iy）: c = y / artanh(y)。釣り合いから D = Σt (c − 1) / B²
fn one_group(buckling_over_total2: f64) -> (Medium, f64) {
    let total = 0.8;
    let scatter = 0.5;
    let s: f64 = buckling_over_total2;
    let x = s.abs().sqrt();
    let c = if s > 0.0 { x / x.atan() } else { x / x.atanh() };
    let medium = Medium {
        total: vec![total],
        scatter: vec![vec![scatter]],
        scatter1: vec![vec![0.0]],
        nu_fission: vec![c * total - scatter],
        chi: vec![1.0],
    };
    let buckling = s * total * total;
    (medium, total * (c - 1.0) / buckling)
}

#[test]
fn one_group_critical_mode_matches_the_transport_dispersion_relation() {
    // 1群の等方散乱では B1 は厳密で、臨界のバックリングと D が輸送方程式の解析解に一致する
    for s in [0.25, -0.25, 2.0e-4, -2.0e-4] {
        let (medium, expected) = one_group(s);
        let mode = medium.critical().unwrap();
        let buckling = s * 0.8 * 0.8;
        assert!(
            close(mode.buckling, buckling, 1e-8),
            "{} vs {buckling}",
            mode.buckling
        );
        assert!(close(mode.diffusion()[0], expected, 1e-8));
        let gamma_form = 1.0 / (3.0 * gamma(s) * 0.8);
        assert!(close(mode.diffusion()[0], gamma_form, 1e-8));
    }
}

#[test]
fn critical_buckling_sign_follows_k_inf() {
    let p1 = [[0.15, 0.01], [0.0, 0.5]];
    for (nu_fission, positive) in [([0.01, 0.4], true), ([0.005, 0.2], false)] {
        let medium = two_group(p1, nu_fission);
        let k_inf = medium.mode(0.0).unwrap().k;
        assert_eq!(k_inf > 1.0, positive, "k∞ {k_inf}");
        let mode = medium.critical().unwrap();
        assert!((mode.k - 1.0).abs() < 1e-10, "k {}", mode.k);
        assert_eq!(mode.buckling > 0.0, positive, "B² {}", mode.buckling);
    }
}

#[test]
fn no_fission_is_an_error() {
    let medium = two_group([[0.0; 2]; 2], [0.0, 0.0]);
    assert!(matches!(medium.critical(), Err(B1Error::NoFission)));
}

#[test]
fn nan_is_rejected() {
    let mut medium = two_group([[0.0; 2]; 2], [0.01, 0.3]);
    assert!(matches!(
        medium.mode(f64::NAN),
        Err(B1Error::Imaginary { .. })
    ));
    medium.scatter1[1][0] = f64::NAN;
    assert!(matches!(medium.critical(), Err(B1Error::NotFinite)));
}

#[test]
fn negative_buckling_without_a_bracket_is_an_error() {
    // 群 0 は Σt が最小だが源がなく中性子束 0、核分裂は群 1 だけ（k∞ = 0.5）。−min(Σt)² の近くでも
    // モードは物理的で k < 1 のままで、区間 (−min(Σt)², 0) に臨界のバックリングがない
    let medium = Medium {
        total: vec![0.5, 2.0],
        scatter: vec![vec![0.1, 0.0], vec![0.0, 0.0]],
        scatter1: vec![vec![0.0; 2]; 2],
        nu_fission: vec![0.0, 1.0],
        chi: vec![0.0, 1.0],
    };
    assert!(matches!(
        medium.critical(),
        Err(B1Error::NoNegativeBracket { k, .. }) if k < 1.0
    ));
}

#[test]
fn modes_beyond_the_pole_count_as_too_negative() {
    // 1群（k∞ < 1）で、−min(Σt)² と臨界の B² の間には極があり、その先のモードは物理的でない。
    // 探索はその点を下端として扱い、臨界に収束する
    let (medium, _) = one_group(-0.25);
    let limit = NEGATIVE_LIMIT * medium.negative_limit();
    let beyond = medium.solve(0.9 * limit).unwrap();
    assert!(!beyond.physical(), "k {}", beyond.k);
    let mode = medium.critical().unwrap();
    assert!((mode.k - 1.0).abs() < 1e-10);
    assert!(mode.buckling > 0.9 * limit);
}

#[test]
fn zero_buckling_current_solves_the_p1_equation() {
    // 核分裂のない媒質でも使える。3 Σt d − 3 Σs1ᵀ d = φ
    let medium = two_group([[0.2, 0.05], [0.01, 0.6]], [0.0, 0.0]);
    let flux = [1.0, 2.0];
    let d = medium.zero_buckling_current(&flux).unwrap();
    for g in 0..2 {
        let lhs = 3.0 * medium.total[g] * d[g]
            - 3.0 * (0..2).map(|f| medium.scatter1[f][g] * d[f]).sum::<f64>();
        assert!((lhs - flux[g]).abs() < 1e-14, "{lhs} vs {}", flux[g]);
    }
}
