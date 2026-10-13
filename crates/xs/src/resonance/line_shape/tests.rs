use super::*;
use std::f64::consts::PI;

#[test]
fn parallel_rows_come_back_in_index_order() {
    let rows: Result<Vec<usize>, ()> = in_parallel(1000, |i| Ok(i * i));
    assert_eq!(rows.unwrap(), (0..1000).map(|i| i * i).collect::<Vec<_>>());
}

#[test]
fn parallel_rows_return_the_error_of_the_lowest_index() {
    for _ in 0..20 {
        let rows: Result<Vec<usize>, usize> =
            in_parallel(1000, |i| if i % 97 == 50 { Err(i) } else { Ok(i) });
        assert_eq!(rows, Err(50));
    }
}

#[test]
fn parallel_rows_stop_taking_rows_after_a_failure() {
    let calls = AtomicUsize::new(0);
    let rows: Result<Vec<usize>, usize> = in_parallel(100_000, |i| {
        calls.fetch_add(1, Ordering::Relaxed);
        if i == 0 {
            return Err(i);
        }
        std::thread::sleep(std::time::Duration::from_millis(1));
        Ok(i)
    });
    assert_eq!(rows, Err(0));
    let calls = calls.into_inner();
    assert!(calls < 1000, "{calls} rows after the first failed");
}

fn assert_rel(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol * expected.abs(),
        "{actual} vs {expected}"
    );
}

#[test]
fn jk_integrals_need_positive_beta() {
    for beta in [0.0, -1.0, f64::NAN] {
        assert!(matches!(
            jk_integrals(&[1.0, beta], 1.0),
            Err(NumericError::Domain(_))
        ));
    }
}

#[test]
fn natural_line_shape_without_doppler() {
    // θ → ∞ で ψ = 1/(1+x²)。a² = (1+β)/β として
    // J = π/(2aβ)、K = (1/β)[π/(2a) + π(a²+1)/(4a³)]
    let theta = 1.0e5;
    let betas = [1.0e-4, 1.0e-2, 1.0, 1.0e2];
    for (&beta, (j, k)) in betas.iter().zip(jk_integrals(&betas, theta).unwrap()) {
        let a = ((1.0 + beta) / beta).sqrt();
        assert_rel(j, PI / (2.0 * a * beta), 1e-8);
        assert_rel(
            k,
            (PI / (2.0 * a) + PI * (a * a + 1.0) / (4.0 * a * a * a)) / beta,
            1e-8,
        );
    }
    assert_rel(psi(0.7, theta), 1.0 / (1.0 + 0.49), 1e-8);
}

#[test]
fn dilute_limit() {
    // β → ∞ で J → (1/β)∫ψ dx = π/(2β)、K → 2J
    for theta in [0.01, 0.3, 10.0] {
        let [(j, k)] = jk_integrals(&[1.0e8], theta).unwrap()[..] else {
            panic!("one β")
        };
        assert_rel(j, PI / 2.0e8, 1e-7);
        assert_rel(k, PI / 1.0e8, 1e-7);
    }
}

/// J を x について Simpson 則で直接積分した値（x ≤ X、X より先は ψ ≈ 1/x² で足す）。
fn j_by_simpson(beta: f64, theta: f64) -> f64 {
    let x_max = 1.0e4;
    let n = 2_000_000;
    let h = x_max / n as f64;
    let f = |x: f64| {
        let p = psi(x, theta);
        p / (p + beta)
    };
    let mut sum = f(0.0) + f(x_max);
    for i in 1..n {
        sum += if i % 2 == 1 { 4.0 } else { 2.0 } * f(h * i as f64);
    }
    sum * h / 3.0 + 1.0 / (beta * x_max)
}

#[test]
fn doppler_dominated_line_matches_direct_integration() {
    // ドップラー幅が自然幅の 100 倍
    let theta = 0.01;
    let betas = [1.0e-3, 0.1];
    for (&beta, (j, _)) in betas.iter().zip(jk_integrals(&betas, theta).unwrap()) {
        assert_rel(j, j_by_simpson(beta, theta), 1e-6);
    }
}

#[test]
fn table_matches_direct_integration() {
    // 表の範囲から決まった疑似乱数で選んだ点で、補間と直接の積分の相対差が 1e-6 以下
    let mut seed: u64 = 12345;
    let mut uniform = || {
        seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    for _ in 0..400 {
        let beta = 10f64.powf(LOG_BETA.0 + (LOG_BETA.1 - LOG_BETA.0) * uniform());
        let theta = 10f64.powf(LOG_THETA.0 + (LOG_THETA.1 - LOG_THETA.0) * uniform());
        let (j, k) = jk_interpolated(beta, theta).unwrap();
        let [(jd, kd)] = jk_integrals(&[beta], theta).unwrap()[..] else {
            panic!("one β")
        };
        assert_rel(j, jd, 1e-6);
        assert_rel(k, kd, 1e-6);
    }
}

#[test]
fn table_uses_dilute_limit_above_its_range() {
    let (j, k) = jk_interpolated(1.0e12, 0.5).unwrap();
    assert_rel(j, PI / 2.0e12, 1e-12);
    assert_rel(k, PI / 1.0e12, 1e-12);
}

#[test]
fn table_uses_natural_line_shape_above_its_theta_range() {
    // θ = Γ/Δ が表より大きいとき（低温の構造材の非分離共鳴。Fe-58 の 550 K で θ ≈ 1.0e4）は、ドップラー
    // 拡がりのない Lorentz の線形状の J, K。表の端で値が跳ぶと、幅のゆらぎの適応積分（許容誤差 1e-10）が
    // 収束しにくいので、表の端の内と外で直接の積分との差が 1e-10 より小さいこと
    let top = 10f64.powf(LOG_THETA.1);
    for theta in [1.0025e4, 0.999 * top, 1.001 * top, 1.0e8] {
        for beta in [1.0e-6, 1.0e-4, 1.0e-2, 1.0, 1.0e2, 1.0e10] {
            let (j, k) = jk_interpolated(beta, theta).unwrap();
            let direct = jk_integrals(&[beta], theta).unwrap()[0];
            assert_rel(j, direct.0, 1e-10);
            assert_rel(k, direct.1, 1e-10);
        }
    }
    // 表の β より小さくても、Lorentz の線形状の J, K は使える
    let (j, k) = jk_interpolated(1.0e-8, 1.0e8).unwrap();
    let direct = jk_integrals(&[1.0e-8], 1.0e8).unwrap()[0];
    assert_rel(j, direct.0, 1e-10);
    assert_rel(k, direct.1, 1e-10);
}

#[test]
fn table_rejects_points_outside_its_range() {
    assert!(jk_interpolated(1.0e-8, 1.0).is_err());
    assert!(jk_interpolated(1.0e-8, 1.0e-13).is_err());
}

#[test]
fn table_covers_resonances_narrower_than_a_thousandth_of_the_doppler_width() {
    // θ = Γ/Δ が 1e-3 より小さい非分離共鳴（Nd-150、Sm-154、Gd-160 の 565〜1200 K、Nd-148 の 900 K と 1200 K。#0119）でも、
    // 補間が直接の積分と表の精度で合う
    for theta in [9.0e-4, 1.0e-4, 1.0e-6, 2.0e-9] {
        for beta in [1.0e-6, 1.0e-2, 1.0, 1.2e5, 1.9e8] {
            let (j, k) = jk_interpolated(beta, theta).unwrap();
            let direct = jk_integrals(&[beta], theta).unwrap()[0];
            assert_rel(j, direct.0, 1e-6);
            assert_rel(k, direct.1, 1e-6);
        }
    }
}

/// θ → 0 の線形状（ガウス形 ψ = (√π θ/2) e^{-u²}、u = θx/2）の (J, K)。b = 2β/(√π θ) > 1 で、g = e^{-u²} について
/// g/(g + b) と g(g + 2b)/(g + b)² を 1/b の冪に展開し、∫₀^∞ gⁿ du = √π/(2√n) で項ごとに積分した級数
fn gaussian_limit(beta: f64, theta: f64) -> (f64, f64) {
    let b = 2.0 * beta / (PI.sqrt() * theta);
    assert!(b > 1.0, "the series needs b > 1, got {b}");
    let scale = 2.0 / theta * 0.5 * PI.sqrt();
    let (mut j, mut k) = (0.0, 0.0);
    for n in 1..200 {
        let nf = n as f64;
        let sign = if n % 2 == 1 { 1.0 } else { -1.0 };
        let power = b.powi(-n);
        j += sign * power / nf.sqrt();
        k += sign * (nf + 1.0) * power / nf.sqrt();
        if (nf + 2.0) * power < 1.0e-18 {
            break;
        }
    }
    (scale * j, scale * k)
}

#[test]
fn quadrature_approaches_the_gaussian_limit_for_narrow_lines() {
    // ローレンツの裾による差は θ の程度。b = 2β/(√π θ) が 1 に近い芯の形まで、積分と級数が θ の数倍の内で合う
    for theta in [1.0e-6, 1.0e-8, 1.0e-10] {
        for beta in [2.5 * theta, 10.0 * theta, 1.0e3 * theta, 1.0] {
            let direct = jk_integrals(&[beta], theta).unwrap()[0];
            let (j, k) = gaussian_limit(beta, theta);
            assert_rel(direct.0, j, 10.0 * theta + 1.0e-11);
            assert_rel(direct.1, k, 10.0 * theta + 1.0e-11);
        }
    }
}

#[test]
fn table_takes_the_gaussian_limit_below_its_theta_range() {
    // Au-197 の非分離共鳴は Γγ = 1e-9 eV の系列を持ち、θ ≈ 1e-10 になる（#0119）。それより狭い線にも下限を置かない
    for theta in [9.0e-11, 1.0e-13, 1.0e-20] {
        for beta in [1.0e-6, 1.0, 1.0e8] {
            let (j, k) = jk_interpolated(beta, theta).unwrap();
            let (jg, kg) = gaussian_limit(beta, theta);
            assert_rel(j, jg, 1.0e-12);
            assert_rel(k, kg, 1.0e-12);
        }
    }
}
