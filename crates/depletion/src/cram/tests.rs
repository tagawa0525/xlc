//! 解析解との比較。閾値は ADR-0161 の「確認」の値: 誤差の絶対値を初期の数密度の最大で割ったものが、
//! 実の固有値のチェーンで 1e-13 以下（両次数）、複素の固有値を持つ環で次数 48 は 1e-13、次数 16 は 1e-10 以下。
//! 親の 1e-9 倍しかない娘は、相対誤差で 1e-12 以下。

use super::*;

const ORDERS: [CramOrder; 2] = [CramOrder::Order16, CramOrder::Order48];
const REAL_SPECTRUM_TOL: f64 = 1e-13;
const COMPLEX_SPECTRUM_TOL_16: f64 = 1e-10;
const COMPLEX_SPECTRUM_TOL_48: f64 = 1e-13;
const SMALL_DAUGHTER_RTOL: f64 = 1e-12;

fn matrix(n: usize, entries: &[(usize, usize, f64)]) -> DepletionMatrix {
    DepletionMatrix::from_triplets(n, entries).unwrap()
}

/// 誤差の絶対値の最大を、初期の数密度の最大で割ったもの
fn scaled_error(got: &[f64], exact: &[f64], n0: &[f64]) -> f64 {
    let scale = n0.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    got.iter()
        .zip(exact)
        .map(|(g, e)| (g - e).abs())
        // f64::max は NaN を捨てるので、NaN の項（参照か結果の誤り）を誤差 0 と数えないよう、比べて畳む
        .fold(0.0, |m: f64, d| if d.is_nan() || d > m { d } else { m })
        / scale
}

/// 親（λ1）→ 娘（λ2）→ 安定な孫の、初めに親だけが 1 のときの解。娘と孫は打ち消しのない形
/// （孫は λ2 t か λ1 t が 10 以上のときに打ち消しが小さい）
fn parent_daughter_stable(l1: f64, l2: f64, t: f64) -> [f64; 3] {
    let parent = (-l1 * t).exp();
    // 遅い方の指数を外に出し、exp_m1 の引数を負に保つ（速い方を外に出すと 0 × ∞ になる）
    let daughter = if l2 > l1 {
        l1 / (l2 - l1) * parent * -(-(l2 - l1) * t).exp_m1()
    } else {
        l1 / (l1 - l2) * (-l2 * t).exp() * -(-(l1 - l2) * t).exp_m1()
    };
    let stable = (l2 * -(-l1 * t).exp_m1() - l1 * -(-l2 * t).exp_m1()) / (l2 - l1);
    [parent, daughter, stable]
}

fn chain3(l1: f64, l2: f64) -> DepletionMatrix {
    matrix(3, &[(0, 0, -l1), (1, 0, l1), (1, 1, -l2), (2, 1, l2)])
}

#[test]
fn one_nuclide_decays_exponentially() {
    let a = matrix(1, &[(0, 0, -1.0)]);
    for order in ORDERS {
        for step in -60..=50 {
            let x = 10f64.powf(f64::from(step) / 10.0);
            let got = deplete(&a, &[1.0], x, order).unwrap();
            let error = scaled_error(&got, &[(-x).exp()], &[1.0]);
            assert!(error <= REAL_SPECTRUM_TOL, "{order:?} x={x}: {error:e}");
        }
    }
}

#[test]
fn zero_time_keeps_the_densities() {
    let a = chain3(1e-3, 1e-5);
    let n0 = [2.0, 1.0, 0.5];
    for order in ORDERS {
        let got = deplete(&a, &n0, 0.0, order).unwrap();
        assert!(
            scaled_error(&got, &n0, &n0) <= REAL_SPECTRUM_TOL,
            "{order:?}: {got:?}"
        );
    }
}

#[test]
fn a_chain_matches_bateman() {
    // (λ1, λ2, t): 娘が速い、娘が遅い、λ‖A‖t が 1e9 の硬いチェーン
    let cases = [
        (1e-5, 1e-3, 1e4),
        (1e-5, 1e-3, 1e6),
        (1e-3, 1e-5, 1e4),
        (1e3, 1e-6, 1e6),
    ];
    for order in ORDERS {
        for (l1, l2, t) in cases {
            let got = deplete(&chain3(l1, l2), &[1.0, 0.0, 0.0], t, order).unwrap();
            let exact = parent_daughter_stable(l1, l2, t);
            let error = scaled_error(&got, &exact, &[1.0]);
            assert!(
                error <= REAL_SPECTRUM_TOL,
                "{order:?} {l1} {l2} {t}: {error:e}"
            );
        }
    }
}

#[test]
fn a_small_daughter_keeps_its_relative_accuracy() {
    // 長寿命の親が短寿命の娘を作り、娘は親の 1e-9 倍で釣り合う（Xe-135 のような核種）
    let (l1, l2, t) = (1e-12, 1e-3, 1e6);
    for order in ORDERS {
        let got = deplete(&chain3(l1, l2), &[1.0, 0.0, 0.0], t, order).unwrap();
        let exact = parent_daughter_stable(l1, l2, t);
        let relative = ((got[1] - exact[1]) / exact[1]).abs();
        assert!(relative <= SMALL_DAUGHTER_RTOL, "{order:?}: {relative:e}");
    }
}

#[test]
fn a_branching_decay_splits_by_the_ratios() {
    // 親 → 娘（割合 b、λ2 で安定な孫へ）と、安定な核種（割合 1 − b）
    let (l1, l2, b, t) = (1e-4, 1e-2, 0.3, 1e4);
    let a = matrix(
        4,
        &[
            (0, 0, -l1),
            (1, 0, b * l1),
            (3, 0, (1.0 - b) * l1),
            (1, 1, -l2),
            (2, 1, l2),
        ],
    );
    let [parent, daughter, stable] = parent_daughter_stable(l1, l2, t);
    let other = (1.0 - b) * -(-l1 * t).exp_m1();
    let exact = [parent, b * daughter, b * stable, other];
    for order in ORDERS {
        let got = deplete(&a, &[1.0, 0.0, 0.0, 0.0], t, order).unwrap();
        let error = scaled_error(&got, &exact, &[1.0]);
        assert!(error <= REAL_SPECTRUM_TOL, "{order:?}: {error:e}");
    }
}

#[test]
fn a_cycle_with_complex_eigenvalues_matches_the_closed_form() {
    // 0 → 1 → 2 → 0 を同じ割合 k で回る環。固有値は 0 と k(−3/2 ± i √3/2)。初めに 0 だけが 1 のとき
    // n_i = (1 + 2 e^{−3kt/2} cos(√3 kt/2 − 2πi/3)) / 3
    let a = matrix(
        3,
        &[
            (0, 0, -1.0),
            (1, 0, 1.0),
            (1, 1, -1.0),
            (2, 1, 1.0),
            (2, 2, -1.0),
            (0, 2, 1.0),
        ],
    );
    for (order, tol) in [
        (CramOrder::Order16, COMPLEX_SPECTRUM_TOL_16),
        (CramOrder::Order48, COMPLEX_SPECTRUM_TOL_48),
    ] {
        for kt in [0.1, 1.0, 10.0, 100.0] {
            let got = deplete(&a, &[1.0, 0.0, 0.0], kt, order).unwrap();
            let exact: Vec<f64> = (0..3)
                .map(|i| {
                    let phase =
                        3f64.sqrt() / 2.0 * kt - 2.0 * std::f64::consts::PI * f64::from(i) / 3.0;
                    (1.0 + 2.0 * (-1.5 * kt).exp() * phase.cos()) / 3.0
                })
                .collect();
            let error = scaled_error(&got, &exact, &[1.0]);
            assert!(error <= tol, "{order:?} kt={kt}: {error:e}");
        }
    }
}

#[test]
fn the_densities_scale_with_the_initial_densities() {
    // 2 の冪で拡大すると丸めが変わらないので、結果はちょうど同じ倍になる（単位を問わない）
    let a = chain3(1e-5, 1e-3);
    let scale = 2f64.powi(80);
    for order in ORDERS {
        let one = deplete(&a, &[1.0, 0.5, 0.0], 1e5, order).unwrap();
        let big = deplete(&a, &[scale, 0.5 * scale, 0.0], 1e5, order).unwrap();
        let scaled: Vec<f64> = one.iter().map(|v| v * scale).collect();
        assert_eq!(big, scaled, "{order:?}");
    }
}

#[test]
fn rejects_a_density_vector_of_another_length() {
    let a = chain3(1e-5, 1e-3);
    let error = deplete(&a, &[1.0, 0.0], 1.0, CramOrder::Order16).unwrap_err();
    assert!(
        matches!(
            error,
            DepletionError::LengthMismatch {
                expected: 3,
                found: 2
            }
        ),
        "{error}"
    );
}

#[test]
fn rejects_negative_and_non_finite_times() {
    let a = chain3(1e-5, 1e-3);
    for t in [-1.0, f64::NAN, f64::INFINITY] {
        let error = deplete(&a, &[1.0, 0.0, 0.0], t, CramOrder::Order48).unwrap_err();
        assert!(
            matches!(error, DepletionError::InvalidTime(_)),
            "{t}: {error}"
        );
    }
}

#[test]
fn rejects_non_finite_densities() {
    let a = chain3(1e-5, 1e-3);
    let error = deplete(&a, &[1.0, f64::NAN, 0.0], 1.0, CramOrder::Order48).unwrap_err();
    assert!(
        matches!(error, DepletionError::NonFiniteDensity { nuclide: 1, .. }),
        "{error}"
    );
}

#[test]
fn rejects_a_time_step_that_overflows_the_matrix() {
    let a = matrix(1, &[(0, 0, -1e300)]);
    let error = deplete(&a, &[1.0], 1e10, CramOrder::Order16).unwrap_err();
    assert!(
        matches!(error, DepletionError::Overflow { row: 0, col: 0 }),
        "{error}"
    );
}
