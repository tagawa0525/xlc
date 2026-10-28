use super::*;
use xlc_numeric::error_function::erfc;

/// エラーにならない関数の期待値
fn mean(dof: WidthDof, tolerance: f64, f: impl Fn(f64) -> f64) -> f64 {
    expectation::<NumericError>(dof, tolerance, |x| Ok(f(x))).unwrap()
}

const SQRT_PI: f64 = 1.772_453_850_905_516;
const INNER: f64 = 1.0e-10;
const OUTER: f64 = 1.0e-8;

#[test]
fn reproduces_moments_of_chi_square() {
    // 平均 1 の χ²_ν の k 次の積率は Π_{i<k} (1 + 2i/ν)
    for dof in [
        WidthDof::One,
        WidthDof::Two,
        WidthDof::Three,
        WidthDof::Four,
    ] {
        let nu = f64::from(dof.number());
        let mut expected = 1.0;
        for k in 0..=6 {
            let moment = mean(dof, INNER, |x| x.powi(k));
            assert!(
                (moment - expected).abs() <= 1e-9 * expected,
                "ν={dof:?}, k={k}: {moment} vs {expected}"
            );
            expected *= 1.0 + 2.0 * f64::from(k) / nu;
        }
    }
}

#[test]
fn matches_closed_form_near_the_pole() {
    // E[1/(x+c)] = k^k c^(k-1) e^(kc) Γ(1-k, kc)。
    // Γ(1/2, z) = √π erfc(√z)、Γ(-1/2, z) = 2 e^(-z)/√z - 2√π erfc(√z)
    for c in [1.0e-4_f64, 1.0e-2, 1.0, 10.0] {
        let z = 0.5 * c;
        let nu1 = 0.5_f64.sqrt() / c.sqrt() * z.exp() * SQRT_PI * erfc(z.sqrt());
        let z = 1.5 * c;
        let gamma = 2.0 * (-z).exp() / z.sqrt() - 2.0 * SQRT_PI * erfc(z.sqrt());
        let nu3 = 1.5_f64.powf(1.5) * c.sqrt() * z.exp() * gamma;
        for (dof, exact) in [(WidthDof::One, nu1), (WidthDof::Three, nu3)] {
            let actual = mean(dof, INNER, |x| 1.0 / (x + c));
            assert!(
                (actual - exact).abs() <= 1e-8 * exact,
                "ν={dof:?}, c={c}: {actual} vs {exact}"
            );
        }
    }
}

#[test]
fn nested_expectations_converge_in_either_order() {
    // 内側の期待値は許容誤差の範囲でばらつくので、外側は許容誤差を緩めないと収束しない
    let f = |x: f64, y: f64| x * y / (x + 0.3 * y + 0.01);
    let x_inside = mean(WidthDof::One, OUTER, |y| {
        mean(WidthDof::One, INNER, |x| f(x, y))
    });
    let y_inside = mean(WidthDof::One, OUTER, |x| {
        mean(WidthDof::One, INNER, |y| f(x, y))
    });
    assert!(
        (x_inside - y_inside).abs() <= 1e-7 * x_inside,
        "{x_inside} vs {y_inside}"
    );
}

#[test]
fn width_ratios_match_integration() {
    // c が大きいところ（Γn ≪ Γγ、たとえば低エネルギーの l=2）まで確かめる
    for dof in [
        WidthDof::Fixed,
        WidthDof::One,
        WidthDof::Two,
        WidthDof::Three,
        WidthDof::Four,
    ] {
        for c in [
            1.0e-6, 1.0e-3, 0.1, 1.0, 10.0, 99.0, 101.0, 1.0e3, 1.0e6, 1.0e9,
        ] {
            let numeric = [
                mean(dof, 1.0e-12, |x| x / (x + c)),
                mean(dof, 1.0e-12, |x| x * x / (x + c)),
            ];
            let closed = width_ratios(dof, c).unwrap();
            for i in 0..2 {
                assert!(
                    (closed[i] - numeric[i]).abs() <= 1e-9 * numeric[i],
                    "ν={dof:?}, c={c}, [{i}]: {} vs {}",
                    closed[i],
                    numeric[i]
                );
            }
        }
    }
}

#[test]
fn zero_degrees_of_freedom_do_not_fluctuate() {
    assert_eq!(mean(WidthDof::Fixed, INNER, |x| 3.0 * x + 1.0), 4.0);
}
