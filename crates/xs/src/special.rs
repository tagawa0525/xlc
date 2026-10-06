//! 特殊関数（標準ライブラリにない誤差関数）。

/// 誤差関数 erf(x)。
pub fn erf(x: f64) -> f64 {
    if x < 0.0 {
        return -erf(-x);
    }
    if x < SERIES_LIMIT {
        erf_series(x)
    } else {
        1.0 - erfc(x)
    }
}

/// これより小さい x では級数、大きい x では連分数を使う
const SERIES_LIMIT: f64 = 2.5;

/// erf(x) = 2/√π exp(-x²) Σ 2ⁿ x^(2n+1) / (1·3·…·(2n+1))。項がすべて正なので桁落ちしない。
fn erf_series(x: f64) -> f64 {
    let x2 = x * x;
    let mut term = x;
    let mut sum = x;
    let mut n = 0.0;
    while term > sum * f64::EPSILON {
        n += 1.0;
        term *= 2.0 * x2 / (2.0 * n + 1.0);
        sum += term;
    }
    2.0 / std::f64::consts::PI.sqrt() * (-x2).exp() * sum
}

/// 相補誤差関数 erfc(x) = 1 - erf(x)。
pub fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        2.0 - erfc(-x)
    } else if x < SERIES_LIMIT {
        1.0 - erf_series(x)
    } else {
        erfcx(x) * (-x * x).exp()
    }
}

/// スケールした相補誤差関数 erfcx(x) = exp(x²) erfc(x)（x ≥ 0）。
/// exp(x²) と erfc(x) を別に計算すると、大きな x で桁あふれ・桁落ちする。
pub fn erfcx(x: f64) -> f64 {
    if x < SERIES_LIMIT {
        return (x * x).exp() * (1.0 - erf_series(x));
    }
    // 連分数 erfcx(x) = 1/√π · 1/(x + (1/2)/(x + 1/(x + (3/2)/(x + …)))) を後ろから評価する。
    // x ≥ 2.5 では 40 項で倍精度の精度になる
    const TERMS: u32 = 40;
    let mut t = x;
    for k in (1..=TERMS).rev() {
        t = x + 0.5 * f64::from(k) / t;
    }
    1.0 / (std::f64::consts::PI.sqrt() * t)
}

/// スケールした指数積分 e^z E1(z)（z > 0）。E1(z) = ∫_z^∞ e^(-t)/t dt。
///
/// z ≤ 1 では級数 E1 = -γ - ln z - Σ (-z)^n/(n n!)、z > 1 では連分数
/// E1 = e^(-z) / (z + 1 - 1²/(z + 3 - 2²/(z + 5 - …)))（修正 Lentz 法）で計算する。
pub fn exp_e1(z: f64) -> f64 {
    assert!(z > 0.0, "E1 needs a positive argument, got {z}");
    const EULER: f64 = 0.577_215_664_901_532_9;
    const MAX_TERMS: u32 = 1000;
    if z <= 1.0 {
        let (mut sum, mut term) = (0.0, 1.0);
        for n in 1..=MAX_TERMS {
            let n = f64::from(n);
            term *= -z / n;
            sum += term / n;
            if (term / n).abs() <= f64::EPSILON * sum.abs() {
                return (-EULER - z.ln() - sum) * z.exp();
            }
        }
    } else {
        let tiny = f64::MIN_POSITIVE;
        let mut b = z + 1.0;
        let (mut c, mut d) = (1.0 / tiny, 1.0 / b);
        let mut h = d;
        for i in 1..=MAX_TERMS {
            let a = -f64::from(i * i);
            b += 2.0;
            d = 1.0 / (a * d + b);
            c = b + a / c;
            let delta = c * d;
            h *= delta;
            if (delta - 1.0).abs() <= f64::EPSILON {
                return h;
            }
        }
    }
    panic!("E1({z}) did not converge");
}

/// Faddeeva 関数 w(z) = e^(-z²) erfc(-iz) の (実部, 虚部)（z = x + iy, y ≥ 0）。
///
/// Weideman (1994, SIAM J. Numer. Anal. 31) の有理近似 w(z) ≈ 2 p(Z)/(L - iz)² + 1/(√π (L - iz))、
/// Z = (L + iz)/(L - iz) を使う。p は `FADDEEVA_TERMS` 次の多項式で、係数は `faddeeva_coefficients`。
pub fn faddeeva(x: f64, y: f64) -> (f64, f64) {
    assert!(y >= 0.0, "faddeeva is implemented for Im z >= 0, got {y}");
    let a = faddeeva_coefficients();
    let l = faddeeva_scale();
    // L - iz = (L + y) - ix、L + iz = (L - y) + ix
    let den = (l + y, -x);
    let z = complex_div((l - y, x), den);
    // Horner 法で p(Z) = Σ a_m Z^(m-1)
    let mut p = (0.0, 0.0);
    for &coefficient in a.iter().rev() {
        p = complex_mul(p, z);
        p.0 += coefficient;
    }
    let den2 = complex_mul(den, den);
    let first = complex_div((2.0 * p.0, 2.0 * p.1), den2);
    let second = complex_div((1.0 / std::f64::consts::PI.sqrt(), 0.0), den);
    (first.0 + second.0, first.1 + second.1)
}

/// Weideman の近似の項数 N。N = 32 で、テストの範囲の相対誤差が 1e-12 以下になる
const FADDEEVA_TERMS: usize = 32;

/// Weideman の近似の尺度 L = (N/√2)^(1/2)
fn faddeeva_scale() -> f64 {
    (FADDEEVA_TERMS as f64 / std::f64::consts::SQRT_2).sqrt()
}

/// p の係数 a_1..a_N。f(t) = e^(-t²)(L² + t²) を t_k = L tan(kπ/(2M))（M = 2N, |k| < M）で標本化し、
/// a_m = (1/2M) Σ_k f(t_k) cos(π m k / M)（Weideman の FFT を偶関数として直接計算したもの）。
fn faddeeva_coefficients() -> &'static [f64; FADDEEVA_TERMS] {
    static COEFFICIENTS: std::sync::OnceLock<[f64; FADDEEVA_TERMS]> = std::sync::OnceLock::new();
    COEFFICIENTS.get_or_init(|| {
        let m = 2 * FADDEEVA_TERMS;
        let l = faddeeva_scale();
        let samples: Vec<(f64, f64)> = (1 - m as i64..m as i64)
            .map(|k| {
                let theta = k as f64 * std::f64::consts::PI / m as f64;
                let t = l * (0.5 * theta).tan();
                (k as f64, (-t * t).exp() * (l * l + t * t))
            })
            .collect();
        std::array::from_fn(|i| {
            let order = (i + 1) as f64;
            let sum: f64 = samples
                .iter()
                .map(|&(k, f)| f * (std::f64::consts::PI * order * k / m as f64).cos())
                .sum();
            sum / (2 * m) as f64
        })
    })
}

fn complex_mul(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}

fn complex_div(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    let d = b.0 * b.0 + b.1 * b.1;
    ((a.0 * b.0 + a.1 * b.1) / d, (a.1 * b.0 - a.0 * b.1) / d)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_rel(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() <= tol * expected.abs(),
            "actual {actual}, expected {expected}"
        );
    }

    #[test]
    fn erf_values() {
        assert_eq!(erf(0.0), 0.0);
        assert_rel(erf(0.5), 0.520_499_877_813_046_5, 1e-14);
        assert_rel(erf(1.0), 0.842_700_792_949_714_9, 1e-14);
        assert_rel(erf(2.0), 0.995_322_265_018_952_7, 1e-14);
        assert_rel(erf(-1.0), -0.842_700_792_949_714_9, 1e-14);
        assert_rel(erf(5.0), 0.999_999_999_998_462_6, 1e-14);
    }

    #[test]
    fn erfc_values() {
        assert_rel(erfc(0.5), 0.479_500_122_186_953_5, 1e-13);
        assert_rel(erfc(2.0), 0.004_677_734_981_047_266, 1e-12);
        assert_rel(erfc(3.0), 2.209_049_699_858_544e-5, 1e-12);
        assert_rel(erfc(5.0), 1.537_459_794_428_035e-12, 1e-12);
        assert_rel(erfc(-1.0), 1.842_700_792_949_715, 1e-14);
    }

    #[test]
    fn erfcx_large_argument() {
        assert_rel(erfcx(10.0), 0.056_140_992_743_822_594, 1e-12);
        // 漸近形 1/(x√π)
        let x = 1.0e4;
        assert_rel(erfcx(x), 1.0 / (x * std::f64::consts::PI.sqrt()), 1e-8);
        assert_rel(erfcx(1.0), 1.0_f64.exp() * 0.157_299_207_050_285_13, 1e-13);
    }

    #[test]
    fn scaled_exponential_integral() {
        // decimal（60 桁）で級数を足して求めた値
        for (z, expected) in [
            (0.1, 2.014_642_544_708_451_7),
            (0.5, 0.922_910_632_483_730_5),
            (1.0, 0.596_347_362_323_194_1),
            (2.0, 0.361_328_616_888_222_6),
            (5.0, 0.170_422_176_284_732_2),
            (10.0, 0.091_563_333_939_788_08),
            (30.0, 0.032_289_738_758_980_13),
        ] {
            assert_rel(exp_e1(z), expected, 1e-13);
        }
        // 漸近形 1/z (1 - 1/z + 2/z² - ...)
        let z = 1.0e6;
        assert_rel(exp_e1(z), (1.0 - 1.0 / z + 2.0 / (z * z)) / z, 1e-15);
    }

    /// Re w(x + iy) = (y/π) ∫ e^(-t²) / ((x - t)² + y²) dt（y > 0）を Simpson 則で積分した値。
    fn re_w_by_integration(x: f64, y: f64) -> f64 {
        // Lorentz 形の幅 y と Gauss 形の幅 1 を十分細かく分ける
        let h = y.min(1.0) / 200.0;
        let (lo, hi) = ((-9.0_f64).min(x - 60.0 * y), 9.0_f64.max(x + 60.0 * y));
        let n = (((hi - lo) / h).ceil() as usize).next_multiple_of(2);
        let h = (hi - lo) / n as f64;
        let f = |t: f64| (-t * t).exp() / ((x - t) * (x - t) + y * y);
        let mut sum = f(lo) + f(hi);
        for i in 1..n {
            sum += if i % 2 == 1 { 4.0 } else { 2.0 } * f(lo + h * i as f64);
        }
        sum * h / 3.0 * y / std::f64::consts::PI
    }

    #[test]
    fn faddeeva_real_part_matches_definition() {
        for (x, y) in [(0.5, 0.05), (3.0, 0.1), (1.0, 1.0), (10.0, 0.2), (2.0, 5.0)] {
            assert_rel(faddeeva(x, y).0, re_w_by_integration(x, y), 1e-10);
        }
    }

    #[test]
    fn faddeeva_on_imaginary_axis_is_erfcx() {
        for y in [1.0e-3, 0.1, 1.0, 10.0, 100.0] {
            let (re, im) = faddeeva(0.0, y);
            assert_rel(re, erfcx(y), 1e-12);
            assert!(im.abs() <= 1e-15, "{im}");
        }
    }

    #[test]
    fn faddeeva_on_real_axis_is_gaussian() {
        for x in [0.0, 0.5, 2.0, 5.0] {
            let re = faddeeva(x, 0.0).0;
            assert!((re - (-x * x).exp()).abs() <= 1e-13, "x={x}: {re}");
        }
    }

    #[test]
    fn faddeeva_far_from_origin_is_asymptotic() {
        // w(z) ≈ i/(√π z) (1 + 1/(2z²) + 3/(4z⁴))
        let (x, y) = (300.0, 0.07);
        let r2 = x * x + y * y;
        let sqrt_pi = std::f64::consts::PI.sqrt();
        // i/z = (y + ix)/|z|²。実部は y/(√π |z|²) に高次の補正が付く
        let (re, im) = faddeeva(x, y);
        assert_rel(re, y / (sqrt_pi * r2) * (1.0 + 3.0 / (2.0 * x * x)), 1e-6);
        assert_rel(im, x / (sqrt_pi * r2) * (1.0 + 1.0 / (2.0 * x * x)), 1e-9);
    }
}
