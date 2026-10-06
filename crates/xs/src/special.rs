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
pub fn exp_e1(z: f64) -> f64 {
    let _ = z;
    todo!()
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
}
