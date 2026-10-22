use super::*;

fn assert_rel(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol * expected.abs(),
        "actual {actual}, expected {expected}"
    );
}

#[test]
fn faddeeva_below_the_real_axis() {
    // 虚軸の上では w(iy) = e^(y²) erfc(y)。y < 0 でも成り立つ
    for y in [-0.1, -0.5, -2.0] {
        let Complex64 { re, im } = faddeeva(Complex64::new(0.0, y));
        assert_rel(re, (y * y).exp() * erfc(y), 1e-10);
        assert!(im.abs() < 1e-12, "{im}");
    }
}

#[test]
fn e1_outside_its_domain_is_error() {
    for z in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(exp_e1(z), Err(NumericError::Domain(_))), "E1({z})");
    }
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
        assert_rel(exp_e1(z).unwrap(), expected, 1e-13);
    }
    // 漸近形 1/z (1 - 1/z + 2/z² - ...)
    let z = 1.0e6;
    assert_rel(
        exp_e1(z).unwrap(),
        (1.0 - 1.0 / z + 2.0 / (z * z)) / z,
        1e-15,
    );
}

#[test]
fn regularized_incomplete_gamma() {
    // decimal（50 桁）で級数 P(a, x) = x^a e^(-x)/Γ(a+1) Σ xⁿ/((a+1)…(a+n)) を足して求めた値
    for (x, p32, p2) in [
        (1.0e-8, 7.522_527_735_501_584e-13, 4.999_999_966_666_667e-17),
        (0.01, 7.477_553_393_911_979e-4, 4.966_791_334_026_589e-5),
        (0.5, 0.198_748_043_098_799_2, 0.090_204_010_431_049_86),
        (1.0, 0.427_593_295_529_120_2, 0.264_241_117_657_115_3),
        (3.0, 0.888_389_774_905_287_5, 0.800_851_726_528_544_2),
        (20.0, 0.999_999_989_344_909_6, 0.999_999_956_715_773_9),
    ] {
        assert_rel(gamma_p_3_2(x), p32, 1e-14);
        assert_rel(gamma_p_2(x), p2, 1e-14);
    }
    assert_eq!((gamma_p_3_2(0.0), gamma_p_2(0.0)), (0.0, 0.0));
    assert_eq!((gamma_p_3_2(1.0e3), gamma_p_2(1.0e3)), (1.0, 1.0));
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
        assert_rel(
            faddeeva(Complex64::new(x, y)).re,
            re_w_by_integration(x, y),
            1e-10,
        );
    }
}

#[test]
fn faddeeva_on_imaginary_axis_is_erfcx() {
    for y in [1.0e-3, 0.1, 1.0, 10.0, 100.0] {
        let Complex64 { re, im } = faddeeva(Complex64::new(0.0, y));
        assert_rel(re, erfcx(y), 1e-12);
        assert!(im.abs() <= 1e-15, "{im}");
    }
}

#[test]
fn faddeeva_on_real_axis_is_gaussian() {
    for x in [0.0, 0.5, 2.0, 5.0] {
        let re = faddeeva(Complex64::new(x, 0.0)).re;
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
    let Complex64 { re, im } = faddeeva(Complex64::new(x, y));
    assert_rel(re, y / (sqrt_pi * r2) * (1.0 + 3.0 / (2.0 * x * x)), 1e-6);
    assert_rel(im, x / (sqrt_pi * r2) * (1.0 + 1.0 / (2.0 * x * x)), 1e-9);
}

#[test]
fn faddeeva_over_many_x_gives_the_bits_of_each_point() {
    // まとめた計算は、点ごとの `faddeeva` と同じ演算なので、ビットまで同じ。内部の区切りをまたぐ長さと、空の列も確かめる
    for y in [0.0, 1.0e-3, 0.07, 2.0, 99.0] {
        for n in [0, 1, 7, 64, 65, 250] {
            let x: Vec<f64> = (0..n).map(|i| -120.0 + 240.0 * i as f64 / 249.0).collect();
            let mut w = vec![Complex64::ZERO; n];
            faddeeva_upper(y, &x, &mut w);
            for (&x, w) in x.iter().zip(&w) {
                let expected = faddeeva(Complex64::new(x, y));
                assert_eq!(
                    (w.re.to_bits(), w.im.to_bits()),
                    (expected.re.to_bits(), expected.im.to_bits()),
                    "x={x}, y={y}"
                );
            }
        }
    }
}
