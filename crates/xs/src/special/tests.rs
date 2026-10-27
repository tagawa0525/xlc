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
    // まとめた計算は、点ごとの `faddeeva` と同じ演算なので、ビットまで同じ。内部の区切りをまたぐ長さ、空の列、下半面も確かめる
    for y in [-0.5, 0.0, 1.0e-3, 0.07, 2.0, 99.0] {
        for n in [0, 1, 7, 32, 33, 250] {
            let x: Vec<f64> = (0..n).map(|i| -120.0 + 240.0 * i as f64 / 249.0).collect();
            let mut w = vec![Complex64::new(7.0, 7.0)];
            faddeeva_along_x(y, &x, &mut w);
            assert_eq!(w.len(), n + 1);
            assert_eq!(w[0], Complex64::new(7.0, 7.0));
            for (&x, w) in x.iter().zip(&w[1..]) {
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

/// erfcx(x) の多倍長の値（x, hi, lo。hi + lo が値）。mpmath 1.4.1 の 60 桁で erfc(x)·exp(x²) を求め、hi に倍精度へ丸めた値、
/// lo に残りを丸めた値を置いた。x = 1e200 は mpmath の erfc が扱えないので、漸近展開の 2 項 (1 − 1/(2x²))/(√π x)（次の項は
/// 1e−800 の相対の大きさ）。x は区間の端と、2.5 ≤ x の一様な乱数の 10 万点で誤差が大きかった点を含む
const ERFCX_REFERENCE: &[(f64, f64, f64)] = &[
    (2.5, 0.2108063640611436, -5.6277259093102524e-18),
    (
        2.5124651723206233,
        0.20988347101234978,
        1.3445554062764545e-17,
    ),
    (
        2.5846515745793606,
        0.20468671679751543,
        -4.829544975244344e-18,
    ),
    (
        2.6050776195927656,
        0.20326047815980136,
        -2.2696291215348346e-18,
    ),
    (2.75, 0.1936620962790687, -1.2015846532739174e-17),
    (3.0, 0.17900115118138996, -5.4272175920200274e-18),
    (3.5, 0.1552936556088943, -1.355844542216092e-18),
    (
        3.9999999999999996,
        0.1369994576250614,
        -6.177841842557427e-18,
    ),
    (4.0, 0.13699945762506138, 7.196568139158719e-18),
    (4.000000000000001, 0.13699945762506135, 6.18981248696211e-18),
    (
        4.44892545588937,
        0.12382782538384467,
        6.6209316208332245e-18,
    ),
    (
        4.463089702624895,
        0.12345249057391196,
        1.6302355169082588e-18,
    ),
    (5.0, 0.11070463773306863, -1.832347493639739e-18),
    (6.0, 0.09277656780053835, 6.215364755528485e-18),
    (8.0, 0.06998516620088092, 3.2863406596468746e-18),
    (
        9.57370480820499,
        0.05861480779876944,
        2.9943406702944013e-18,
    ),
    (10.0, 0.05614099274382259, -1.6720611399896374e-18),
    (20.0, 0.02817434874105132, -4.829674968674034e-20),
    (
        36.82194217585145,
        0.015316458700546182,
        5.269510661225192e-20,
    ),
    (50.0, 0.011281536265323773, -8.629279370172108e-19),
    (100.0, 0.005641613782989433, -1.186182619499244e-19),
    (
        602.8467785386605,
        0.0009358743007719356,
        4.2957790157627714e-20,
    ),
    (1000.0, 0.0005641893014533876, 2.3118263082589717e-20),
    (10000.0, 5.641895807268084e-05, -1.484202149398048e-21),
    (1000000.0, 5.641895835474742e-07, 1.0958100640939896e-23),
    (100000000.0, 5.641895835477562e-09, 3.2112571487601016e-25),
    (
        1000000000000.0,
        5.641895835477563e-13,
        3.776500696651908e-30,
    ),
    (1e+100, 5.6418958354775625e-101, 2.76431524922385e-117),
    (1e200, 5.641895835477563e-201, 3.5981468851464104e-217),
];

#[test]
fn erfcx_is_within_two_ulp_of_multiprecision_values() {
    for &(x, hi, lo) in ERFCX_REFERENCE {
        let ulp = f64::from_bits(hi.to_bits() + 1) - hi;
        let error = ((erfcx(x) - hi) - lo) / ulp;
        assert!(error.abs() <= 2.0, "x={x}: {error} ulp");
    }
}

/// erfcx(x) のビット。lattice の `ufg::free_gas::tests` の同じ名前のテストに同じ表を置き、2 つのクレートの erfcx が同じ値を返すことを確かめる
/// （2 つは同じ計算の複製。ADR-0203）
const ERFCX_BITS: &[(f64, u64)] = &[
    (0.0, 0x3ff0000000000000),
    (0.5, 0x3fe3b3bc3c98b0f1),
    (1.0, 0x3fdb5d8780f956ab),
    (2.0, 0x3fd058671b52c866),
    (2.4999999999999996, 0x3fcafbb3f3b72f80),
    (2.5, 0x3fcafbb3f3b7343b),
    (2.5124651723206233, 0x3fcadd7629fc0866),
    (2.5846515745793606, 0x3fca332ca14914a4),
    (2.6050776195927656, 0x3fca047079220171),
    (2.75, 0x3fc8c9eb68ff27d6),
    (3.0, 0x3fc6e9827d229d2d),
    (3.5, 0x3fc3e0a99a0ee914),
    (3.9999999999999996, 0x3fc18932bf08e155),
    (4.0, 0x3fc18932bf08e154),
    (4.000000000000001, 0x3fc18932bf08e153),
    (4.44892545588937, 0x3fbfb32e2c5bc0da),
    (4.463089702624895, 0x3fbf9a95199fec59),
    (5.0, 0x3fbc57239e943d1a),
    (6.0, 0x3fb7c0348489d721),
    (8.0, 0x3fb1ea8c4009b459),
    (9.57370480820499, 0x3fae02c2951d3468),
    (10.0, 0x3facbe831f997123),
    (20.0, 0x3f9cd9bc89b73548),
    (36.82194217585145, 0x3f8f5e3c49accfb1),
    (50.0, 0x3f871ac62a773371),
    (100.0, 0x3f771ba92b1d1777),
    (602.8467785386605, 0x3f4eaaaec1eaf0e6),
    (1000.0, 0x3f427cc34c9b7f2d),
    (10000.0, 0x3f0d946ca36c60e7),
    (1000000.0, 0x3ea2ee5a03c7620f),
    (100000000.0, 0x3e383b547fb79514),
    (1000000000000.0, 0x3d63d9c49b6a903e),
    (1e100, 0x2b1f974cd8b00f7d),
    (1e200, 0x165ba394ce53f796),
];

#[test]
fn erfcx_gives_the_same_bits_as_its_copy() {
    for &(x, bits) in ERFCX_BITS {
        assert_eq!(erfcx(x).to_bits(), bits, "x={x}");
    }
}
