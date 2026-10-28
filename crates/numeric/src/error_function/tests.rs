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

mod cases;
use cases::*;

#[test]
fn erfcx_is_within_two_and_a_half_ulp_of_multiprecision_values() {
    for &(x, hi, lo) in ERFCX_REFERENCE {
        let error = ulp_error(erfcx(x), hi, lo);
        assert!(error.abs() <= 2.5, "x={x}: {error} ulp");
    }
}

#[test]
fn erfcx_of_negative_x_is_within_four_ulp_of_multiprecision_values() {
    for &(x, hi, lo) in ERFCX_NEGATIVE_REFERENCE {
        let error = ulp_error(erfcx(x), hi, lo);
        assert!(error.abs() <= 4.0, "x={x}: {error} ulp");
    }
}

#[test]
fn erf_is_within_four_ulp_of_multiprecision_values() {
    for &(x, hi, lo) in ERF_REFERENCE {
        let error = ulp_error(erf(x), hi, lo);
        assert!(error.abs() <= 4.0, "x={x}: {error} ulp");
    }
}

#[test]
fn erfc_and_erfcx_saturate_at_large_arguments() {
    assert_eq!(erfc(f64::INFINITY), 0.0);
    assert_eq!(erfc(1e308), 0.0);
    assert_eq!(erfc(f64::NEG_INFINITY), 2.0);
    assert_eq!(erfcx(f64::NEG_INFINITY), f64::INFINITY);
    assert_eq!(erfcx(-1e308), f64::INFINITY);
    assert_eq!(erf(1e308), 1.0);
}

#[test]
fn erf_and_erfcx_of_nan_are_nan() {
    assert!(erf(f64::NAN).is_nan());
    assert!(erfcx(f64::NAN).is_nan());
}

#[test]
fn erfc_is_within_four_ulp_of_multiprecision_values() {
    for &(x, hi, lo) in ERFC_REFERENCE {
        let error = ulp_error(erfc(x), hi, lo);
        assert!(error.abs() <= 4.0, "x={x}: {error} ulp");
    }
}
