use super::*;

#[test]
fn erf_is_odd_and_erfc_is_two_minus_its_reflection() {
    assert_eq!(erf(0.0), 0.0);
    assert_eq!(erf(-1.0), -erf(1.0));
    assert_eq!(erfc(-1.0), 2.0 - erfc(1.0));
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
fn erf_with_a_given_gaussian_is_erf_bit_for_bit() {
    let samples = ERF_REFERENCE.iter().map(|&(x, _, _)| x).chain([
        0.0,
        CODY_LOWER_LIMIT,
        SATURATION,
        1e308,
        f64::INFINITY,
    ]);
    for x in samples.flat_map(|x| [x, -x]) {
        let gaussian = (-x * x).exp();
        assert_eq!(
            erf_with_gaussian(x, gaussian).to_bits(),
            erf(x).to_bits(),
            "x={x}"
        );
    }
    assert!(erf_with_gaussian(f64::NAN, f64::NAN).is_nan());
}

#[test]
fn erfc_and_erfcx_saturate_at_large_arguments() {
    assert_eq!(erfc(f64::INFINITY), 0.0);
    assert_eq!(erfcx(f64::INFINITY), 0.0);
    assert_eq!(erfcx(1e308), 0.5 * std::f64::consts::FRAC_2_SQRT_PI / 1e308);
    assert_eq!(erf(f64::NEG_INFINITY), -1.0);
    assert_eq!(erfc(-1e308), 2.0);
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
