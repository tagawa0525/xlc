use super::*;

#[test]
fn nan_integrand_does_not_converge() {
    let result = gauss_kronrod::<NumericError>(
        |_, out| {
            out[0] = f64::NAN;
            Ok(())
        },
        1,
        0.0,
        1.0,
        1e-10,
    );
    assert!(matches!(result, Err(NumericError::NotConverged(_))));
}
