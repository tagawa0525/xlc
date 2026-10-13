use super::*;

#[test]
fn continuous_at_breakpoints() {
    let w = Weight::lwr();
    for e in w.breakpoints() {
        let (below, above) = (w.value(e * (1.0 - 1e-12)), w.value(e * (1.0 + 1e-12)));
        assert!((below - above).abs() <= 1e-9 * above, "at {e}");
    }
    assert_eq!(w.value(1.0), 1.0);
}

#[test]
fn shapes() {
    let w = Weight::lwr();
    let r = w.value(0.02) / w.value(0.01);
    assert!((r - 2.0 * (-0.01_f64 / 0.0253).exp()).abs() < 1e-12);
    let r = w.value(4.0e6) / w.value(2.0e6);
    assert!((r - 2.0_f64.sqrt() * (-2.0e6_f64 / 1.4e6).exp()).abs() < 1e-12);
}
