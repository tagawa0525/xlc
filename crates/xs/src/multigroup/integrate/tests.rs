use super::*;

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-10 * expected.abs(),
        "actual {actual}, expected {expected}"
    );
}

#[test]
fn weight_only_in_one_over_e_region() {
    let groups = GroupStructure::new(vec![1.0e5, 1.0e3, 1.0]).unwrap();
    let r = group_integrals(&groups, &Weight::lwr(), &[], |_| 1.0);
    assert_close(r[0], (1.0e5_f64 / 1.0e3).ln());
    assert_close(r[1], (1.0e3_f64).ln());
}

#[test]
fn weight_only_in_maxwellian_region() {
    let w = Weight::lwr();
    let (a, b, tb) = (0.01, 0.05, 0.0253);
    let groups = GroupStructure::new(vec![b, a]).unwrap();
    let r = group_integrals(&groups, &w, &[], |_| 1.0);
    // ∫ ab E exp(-E/tb) dE = ab tb² [(1 + a/tb) e^{-a/tb} - (1 + b/tb) e^{-b/tb}]
    let ab = w.value(a) / (a * (-a / tb).exp());
    let primitive = |e: f64| (1.0 + e / tb) * (-e / tb).exp();
    assert_close(r[0], ab * tb * tb * (primitive(a) - primitive(b)));
}

#[test]
fn linear_function_times_one_over_e() {
    let groups = GroupStructure::new(vec![100.0, 10.0]).unwrap();
    let r = group_integrals(&groups, &Weight::lwr(), &[], |e| e);
    assert_close(r[0], 90.0);
}

#[test]
fn step_at_breakpoint() {
    let groups = GroupStructure::new(vec![100.0, 10.0]).unwrap();
    let r = group_integrals(&groups, &Weight::lwr(), &[30.0], |e| {
        if e < 30.0 { 1.0 } else { 2.0 }
    });
    assert_close(r[0], (3.0_f64).ln() + 2.0 * (100.0_f64 / 30.0).ln());
}

#[test]
fn multi_component_matches_scalar() {
    let groups = GroupStructure::new(vec![100.0, 10.0, 1.0]).unwrap();
    let w = Weight::lwr();
    let r = group_integrals_multi(&groups, &w, &[30.0], 2, |e, out| {
        out[0] = e;
        out[1] = 1.0;
    });
    let e = group_integrals(&groups, &w, &[30.0], |e| e);
    let one = group_integrals(&groups, &w, &[30.0], |_| 1.0);
    for g in 0..2 {
        assert_close(r[g][0], e[g]);
        assert_close(r[g][1], one[g]);
    }
}

#[test]
fn group_spanning_weight_breakpoint() {
    let w = Weight::lwr();
    let groups = GroupStructure::new(vec![1.0, 0.05]).unwrap();
    let r = group_integrals(&groups, &w, &[], |_| 1.0);
    let tb = 0.0253;
    let ab = w.value(0.05) / (0.05 * (-0.05_f64 / tb).exp());
    let primitive = |e: f64| (1.0 + e / tb) * (-e / tb).exp();
    let expected = ab * tb * tb * (primitive(0.05) - primitive(0.1)) + (1.0_f64 / 0.1).ln();
    assert_close(r[0], expected);
}
