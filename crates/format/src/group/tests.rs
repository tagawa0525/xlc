use super::*;

#[test]
fn wims69_bounds() {
    let g = GroupStructure::wims69();
    assert_eq!(g.len(), 69);
    assert_eq!(g.bounds()[0], 1.0e7);
    assert_eq!(g.bounds()[69], 1.0e-5);
    assert_eq!(g.range(0), (6.0655e6, 1.0e7));
    assert_eq!(g.range(68), (1.0e-5, 5.0e-3));
    // 熱群の境界 0.625 eV と 4 eV を含む
    assert!(g.bounds().contains(&0.625));
    assert!(g.bounds().contains(&4.0));
}

#[test]
fn find_group() {
    let g = GroupStructure::new(vec![10.0, 5.0, 1.0]).unwrap();
    assert_eq!(g.find(7.0), Some(0));
    assert_eq!(g.find(5.0), Some(1));
    assert_eq!(g.find(10.0), Some(0));
    assert_eq!(g.find(1.0), Some(1));
    assert_eq!(g.find(0.5), None);
    assert_eq!(g.find(11.0), None);
}

#[test]
fn invalid_bounds_are_errors() {
    // 積分区間を対数で分割するので、0 eV の境界は扱えない
    for bounds in [
        vec![1.0, 0.0],
        vec![1.0, 5.0],
        vec![1.0],
        vec![f64::INFINITY, 1.0],
        vec![1.0, f64::NAN],
    ] {
        assert!(GroupStructure::new(bounds.clone()).is_err(), "{bounds:?}");
    }
}

#[test]
fn wims69_passes_the_checks() {
    let bounds = GroupStructure::wims69().bounds().to_vec();
    assert!(GroupStructure::new(bounds).is_ok());
}
