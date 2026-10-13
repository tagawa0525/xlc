use super::*;

#[test]
fn legendre_matches_closed_forms() {
    let x = 0.3_f64;
    let p = legendre(x, 3);
    assert_eq!(p[0], 1.0);
    assert_eq!(p[1], x);
    assert!((p[2] - 0.5 * (3.0 * x * x - 1.0)).abs() < 1e-15);
    assert!((p[3] - 0.5 * (5.0 * x * x * x - 3.0 * x)).abs() < 1e-15);
}

#[test]
fn tables_match_the_computed_nodes() {
    for table in [
        &GAUSS_LEGENDRE_2[..],
        &GAUSS_LEGENDRE_4[..],
        &GAUSS_LEGENDRE_8[..],
    ] {
        let mut computed = gauss_legendre(table.len()).unwrap();
        computed.sort_by(|a, b| a.0.total_cmp(&b.0));
        for ((x, w), (cx, cw)) in table.iter().zip(&computed) {
            assert!((x - cx).abs() < 1e-15 && (w - cw).abs() < 1e-15);
        }
    }
}

#[test]
fn zero_points_is_an_error() {
    assert!(gauss_legendre(0).is_err());
}

#[test]
fn nodes_integrate_polynomials_exactly() {
    // n 点は 2n − 1 次までの多項式を正確に積分する
    for n in [1, 3, 7, 20] {
        let nodes = gauss_legendre(n).unwrap();
        let degree = 2 * n - 2;
        let integral: f64 = nodes.iter().map(|&(x, w)| w * x.powi(degree as i32)).sum();
        assert!((integral - 2.0 / (degree as f64 + 1.0)).abs() < 1e-13);
    }
}
