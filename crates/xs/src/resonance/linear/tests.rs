use super::*;

fn c(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}

#[test]
fn solves_a_system_that_needs_pivoting() {
    // 先頭の対角が 0 なので行の入れ替えが要る。右辺は2列
    let m = [
        c(0.0, 0.0),
        c(2.0, 1.0),
        c(1.0, 0.0),
        c(1.0, -1.0),
        c(0.5, 0.0),
        c(0.0, 3.0),
        c(4.0, 0.0),
        c(-1.0, 2.0),
        c(2.0, 2.0),
    ];
    let b = [
        c(1.0, 0.0),
        c(0.0, 1.0),
        c(2.0, -1.0),
        c(3.0, 0.0),
        c(-1.0, 1.0),
        c(0.0, 0.5),
        c(1.0, 1.0),
        c(2.0, 0.0),
        c(0.0, -2.0),
    ];
    let mut work = m;
    let mut x = b;
    solve(&mut work, &mut x, 3).unwrap();
    for i in 0..3 {
        for j in 0..3 {
            let mx: Complex64 = (0..3).map(|k| m[i * 3 + k] * x[k * 3 + j]).sum();
            assert!((mx - b[i * 3 + j]).norm() < 1e-14, "({i}, {j}): {mx}");
        }
    }
}

#[test]
fn rejects_a_zero_pivot() {
    let mut m = [c(1.0, 0.0), c(2.0, 0.0), c(2.0, 0.0), c(4.0, 0.0)];
    let mut b = [c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(1.0, 0.0)];
    assert_eq!(solve(&mut m, &mut b, 2), Err(Singular));
}

#[test]
fn rejects_a_nan_pivot() {
    let mut m = [c(f64::NAN, 0.0), c(1.0, 0.0), c(1.0, 0.0), c(1.0, 0.0)];
    let mut b = [c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(1.0, 0.0)];
    assert_eq!(solve(&mut m, &mut b, 2), Err(Singular));
}
