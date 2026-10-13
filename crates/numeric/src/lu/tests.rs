use super::*;

#[test]
fn lu_solves_with_pivoting() {
    // 先頭の対角が 0 なので、行の入れ替えが要る
    let a = vec![0.0, 2.0, 1.0, 1.0, 1.0, 0.0, 3.0, 0.0, 1.0];
    let x = [1.0, -2.0, 0.5];
    let mut b: Vec<f64> = (0..3)
        .map(|i| (0..3).map(|j| a[i * 3 + j] * x[j]).sum())
        .collect();
    Lu::new(a, 3).unwrap().solve(&mut b);
    for (got, want) in b.iter().zip(x) {
        assert!((got - want).abs() < 1e-14, "{b:?}");
    }
    assert!(Lu::new(vec![1.0, 2.0, 2.0, 4.0], 2).is_none());
}

#[test]
fn solve_rejects_singular_and_overflowing_systems() {
    assert_eq!(solve(vec![1.0, 2.0, 2.0, 4.0], vec![1.0, 2.0]), None);
    // ピボットは 0 でないが、解があふれる
    assert_eq!(
        solve(vec![1.0e-300, 0.0, 0.0, 1.0], vec![1.0e300, 1.0]),
        None
    );
}

#[test]
fn nan_pivot_is_singular() {
    // total_cmp では NaN がどの数より大きいので、ピボットに選ばれてしまう
    assert!(Lu::new(vec![f64::NAN, 1.0, 1.0, 1.0], 2).is_none());
}
