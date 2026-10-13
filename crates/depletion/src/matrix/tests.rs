use crate::{CramOrder, DepletionError, DepletionMatrix, deplete};

#[test]
fn sums_entries_at_the_same_position() {
    let split = DepletionMatrix::from_triplets(
        2,
        &[(0, 0, -1e-3), (0, 0, -2e-3), (1, 0, 1e-3), (1, 0, 2e-3)],
    )
    .unwrap();
    let whole = DepletionMatrix::from_triplets(2, &[(0, 0, -3e-3), (1, 0, 3e-3)]).unwrap();
    let a = deplete(&split, &[1.0, 0.0], 100.0, CramOrder::Order16).unwrap();
    let b = deplete(&whole, &[1.0, 0.0], 100.0, CramOrder::Order16).unwrap();
    assert_eq!(a, b);
}

#[test]
fn keeps_a_nuclide_without_entries_constant() {
    // 対角を書かない核種（安定で反応もない）も、行列は対角を持つ
    let a = DepletionMatrix::from_triplets(2, &[(0, 0, -1e-3)]).unwrap();
    let got = deplete(&a, &[1.0, 2.0], 100.0, CramOrder::Order48).unwrap();
    assert!((got[1] - 2.0).abs() <= 1e-13, "{got:?}");
}

#[test]
fn rejects_negative_production() {
    let error = DepletionMatrix::from_triplets(2, &[(0, 0, -1.0), (1, 0, -1.0)]).unwrap_err();
    assert!(
        matches!(
            error,
            DepletionError::NegativeProduction { row: 1, col: 0, .. }
        ),
        "{error}"
    );
}

#[test]
fn rejects_a_positive_diagonal() {
    let error = DepletionMatrix::from_triplets(2, &[(1, 1, 1.0)]).unwrap_err();
    assert!(
        matches!(error, DepletionError::PositiveLoss { nuclide: 1, .. }),
        "{error}"
    );
}

#[test]
fn judges_the_signs_after_summing() {
    // 同じ位置の負の寄与と、それより大きい正の寄与は、足して正なのでよい
    DepletionMatrix::from_triplets(2, &[(0, 0, -1.0), (1, 0, 2.0), (1, 0, -1.0)]).unwrap();
}

#[test]
fn rejects_entries_outside_the_matrix() {
    let error = DepletionMatrix::from_triplets(2, &[(2, 0, 1.0)]).unwrap_err();
    assert!(
        matches!(
            error,
            DepletionError::OutOfRange {
                row: 2,
                col: 0,
                n: 2
            }
        ),
        "{error}"
    );
}

#[test]
fn rejects_non_finite_entries() {
    let error = DepletionMatrix::from_triplets(2, &[(1, 0, f64::INFINITY)]).unwrap_err();
    assert!(
        matches!(error, DepletionError::NonFiniteEntry { row: 1, col: 0, .. }),
        "{error}"
    );
}
