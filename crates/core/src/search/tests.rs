use super::*;

fn values(f: impl Fn(f64) -> f64) -> impl FnMut(f64) -> Result<(f64, f64), String> {
    move |x| Ok((f(x), x))
}

/// k = 1 / (0.8 + 2e-4 ppm): 1000 ppm で 1。範囲の中で強く曲がるので、片方の端が動かない regula falsi は遅い
fn curved(x: f64) -> f64 {
    1.0 / (0.8 + 2.0e-4 * x)
}

#[test]
fn linear_values_reach_the_root_one_solve_after_the_ends() {
    let found = find(
        [800.0, 1300.0],
        1.0,
        1.0e-12,
        10,
        values(|x| 1.1 - 1.0e-4 * x),
    )
    .unwrap();
    assert_eq!(found.solves, 3);
    assert!((found.x - 1000.0).abs() < 1.0e-8, "{}", found.x);
    assert!((found.value - 1.0).abs() <= 1.0e-12);
    // 最後に評価した点のものを返す
    assert_eq!(found.at, found.x);
}

#[test]
fn curved_values_converge_within_the_tolerance() {
    let found = find([0.0, 3000.0], 1.0, 1.0e-10, 30, values(curved)).unwrap();
    assert!((found.value - 1.0).abs() <= 1.0e-10, "{}", found.value);
    assert!((found.x - 1000.0).abs() < 1.0e-6, "{}", found.x);
    assert_eq!(found.value, curved(found.x));
    // Illinois 法は超 1 次で収束する（この関数と範囲で、片方の端が残る regula falsi は 16 回、二分法は 35 回）
    assert!(found.solves <= 12, "{} solves", found.solves);
}

#[test]
fn decreasing_and_increasing_values_both_converge() {
    let down = find([0.0, 3000.0], 1.0, 1.0e-10, 30, values(curved)).unwrap();
    let up = find([0.0, 3000.0], -1.0, 1.0e-10, 30, values(|x| -curved(x))).unwrap();
    assert!((down.x - 1000.0).abs() < 1.0e-6);
    assert!((up.x - 1000.0).abs() < 1.0e-6);
}

#[test]
fn an_end_within_the_tolerance_is_the_root() {
    let found = find([1000.0, 1300.0], 1.0, 1.0e-9, 10, values(curved)).unwrap();
    assert_eq!((found.x, found.solves), (1000.0, 1));
    let found = find([500.0, 1000.0], 1.0, 1.0e-9, 10, values(curved)).unwrap();
    assert_eq!((found.x, found.solves), (1000.0, 2));
}

#[test]
fn a_bracket_without_the_target_is_an_error() {
    let error = find([1100.0, 1300.0], 1.0, 1.0e-9, 10, values(curved)).unwrap_err();
    assert!(
        matches!(error, SearchError::NotBracketed { low, high, .. } if low == 1100.0 && high == 1300.0),
        "{error}"
    );
}

#[test]
fn the_last_bracket_reports_the_evaluated_values() {
    // Illinois 法が残った端の値を半分にしたあとでも、エラーはその端で評価した値を示す
    let error = find([0.0, 3000.0], 1.0, 1.0e-10, 6, values(curved)).unwrap_err();
    let SearchError::NotConverged {
        low,
        high,
        at_low,
        at_high,
        ..
    } = error
    else {
        panic!("{error}");
    };
    assert_eq!(at_low, curved(low));
    assert_eq!(at_high, curved(high));
}

#[test]
fn exceeding_the_solves_is_an_error() {
    let error = find([0.0, 3000.0], 1.0, 1.0e-10, 4, values(curved)).unwrap_err();
    assert!(
        matches!(error, SearchError::NotConverged { solves: 4, .. }),
        "{error}"
    );
}

#[test]
fn errors_of_the_evaluation_pass_through() {
    let error = find([0.0, 3000.0], 1.0, 1.0e-10, 10, |x| {
        if x > 0.0 {
            Err("outside the table".to_string())
        } else {
            Ok((curved(x), ()))
        }
    })
    .unwrap_err();
    assert!(
        matches!(error, SearchError::Solve(ref e) if e == "outside the table"),
        "{error}"
    );
}
