use super::*;

fn table(law: Interpolation, x: &[f64], y: &[f64]) -> Tab1 {
    Tab1 {
        breaks: vec![x.len()],
        laws: vec![law],
        x: x.to_vec(),
        y: y.to_vec(),
    }
}

fn assert_close(actual: Option<f64>, expected: f64) {
    let actual = actual.expect("value inside the table range");
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs(),
        "actual {actual}, expected {expected}"
    );
}

#[test]
fn non_finite_values_are_rejected() {
    let one = |x: f64, y: f64| Tab1::linear(vec![x], vec![y]);
    assert!(one(f64::NAN, 1.0).check().is_err());
    assert!(one(1.0, f64::INFINITY).check().is_err());
    assert!(
        Tab1::linear(vec![1.0, f64::INFINITY], vec![1.0, 2.0])
            .check()
            .is_err()
    );
    assert!(one(1.0, 2.0).check().is_ok());
    assert!(ClampedTab1::new(one(f64::NAN, 1.0)).is_none());
}

#[test]
fn inconsistent_tables_are_rejected() {
    let mut short_y = Tab1::linear(vec![1.0, 2.0], vec![1.0, 1.0]);
    short_y.y.pop();
    assert!(short_y.check().is_err());
    let mut wrong_break = Tab1::linear(vec![1.0, 2.0], vec![1.0, 1.0]);
    wrong_break.breaks = vec![3];
    assert!(wrong_break.check().is_err());
    let mut extra_law = Tab1::linear(vec![1.0, 2.0], vec![1.0, 1.0]);
    extra_law.laws.push(Interpolation::LinLin);
    assert!(extra_law.check().is_err());
    // 範囲の端に丸める表は、検査に通る表からしか作れない
    let descending = Tab1::linear(vec![3.0, 1.0], vec![1.0, 1.0]);
    assert!(ClampedTab1::new(descending).is_none());
    assert!(ClampedTab1::new(wrong_break).is_none());
}

#[test]
fn clamped_table_extends_its_end_values() {
    assert!(ClampedTab1::new(Tab1::linear(vec![], vec![])).is_none());
    let t = ClampedTab1::new(Tab1::linear(vec![1.0, 3.0], vec![10.0, 30.0])).unwrap();
    assert_eq!(t.eval(0.5), 10.0);
    assert_eq!(t.eval(2.0), 20.0);
    assert_eq!(t.eval(9.0), 30.0);
}

#[test]
fn interpolation_codes() {
    assert_eq!(Interpolation::from_code(1), Some(Interpolation::Histogram));
    assert_eq!(Interpolation::from_code(5), Some(Interpolation::LogLog));
    assert_eq!(Interpolation::from_code(0), None);
    assert_eq!(Interpolation::from_code(6), None);
}

#[test]
fn histogram_takes_left_value() {
    let t = table(
        Interpolation::Histogram,
        &[1.0, 2.0, 3.0],
        &[10.0, 20.0, 30.0],
    );
    assert_close(t.eval(1.5), 10.0);
    assert_close(t.eval(2.0), 20.0);
}

#[test]
fn lin_lin() {
    let t = table(Interpolation::LinLin, &[1.0, 3.0], &[10.0, 30.0]);
    assert_close(t.eval(2.0), 20.0);
}

#[test]
fn lin_log() {
    let t = table(Interpolation::LinLog, &[1.0, 100.0], &[0.0, 2.0]);
    assert_close(t.eval(10.0), 1.0);
}

#[test]
fn log_lin() {
    let t = table(Interpolation::LogLin, &[0.0, 2.0], &[1.0, 100.0]);
    assert_close(t.eval(1.0), 10.0);
}

#[test]
fn log_log() {
    // 1/v 則: σ ∝ E^(-1/2)
    let t = table(Interpolation::LogLog, &[1.0, 100.0], &[10.0, 1.0]);
    assert_close(t.eval(10.0), 10.0_f64.sqrt());
}

#[test]
fn endpoints_are_inside_and_outside_is_none() {
    let t = table(Interpolation::LinLin, &[1.0, 3.0], &[10.0, 30.0]);
    assert_close(t.eval(1.0), 10.0);
    assert_close(t.eval(3.0), 30.0);
    assert_eq!(t.eval(0.999), None);
    assert_eq!(t.eval(3.001), None);
}

#[test]
fn law_changes_at_break() {
    let t = Tab1 {
        breaks: vec![2, 3],
        laws: vec![Interpolation::Histogram, Interpolation::LinLin],
        x: vec![1.0, 2.0, 4.0],
        y: vec![5.0, 10.0, 30.0],
    };
    assert_close(t.eval(1.5), 5.0);
    assert_close(t.eval(3.0), 20.0);
}

#[test]
fn discontinuity_returns_right_limit() {
    let t = table(
        Interpolation::LinLin,
        &[1.0, 2.0, 2.0, 3.0],
        &[1.0, 2.0, 7.0, 8.0],
    );
    assert_close(t.eval(2.0), 7.0);
    assert_close(t.eval(1.5), 1.5);
    assert_close(t.eval(2.5), 7.5);
}

#[test]
fn log_y_with_zero_endpoint_is_zero_inside() {
    // ln y → -∞ の極限として、区間内は 0
    for law in [Interpolation::LogLin, Interpolation::LogLog] {
        let t = table(law, &[1.0, 2.0, 3.0], &[0.0, 1.0, 0.0]);
        assert_eq!(t.eval(1.5), Some(0.0));
        assert_eq!(t.eval(2.5), Some(0.0));
        assert_close(t.eval(2.0), 1.0);
    }
}
