//! TAB1 レコードが表す1次元の表と、その補間。

/// ENDF の補間則（INT）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Interpolation {
    /// INT=1: y は区間内で左端の値をとる
    Histogram,
    /// INT=2: y は x に線形
    LinLin,
    /// INT=3: y は ln x に線形
    LinLog,
    /// INT=4: ln y は x に線形
    LogLin,
    /// INT=5: ln y は ln x に線形
    LogLog,
}

impl Interpolation {
    /// INT の値から変換する。荷電粒子用の INT=6 などは扱わない。
    pub fn from_code(code: i64) -> Option<Self> {
        todo!()
    }
}

/// 区間ごとに補間則を持つ (x, y) の表。
///
/// `breaks[k]` は補間則 `laws[k]` を適用する最後の点の番号（1始まり、ENDF の NBT）。
#[derive(Debug, Clone, PartialEq)]
pub struct Tab1 {
    pub breaks: Vec<usize>,
    pub laws: Vec<Interpolation>,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

impl Tab1 {
    /// x における値を補間する。表の範囲外では `None` を返す。
    ///
    /// x が重複する点（不連続点）上では、右側の極限値を返す。
    pub fn eval(&self, x: f64) -> Option<f64> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
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
}
