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
        match code {
            1 => Some(Self::Histogram),
            2 => Some(Self::LinLin),
            3 => Some(Self::LinLog),
            4 => Some(Self::LogLin),
            5 => Some(Self::LogLog),
            _ => None,
        }
    }

    fn log_x(self) -> bool {
        matches!(self, Self::LinLog | Self::LogLog)
    }

    fn log_y(self) -> bool {
        matches!(self, Self::LogLin | Self::LogLog)
    }

    pub fn interpolate(self, x: f64, (x0, y0): (f64, f64), (x1, y1): (f64, f64)) -> f64 {
        if x == x0 {
            return y0;
        }
        // 端点の y が 0 なら ln y → -∞ となり、区間内の値は極限として 0
        if self.log_y() && (y0 == 0.0 || y1 == 0.0) {
            return 0.0;
        }
        match self {
            Self::Histogram => y0,
            Self::LinLin => y0 + (y1 - y0) * (x - x0) / (x1 - x0),
            Self::LinLog => y0 + (y1 - y0) * (x / x0).ln() / (x1 / x0).ln(),
            Self::LogLin => y0 * (y1 / y0).powf((x - x0) / (x1 - x0)),
            Self::LogLog => y0 * (y1 / y0).powf((x / x0).ln() / (x1 / x0).ln()),
        }
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
    /// 全区間を lin-lin で補間する表。
    pub fn linear(x: Vec<f64>, y: Vec<f64>) -> Self {
        Self {
            breaks: vec![x.len()],
            laws: vec![Interpolation::LinLin],
            x,
            y,
        }
    }

    /// 区間 (k, k+1)（0始まり）に適用する補間則。
    pub fn law_of_interval(&self, k: usize) -> Interpolation {
        self.laws[self.breaks.partition_point(|&b| b < k + 2)]
    }

    /// 値がすべて有限で、x が昇順であることと、各区間の値が補間則の定義域にあることを検査する。
    /// `breaks` と `laws` の整合は呼び出し側で保証されているものとする。
    pub(super) fn check(&self) -> Result<(), String> {
        if let Some(k) = self
            .x
            .iter()
            .zip(&self.y)
            .position(|(x, y)| !x.is_finite() || !y.is_finite())
        {
            return Err(format!(
                "non-finite value at point {}: x = {}, y = {}",
                k + 1,
                self.x[k],
                self.y[k]
            ));
        }
        for k in 0..self.x.len().saturating_sub(1) {
            let (x0, x1) = (self.x[k], self.x[k + 1]);
            if x0.partial_cmp(&x1).is_none_or(|o| o.is_gt()) {
                return Err(format!(
                    "x values are not in ascending order at point {}",
                    k + 1
                ));
            }
            let law = self.law_of_interval(k);
            if law.log_x() && x0 <= 0.0 {
                return Err(format!(
                    "non-positive x {x0} with logarithmic x interpolation"
                ));
            }
            let (y0, y1) = (self.y[k], self.y[k + 1]);
            if law.log_y() && (y0 < 0.0 || y1 < 0.0) {
                return Err(format!(
                    "negative y with logarithmic y interpolation at point {}",
                    k + 1
                ));
            }
        }
        Ok(())
    }

    /// x における値を補間する。表の範囲外では `None` を返す。
    ///
    /// x が重複する点（不連続点）上では、右側の極限値を返す。
    pub fn eval(&self, x: f64) -> Option<f64> {
        let (&first, &last) = (self.x.first()?, self.x.last()?);
        if !(first..=last).contains(&x) {
            return None;
        }
        Some(self.value_within(x))
    }

    /// 範囲の中（first ≤ x ≤ last）の x で補間する。
    fn value_within(&self, x: f64) -> f64 {
        // x[k] <= x を満たす最後の点。重複点では後ろの点が選ばれ、右側の極限になる
        let k = self.x.partition_point(|&xi| xi <= x).saturating_sub(1);
        if k == self.x.len() - 1 {
            return self.y[k];
        }
        // 区間 (k, k+1) の右端の点番号（1始まり）を含む範囲の補間則を使う
        let law = self.law_of_interval(k);
        law.interpolate(x, (self.x[k], self.y[k]), (self.x[k + 1], self.y[k + 1]))
    }
}

/// 範囲の外を端の値で延長して補間する、空でない表（MF5 の θ(E) など）。
#[derive(Debug, Clone, PartialEq)]
pub struct ClampedTab1 {
    table: Tab1,
}

impl ClampedTab1 {
    /// 空の表、または有限でない x がある表なら `None`（範囲の端に丸める clamp が NaN で panic するため）。
    pub fn new(table: Tab1) -> Option<Self> {
        let valid = !table.x.is_empty() && table.x.iter().all(|x| x.is_finite());
        valid.then_some(Self { table })
    }

    pub fn table(&self) -> &Tab1 {
        &self.table
    }

    /// x における値。表の範囲の外では端の値。
    pub fn eval(&self, x: f64) -> f64 {
        let t = &self.table;
        // 空でないことは作るときに確かめている
        let (first, last) = (t.x[0], t.x[t.x.len() - 1]);
        t.value_within(x.clamp(first, last))
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
}
