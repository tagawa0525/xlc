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

    /// x と y の長さ、区切りと補間則の対応（区切りは狭義の昇順で、最後が点の数）、値がすべて有限で
    /// x が昇順であること、各区間の値が補間則の定義域にあることを検査する。
    pub(super) fn check(&self) -> Result<(), String> {
        if self.x.len() != self.y.len() {
            return Err(format!(
                "{} x values but {} y values",
                self.x.len(),
                self.y.len()
            ));
        }
        if self.breaks.len() != self.laws.len() {
            return Err(format!(
                "{} interpolation ranges but {} laws",
                self.breaks.len(),
                self.laws.len()
            ));
        }
        let ascending = self.breaks.windows(2).all(|w| w[0] < w[1]);
        if !ascending || self.breaks.last() != Some(&self.x.len()) {
            return Err(format!(
                "interpolation ranges {:?} must ascend to the {} points",
                self.breaks,
                self.x.len()
            ));
        }
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
    /// 空の表、または検査（`Tab1::check`）に通らない表なら `None`。範囲の端に丸める `clamp` と
    /// 補間則の参照が、昇順で有限な x と整合した区切りを前提にするため。
    pub fn new(table: Tab1) -> Option<Self> {
        let valid = !table.x.is_empty() && table.check().is_ok();
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
mod tests;
