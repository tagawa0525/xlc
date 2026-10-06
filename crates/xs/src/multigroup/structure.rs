//! 群構造。

/// 群構造。境界はエネルギー (eV) の降順で、群 0 が最も高エネルギー。
#[derive(Debug, Clone, PartialEq)]
pub struct GroupStructure {
    bounds: Vec<f64>,
}

/// WIMS 69群の境界（NJOY groupr.f90 の eg9、エネルギーの昇順）
const WIMS69_ASCENDING: [f64; 70] = [
    1.0e-5, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.042, 0.05, 0.058, 0.067, 0.08, 0.1,
    0.14, 0.18, 0.22, 0.25, 0.28, 0.3, 0.32, 0.35, 0.4, 0.5, 0.625, 0.78, 0.85, 0.91, 0.95, 0.972,
    0.996, 1.02, 1.045, 1.071, 1.097, 1.123, 1.15, 1.3, 1.5, 2.1, 2.6, 3.3, 4.0, 9.877, 15.968,
    27.7, 48.052, 75.501, 148.728, 367.262, 906.898, 1425.1, 2239.45, 3519.1, 5530.0, 9118.0,
    1.503e4, 2.478e4, 4.085e4, 6.734e4, 1.11e5, 1.83e5, 3.025e5, 5.0e5, 8.21e5, 1.353e6, 2.231e6,
    3.679e6, 6.0655e6, 1.0e7,
];

/// 群構造の境界の誤り。
#[derive(Debug, thiserror::Error)]
#[error(
    "group bounds must be at least two positive, finite, strictly descending values, \
     got {0:?}"
)]
pub struct GroupError(Vec<f64>);

impl GroupStructure {
    /// 降順の境界から作る。境界は2つ以上、正で有限、厳密に降順であること
    /// （群積分で区間を対数的に分割するため、0 eV は扱えない）。
    pub fn new(bounds: Vec<f64>) -> Result<Self, GroupError> {
        let valid = bounds.len() >= 2
            && bounds.iter().all(|&b| b > 0.0 && b.is_finite())
            && bounds.windows(2).all(|w| w[0] > w[1]);
        if valid {
            Ok(Self { bounds })
        } else {
            Err(GroupError(bounds))
        }
    }

    /// WIMS 69群（NJOY GROUPR の ign=9 と同じ境界）。定数の表なので検査せずに作る
    /// （検査に通ることはテストで確かめている）。
    pub fn wims69() -> Self {
        Self {
            bounds: WIMS69_ASCENDING.iter().rev().copied().collect(),
        }
    }

    pub fn len(&self) -> usize {
        self.bounds.len() - 1
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn bounds(&self) -> &[f64] {
        &self.bounds
    }

    /// 群 g の (下限, 上限) エネルギー。
    pub fn range(&self, g: usize) -> (f64, f64) {
        (self.bounds[g + 1], self.bounds[g])
    }

    /// エネルギー e を含む群。境界上の点は低い側の群に入れる。範囲外なら `None`。
    pub fn find(&self, e: f64) -> Option<usize> {
        let (lowest, highest) = (self.bounds[self.len()], self.bounds[0]);
        if !(lowest..=highest).contains(&e) {
            return None;
        }
        let above_or_equal = self.bounds.partition_point(|&b| b >= e);
        Some((above_or_equal - 1).min(self.len() - 1))
    }
}

#[cfg(test)]
mod tests {
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
}
