//! 群構造。

/// 群構造。境界はエネルギー (eV) の降順で、群 0 が最も高エネルギー。
#[derive(Debug, Clone, PartialEq)]
pub struct GroupStructure {
    bounds: Vec<f64>,
}

impl GroupStructure {
    /// 降順の境界から作る。
    ///
    /// # Panics
    /// 境界が2つ未満のとき、または降順でないとき。
    pub fn new(bounds: Vec<f64>) -> Self {
        todo!()
    }

    /// WIMS 69群（NJOY GROUPR の ign=9 と同じ境界）。
    pub fn wims69() -> Self {
        todo!()
    }

    pub fn len(&self) -> usize {
        todo!()
    }

    pub fn is_empty(&self) -> bool {
        todo!()
    }

    pub fn bounds(&self) -> &[f64] {
        todo!()
    }

    /// 群 g の (下限, 上限) エネルギー。
    pub fn range(&self, g: usize) -> (f64, f64) {
        todo!()
    }

    /// エネルギー e を含む群。境界上の点は低い側の群に入れる。範囲外なら `None`。
    pub fn find(&self, e: f64) -> Option<usize> {
        todo!()
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
        let g = GroupStructure::new(vec![10.0, 5.0, 1.0]);
        assert_eq!(g.find(7.0), Some(0));
        assert_eq!(g.find(5.0), Some(1));
        assert_eq!(g.find(10.0), Some(0));
        assert_eq!(g.find(1.0), Some(1));
        assert_eq!(g.find(0.5), None);
        assert_eq!(g.find(11.0), None);
    }

    #[test]
    #[should_panic]
    fn ascending_bounds_panic() {
        GroupStructure::new(vec![1.0, 5.0]);
    }
}
