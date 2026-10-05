//! ピンセルの形状。

use std::f64::consts::PI;

/// 中心が原点の正方形のセル（一辺 `pitch` cm）と、中心を共有する同心円。
///
/// 領域 i（0 始まり）は半径 radii[i-1] と radii[i] の間の円環で、最後の領域は
/// 最も外側の円より外のセルの残り。
#[derive(Debug, Clone, PartialEq)]
pub struct PinCell {
    pitch: f64,
    radii: Vec<f64>,
}

impl PinCell {
    /// 半径は昇順で、セルの内接円より小さいこと。
    pub fn new(pitch: f64, radii: Vec<f64>) -> Result<Self, String> {
        todo!()
    }

    pub fn pitch(&self) -> f64 {
        self.pitch
    }

    pub fn radii(&self) -> &[f64] {
        &self.radii
    }

    pub fn regions(&self) -> usize {
        self.radii.len() + 1
    }

    /// 点 (x, y) を含む領域。
    pub fn region_at(&self, x: f64, y: f64) -> usize {
        todo!()
    }

    /// 領域ごとの体積（単位高さあたりの面積, cm²）の解析値。
    pub fn volumes(&self) -> Vec<f64> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validates_radii() {
        assert!(PinCell::new(1.26, vec![0.2, 0.4]).is_ok());
        assert!(PinCell::new(1.26, vec![0.4, 0.2]).is_err());
        assert!(PinCell::new(1.26, vec![0.63]).is_err());
        assert!(PinCell::new(1.26, vec![0.0]).is_err());
        assert!(PinCell::new(0.0, vec![]).is_err());
    }

    #[test]
    fn regions_and_volumes() {
        let cell = PinCell::new(2.0, vec![0.5, 0.8]).unwrap();
        assert_eq!(cell.regions(), 3);
        assert_eq!(cell.region_at(0.1, 0.1), 0);
        assert_eq!(cell.region_at(0.6, 0.0), 1);
        assert_eq!(cell.region_at(0.9, 0.9), 2);
        let v = cell.volumes();
        assert!((v[0] - PI * 0.25).abs() < 1e-15);
        assert!((v[1] - PI * (0.64 - 0.25)).abs() < 1e-15);
        assert!((v.iter().sum::<f64>() - 4.0).abs() < 1e-14);
    }
}
