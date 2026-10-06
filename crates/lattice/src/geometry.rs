//! 格子計算の形状: 同心円で区切った正方形のセル（ピンセル）と、それを正方格子に並べたもの。

use std::f64::consts::PI;

/// 中心が原点の正方形のセル（一辺 `pitch` cm）と、中心を共有する同心円。
///
/// 領域 i（0 始まり）は半径 radii[i-1] と radii[i] の間の円環で、最後の領域は
/// 最も外側の円より外のセルの残り。
///
/// 同心円のないセルは、m × m の小さな正方形に分けられる（`grid`）。領域は上の行から、各行は左から
/// 番号を付ける。真空境界に接する反射体のように中性子束が大きく変わる均質な領域を、
/// 平坦線源の近似で表せる大きさに分けるため。
#[derive(Debug, Clone, PartialEq)]
pub struct PinCell {
    pitch: f64,
    radii: Vec<f64>,
    /// 1辺の小さな正方形の数（同心円のあるセルでは 1）
    grid: usize,
}

impl PinCell {
    /// 半径は昇順で、セルの内接円より小さいこと。
    pub fn new(pitch: f64, radii: Vec<f64>) -> Result<Self, String> {
        if !(pitch > 0.0 && pitch.is_finite()) {
            return Err(format!("pitch must be positive, got {pitch}"));
        }
        let ascending = radii.windows(2).all(|w| w[0] < w[1]);
        let inside = radii.iter().all(|&r| r > 0.0 && r < 0.5 * pitch);
        if !(ascending && inside) {
            return Err(format!(
                "radii must be ascending and between 0 and pitch/2, got {radii:?}"
            ));
        }
        Ok(Self {
            pitch,
            radii,
            grid: 1,
        })
    }

    pub fn pitch(&self) -> f64 {
        self.pitch
    }

    pub fn radii(&self) -> &[f64] {
        &self.radii
    }

    pub fn regions(&self) -> usize {
        self.radii.len() + self.grid * self.grid
    }

    /// 同心円がなく、m × m の小さな正方形に分けたセル。
    pub fn grid(pitch: f64, m: usize) -> Result<Self, String> {
        if m == 0 {
            return Err("a grid needs at least one square per side".into());
        }
        Ok(Self {
            grid: m,
            ..Self::new(pitch, Vec::new())?
        })
    }

    /// 点 (x, y) を含む領域。
    pub fn region_at(&self, x: f64, y: f64) -> usize {
        if self.grid > 1 {
            let (m, h) = (self.grid, self.pitch / self.grid as f64);
            let index = |v: f64| ((v / h).floor().max(0.0) as usize).min(m - 1);
            return index(0.5 * self.pitch - y) * m + index(x + 0.5 * self.pitch);
        }
        let r = x.hypot(y);
        self.radii.partition_point(|&ri| ri <= r)
    }

    /// セルの中心を原点とする点 start から方向 direction に長さ length 進む直線を、
    /// 同心円で区間に分ける。
    pub fn segments(&self, start: (f64, f64), direction: (f64, f64), length: f64) -> Vec<Segment> {
        // |start + t d|² = r² の解のうち (0, length) にあるもの
        let b = start.0 * direction.0 + start.1 * direction.1;
        let c0 = start.0 * start.0 + start.1 * start.1;
        let mut cuts = vec![0.0, length];
        // 小さな正方形の境界の線 x, y = -p/2 + k p/m
        let h = self.pitch / self.grid as f64;
        for (origin, d) in [(start.0, direction.0), (start.1, direction.1)] {
            if d != 0.0 {
                for k in 1..self.grid {
                    let t = (-0.5 * self.pitch + k as f64 * h - origin) / d;
                    if t > 0.0 && t < length {
                        cuts.push(t);
                    }
                }
            }
        }
        for &r in &self.radii {
            let disc = b * b - (c0 - r * r);
            if disc > 0.0 {
                let root = disc.sqrt();
                cuts.extend(
                    [-b - root, -b + root]
                        .into_iter()
                        .filter(|&t| t > 0.0 && t < length),
                );
            }
        }
        cuts.sort_by(f64::total_cmp);
        cuts.windows(2)
            .filter(|w| w[1] > w[0])
            .map(|w| {
                let mid = 0.5 * (w[0] + w[1]);
                Segment {
                    region: self
                        .region_at(start.0 + mid * direction.0, start.1 + mid * direction.1),
                    length: w[1] - w[0],
                }
            })
            .collect()
    }

    /// 領域ごとの体積（単位高さあたりの面積, cm²）の解析値。
    pub fn volumes(&self) -> Vec<f64> {
        if self.grid > 1 {
            let h = self.pitch / self.grid as f64;
            return vec![h * h; self.grid * self.grid];
        }
        let mut inner = 0.0;
        let mut volumes: Vec<f64> = self
            .radii
            .iter()
            .map(|&r| {
                let v = PI * (r * r - inner * inner);
                inner = r;
                v
            })
            .collect();
        volumes.push(self.pitch * self.pitch - PI * inner * inner);
        volumes
    }
}

/// トラック上の1つの領域を通る区間。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Segment {
    pub region: usize,
    /// 平面に射影した長さ (cm)
    pub length: f64,
}

/// ピンセルを n × n の正方格子に並べた、一辺 n × pitch の正方形の形状（中心が原点）。
///
/// セルは上の行から順に、各行は左から順に番号を付ける。領域はセルの順に、セルの中では
/// `PinCell` の領域の順に通し番号を付ける。
#[derive(Debug, Clone, PartialEq)]
pub struct Lattice {
    types: Vec<PinCell>,
    /// [セル] のセルの種類（`types` の番号）
    layout: Vec<usize>,
    n: usize,
    /// [セル] のセルの最初の領域の番号
    offsets: Vec<usize>,
}

impl Lattice {
    /// `types` はセルの種類（すべて同じピッチ）、`layout` は上の行から順の各行のセルの種類。
    pub fn new(types: Vec<PinCell>, layout: Vec<Vec<usize>>) -> Result<Self, String> {
        let n = layout.len();
        if n == 0 {
            return Err("the lattice has no cells".into());
        }
        if let Some(row) = layout.iter().position(|row| row.len() != n) {
            return Err(format!(
                "row {} has {} cells, but the lattice has {n} rows",
                row + 1,
                layout[row].len()
            ));
        }
        let pitch = types.first().ok_or("no cell types")?.pitch();
        if types.iter().any(|t| t.pitch() != pitch) {
            return Err("all cell types must have the same pitch".into());
        }
        let layout: Vec<usize> = layout.into_iter().flatten().collect();
        if let Some(&bad) = layout.iter().find(|&&t| t >= types.len()) {
            return Err(format!(
                "cell type {bad} is out of range ({} types)",
                types.len()
            ));
        }
        let mut offsets = Vec::with_capacity(layout.len() + 1);
        let mut next = 0;
        for &t in &layout {
            offsets.push(next);
            next += types[t].regions();
        }
        offsets.push(next);
        Ok(Self {
            types,
            layout,
            n,
            offsets,
        })
    }

    /// 1つのピンセルだけの形状。
    pub fn pin_cell(cell: PinCell) -> Self {
        Self::new(vec![cell], vec![vec![0]]).expect("a single cell is a valid lattice")
    }

    fn pitch(&self) -> f64 {
        self.types[0].pitch()
    }

    /// 1辺のセルの数
    pub fn cells_per_side(&self) -> usize {
        self.n
    }

    /// 一辺の長さ (cm)
    pub fn width(&self) -> f64 {
        self.n as f64 * self.pitch()
    }

    pub fn regions(&self) -> usize {
        self.offsets[self.layout.len()]
    }

    /// セル `cell` の領域の範囲。
    pub fn cell_regions(&self, cell: usize) -> std::ops::Range<usize> {
        self.offsets[cell]..self.offsets[cell + 1]
    }

    /// 領域ごとの体積（単位高さあたりの面積, cm²）の解析値。
    pub fn volumes(&self) -> Vec<f64> {
        self.layout
            .iter()
            .flat_map(|&t| self.types[t].volumes())
            .collect()
    }

    /// 点 start から方向 direction（単位ベクトル）に長さ length 進む直線を、領域ごとの区間に分ける。
    pub fn segments(&self, start: (f64, f64), direction: (f64, f64), length: f64) -> Vec<Segment> {
        let (pitch, half) = (self.pitch(), 0.5 * self.width());
        // 格子線 x, y = -W/2 + k p と交わる点で、セルごとの区間に分ける
        let mut cuts = vec![0.0, length];
        for (origin, d) in [(start.0, direction.0), (start.1, direction.1)] {
            if d != 0.0 {
                for k in 1..self.n {
                    let t = (-half + k as f64 * pitch - origin) / d;
                    if t > 0.0 && t < length {
                        cuts.push(t);
                    }
                }
            }
        }
        cuts.sort_by(f64::total_cmp);
        let mut segments = Vec::new();
        for w in cuts.windows(2).filter(|w| w[1] > w[0]) {
            let mid = 0.5 * (w[0] + w[1]);
            let (x, y) = (start.0 + mid * direction.0, start.1 + mid * direction.1);
            let column = (((x + half) / pitch).floor() as usize).min(self.n - 1);
            let row = (((half - y) / pitch).floor() as usize).min(self.n - 1);
            let cell = row * self.n + column;
            let center = (
                -half + (column as f64 + 0.5) * pitch,
                half - (row as f64 + 0.5) * pitch,
            );
            let local = (
                start.0 + w[0] * direction.0 - center.0,
                start.1 + w[0] * direction.1 - center.1,
            );
            let offset = self.offsets[cell];
            segments.extend(
                self.types[self.layout[cell]]
                    .segments(local, direction, w[1] - w[0])
                    .into_iter()
                    .map(|s| Segment {
                        region: offset + s.region,
                        length: s.length,
                    }),
            );
        }
        segments
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

    fn two_by_two() -> Lattice {
        let pin = PinCell::new(2.0, vec![0.5]).unwrap();
        let water = PinCell::new(2.0, vec![]).unwrap();
        Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 0]]).unwrap()
    }

    #[test]
    fn validates_lattice() {
        let pin = PinCell::new(2.0, vec![0.5]).unwrap();
        let other = PinCell::new(1.0, vec![]).unwrap();
        // ピッチが揃っていない、行の長さが行の数と違う、種類の番号が範囲外
        assert!(Lattice::new(vec![pin.clone(), other], vec![vec![0, 1], vec![1, 0]]).is_err());
        assert!(Lattice::new(vec![pin.clone()], vec![vec![0, 0]]).is_err());
        assert!(Lattice::new(vec![pin.clone()], vec![vec![1]]).is_err());
        assert!(Lattice::new(vec![pin], vec![]).is_err());
    }

    #[test]
    fn lattice_regions_and_volumes() {
        let lattice = two_by_two();
        assert_eq!(lattice.width(), 4.0);
        // 上の行: ピン（2領域）、水（1領域）。下の行: 水、ピン
        assert_eq!(lattice.regions(), 6);
        assert_eq!(lattice.cell_regions(0), 0..2);
        assert_eq!(lattice.cell_regions(1), 2..3);
        assert_eq!(lattice.cell_regions(2), 3..4);
        assert_eq!(lattice.cell_regions(3), 4..6);
        let v = lattice.volumes();
        let disk = PI * 0.25;
        let expected = [disk, 4.0 - disk, 4.0, 4.0, disk, 4.0 - disk];
        for (a, b) in v.iter().zip(expected) {
            assert!((a - b).abs() < 1e-14, "{v:?}");
        }
    }

    #[test]
    fn segments_cross_cells_in_order() {
        // 上の行の中心を左から右へ: ピンの外側, ピンの内側, ピンの外側, 水
        let lattice = two_by_two();
        let s = lattice.segments((-2.0, 1.0), (1.0, 0.0), 4.0);
        let regions: Vec<usize> = s.iter().map(|s| s.region).collect();
        assert_eq!(regions, vec![1, 0, 1, 2]);
        let lengths: Vec<f64> = s.iter().map(|s| s.length).collect();
        for (a, b) in lengths.iter().zip([0.5, 1.0, 0.5, 2.0]) {
            assert!((a - b).abs() < 1e-14, "{lengths:?}");
        }
    }

    #[test]
    fn diagonal_segments_cover_the_line() {
        let lattice = two_by_two();
        let d = std::f64::consts::FRAC_1_SQRT_2;
        let s = lattice.segments((-2.0, 2.0), (d, -d), 4.0 * 2.0_f64.sqrt());
        let length: f64 = s.iter().map(|s| s.length).sum();
        assert!((length - 4.0 * 2.0_f64.sqrt()).abs() < 1e-12);
        // 左上と右下のピンを中心を通って横切る（直径 1.0 の区間がそれぞれある）
        let inner: Vec<f64> = s
            .iter()
            .filter(|s| s.region == 0 || s.region == 4)
            .map(|s| s.length)
            .collect();
        assert_eq!(inner.len(), 2);
        assert!(inner.iter().all(|l| (l - 1.0).abs() < 1e-12), "{inner:?}");
    }

    #[test]
    fn single_pin_cell_lattice() {
        let cell = PinCell::new(1.26, vec![0.4]).unwrap();
        let lattice = Lattice::pin_cell(cell.clone());
        assert_eq!(lattice.width(), 1.26);
        assert_eq!(lattice.regions(), 2);
        assert_eq!(lattice.volumes(), cell.volumes());
    }

    #[test]
    fn grid_cell_regions_and_segments() {
        assert!(PinCell::grid(2.0, 0).is_err());
        let cell = PinCell::grid(2.0, 2).unwrap();
        assert_eq!(cell.regions(), 4);
        assert_eq!(cell.volumes(), vec![1.0; 4]);
        // 左上、右上、左下、右下の順
        assert_eq!(cell.region_at(-0.5, 0.5), 0);
        assert_eq!(cell.region_at(0.5, 0.5), 1);
        assert_eq!(cell.region_at(-0.5, -0.5), 2);
        assert_eq!(cell.region_at(0.5, -0.5), 3);
        // 下の行を左から右へ
        let s = cell.segments((-1.0, -0.5), (1.0, 0.0), 2.0);
        assert_eq!(
            s,
            vec![
                Segment {
                    region: 2,
                    length: 1.0
                },
                Segment {
                    region: 3,
                    length: 1.0
                }
            ]
        );
        // 1 × 1 なら同心円のないセルと同じ
        assert_eq!(
            PinCell::grid(2.0, 1).unwrap(),
            PinCell::new(2.0, vec![]).unwrap()
        );
    }
}
