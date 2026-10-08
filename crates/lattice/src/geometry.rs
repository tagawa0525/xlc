//! 格子計算の形状: 同心円で区切った正方形のセル（ピンセル）と、それを正方格子に並べたもの。

use std::f64::consts::PI;

/// 中心が原点の正方形のセル（一辺 `pitch` cm）と、中心を共有する同心円。
///
/// 領域 i（0 始まり）は半径 radii[i-1] と radii[i] の間の円環で、最後の領域は
/// 最も外側の円より外のセルの残り。内接円より大きい円はセルの外周で切り取る
/// （ピンセルの四隅の減速材を半径方向に分けるため）。
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
    /// 各リングを分ける扇形の数（1, 4, 8）
    sectors: usize,
}

/// 形状の誤り。
#[derive(Debug, thiserror::Error)]
pub enum GeometryError {
    #[error("pitch must be positive, got {0}")]
    Pitch(f64),
    #[error("radii must be ascending and between 0 and pitch/√2, got {0:?}")]
    Radii(Vec<f64>),
    #[error("the number of sectors must be 1, 4 or 8, got {0}")]
    Sectors(usize),
    #[error("a grid needs at least one square per side")]
    EmptyGrid,
    #[error("the lattice has no cells")]
    NoCells,
    #[error("row {row} has {cells} cells, but the lattice has {rows} rows")]
    RowLength {
        row: usize,
        cells: usize,
        rows: usize,
    },
    #[error("no cell types")]
    NoCellTypes,
    #[error("all cell types must have the same pitch")]
    MixedPitch,
    #[error("cell type {index} is out of range ({types} types)")]
    CellType { index: usize, types: usize },
    #[error("the margin around the lattice must be finite and not negative, got {0}")]
    Margin(f64),
    #[error("frame thickness must be positive, below pitch/2 and outside the circles, got {0}")]
    Frame(f64),
}

impl PinCell {
    /// 半径は昇順で、セルの外接円より小さいこと（外接円以上の円の外側は空になる）。
    pub fn new(pitch: f64, radii: Vec<f64>) -> Result<Self, GeometryError> {
        if !(pitch > 0.0 && pitch.is_finite()) {
            return Err(GeometryError::Pitch(pitch));
        }
        let ascending = radii.windows(2).all(|w| w[0] < w[1]);
        let circumscribed = 0.5 * pitch * std::f64::consts::SQRT_2;
        let inside = radii.iter().all(|&r| r > 0.0 && r < circumscribed);
        if !(ascending && inside) {
            return Err(GeometryError::Radii(radii));
        }
        Ok(Self {
            pitch,
            radii,
            grid: 1,
            sectors: 1,
        })
    }

    pub fn pitch(&self) -> f64 {
        self.pitch
    }

    pub fn radii(&self) -> &[f64] {
        &self.radii
    }

    /// 各リングを分ける扇形の数
    pub fn sectors(&self) -> usize {
        self.sectors
    }

    /// 同心円ではなく小さな正方形に分けたセルか
    pub fn is_grid(&self) -> bool {
        self.grid > 1
    }

    pub fn regions(&self) -> usize {
        if self.grid > 1 {
            self.grid * self.grid
        } else {
            (self.radii.len() + 1) * self.sectors
        }
    }

    /// 同心円で区切った各リングを、さらに x 軸から等角度の扇形 `sectors` 個（1, 4, 8）に分けたセル。
    /// 領域はリングの順に、リングの中では角度の順に番号を付ける。
    pub fn with_sectors(
        pitch: f64,
        radii: Vec<f64>,
        sectors: usize,
    ) -> Result<Self, GeometryError> {
        // 正方形の対称性から、1, 4, 8 個なら外側の領域も扇形ごとに同じ面積になる
        if ![1, 4, 8].contains(&sectors) {
            return Err(GeometryError::Sectors(sectors));
        }
        Ok(Self {
            sectors,
            ..Self::new(pitch, radii)?
        })
    }

    /// 同心円がなく、m × m の小さな正方形に分けたセル。
    pub fn grid(pitch: f64, m: usize) -> Result<Self, GeometryError> {
        if m == 0 {
            return Err(GeometryError::EmptyGrid);
        }
        Ok(Self {
            grid: m,
            ..Self::new(pitch, Vec::new())?
        })
    }

    /// セルの外周に厚さ `thickness` cm の正方形の枠（支持格子を体積を保って薄い枠にしたもの）を付けたセル。
    /// 枠は1つの領域で、領域の番号の最後に置く。同心円の領域は、セルの外周の代わりに枠の内側で切り取る。
    pub fn with_frame(self, thickness: f64) -> Result<Self, GeometryError> {
        let _ = thickness;
        Ok(self)
    }

    /// 点 (x, y) を含む領域。
    pub fn region_at(&self, x: f64, y: f64) -> usize {
        if self.grid > 1 {
            let (m, h) = (self.grid, self.pitch / self.grid as f64);
            let index = |v: f64| ((v / h).floor().max(0.0) as usize).min(m - 1);
            return index(0.5 * self.pitch - y) * m + index(x + 0.5 * self.pitch);
        }
        let ring = self.radii.partition_point(|&ri| ri <= x.hypot(y));
        let angle = y.atan2(x).rem_euclid(2.0 * PI);
        let sector = ((angle / (2.0 * PI) * self.sectors as f64) as usize).min(self.sectors - 1);
        ring * self.sectors + sector
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
        // 扇形の境界の半直線 t u（u = (cos θk, sin θk)）: u × (start + t d) = 0
        if self.sectors > 1 {
            for k in 0..self.sectors {
                let theta = 2.0 * PI * k as f64 / self.sectors as f64;
                let u = (theta.cos(), theta.sin());
                let cross = |v: (f64, f64)| u.0 * v.1 - u.1 * v.0;
                let denominator = cross(direction);
                if denominator != 0.0 {
                    let t = -cross(start) / denominator;
                    let point = (start.0 + t * direction.0, start.1 + t * direction.1);
                    // 中心を通る直線では交点が中心そのもの（u·point = 0）になるので、等号も含める
                    if t > 0.0 && t < length && u.0 * point.0 + u.1 * point.1 >= 0.0 {
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
        let disk = |r: f64| clipped_disk_area(r, 0.5 * self.pitch);
        let mut inner = 0.0;
        let mut volumes: Vec<f64> = self
            .radii
            .iter()
            .map(|&r| {
                let v = disk(r) - disk(inner);
                inner = r;
                v
            })
            .collect();
        volumes.push(self.pitch * self.pitch - disk(inner));
        volumes
            .iter()
            .flat_map(|v| std::iter::repeat_n(v / self.sectors as f64, self.sectors))
            .collect()
    }
}

/// 半径 r の円と、中心を共有する半辺 h の正方形の共通部分の面積（r < h√2）。
///
/// r > h では、8 つの合同な 45° の扇形それぞれで、x = h の辺の外を切り取る。扇形の中で円が辺と
/// 交わる角を θc（cos θc = h/r）とすると、扇形の面積は (h² tan θc + r² (π/4 - θc)) / 2。
fn clipped_disk_area(r: f64, h: f64) -> f64 {
    if r <= h {
        return PI * r * r;
    }
    let theta = (h / r).acos();
    4.0 * (h * h * theta.tan() + r * r * (0.25 * PI - theta))
}

/// トラック上の1つの領域を通る区間。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Segment {
    pub region: usize,
    /// 平面に射影した長さ (cm)
    pub length: f64,
}

/// ピンセルを n × n の正方格子に並べ、周りに幅 m の余白（集合体の間の水のギャップなど）を付けた、一辺
/// n × pitch + 2m の正方形の形状（中心が原点）。
///
/// セルは上の行から順に、各行は左から順に番号を付ける。領域はセルの順に、セルの中では
/// `PinCell` の領域の順に通し番号を付ける。余白があれば、その後に余白の領域を、上の帯（セルの列ごとに左から）、
/// 下の帯、左の帯（セルの行ごとに上から）、右の帯、角（左上、右上、左下、右下）の順に置く。帯をセルの幅で
/// 分けるのは、平坦な線源の近似で表せる大きさにするため。
#[derive(Debug, Clone, PartialEq)]
pub struct Lattice {
    types: Vec<PinCell>,
    /// [セル] のセルの種類（`types` の番号）
    layout: Vec<usize>,
    n: usize,
    /// [セル] のセルの最初の領域の番号
    offsets: Vec<usize>,
    /// 格子の周りの余白の幅 (cm)
    margin: f64,
}

impl Lattice {
    /// `types` はセルの種類（すべて同じピッチ）、`layout` は上の行から順の各行のセルの種類。
    pub fn new(types: Vec<PinCell>, layout: Vec<Vec<usize>>) -> Result<Self, GeometryError> {
        Self::with_margin(types, layout, 0.0)
    }

    /// 周りに幅 `margin` (cm) の余白を付けた格子（`new` を参照）。
    pub fn with_margin(
        types: Vec<PinCell>,
        layout: Vec<Vec<usize>>,
        margin: f64,
    ) -> Result<Self, GeometryError> {
        if !(margin >= 0.0 && margin.is_finite()) {
            return Err(GeometryError::Margin(margin));
        }
        let n = layout.len();
        if n == 0 {
            return Err(GeometryError::NoCells);
        }
        if let Some(row) = layout.iter().position(|row| row.len() != n) {
            return Err(GeometryError::RowLength {
                row: row + 1,
                cells: layout[row].len(),
                rows: n,
            });
        }
        let pitch = types.first().ok_or(GeometryError::NoCellTypes)?.pitch();
        if types.iter().any(|t| t.pitch() != pitch) {
            return Err(GeometryError::MixedPitch);
        }
        let layout: Vec<usize> = layout.into_iter().flatten().collect();
        if let Some(&bad) = layout.iter().find(|&&t| t >= types.len()) {
            return Err(GeometryError::CellType {
                index: bad,
                types: types.len(),
            });
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
            margin,
        })
    }

    /// 1つのピンセルだけの形状（`new` と同じものを、検査なしで直接作る）。
    pub fn pin_cell(cell: PinCell) -> Self {
        let regions = cell.regions();
        Self {
            types: vec![cell],
            layout: vec![0],
            n: 1,
            offsets: vec![0, regions],
            margin: 0.0,
        }
    }

    fn pitch(&self) -> f64 {
        self.types[0].pitch()
    }

    /// 1辺のセルの数
    pub fn cells_per_side(&self) -> usize {
        self.n
    }

    /// 一辺の長さ (cm)。余白を含む
    pub fn width(&self) -> f64 {
        self.n as f64 * self.pitch() + 2.0 * self.margin
    }

    /// 格子の周りの余白の幅 (cm)
    pub fn margin(&self) -> f64 {
        self.margin
    }

    pub fn regions(&self) -> usize {
        self.cell_region_count() + self.margin_regions()
    }

    /// セルの領域の数（余白の領域はこの後に続く）
    fn cell_region_count(&self) -> usize {
        self.offsets[self.layout.len()]
    }

    /// 余白の領域の数（余白がなければ 0）
    fn margin_regions(&self) -> usize {
        if self.margin > 0.0 { 4 * self.n + 4 } else { 0 }
    }

    /// 領域の属するセル。余白の領域は、接するセル（角は角のセル）
    pub fn region_cell(&self, region: usize) -> usize {
        let base = self.cell_region_count();
        if region < base {
            return self.offsets.partition_point(|&o| o <= region) - 1;
        }
        let (n, k) = (self.n, region - base);
        match k / n {
            // 上の帯、下の帯
            0 => k,
            1 => (n - 1) * n + (k - n),
            // 左の帯、右の帯
            2 => (k - 2 * n) * n,
            3 => (k - 3 * n) * n + n - 1,
            // 角
            _ => [0, n - 1, (n - 1) * n, n * n - 1][(k - 4 * n).min(3)],
        }
    }

    /// セル `cell` の種類
    pub fn cell_type(&self, cell: usize) -> &PinCell {
        &self.types[self.layout[cell]]
    }

    /// セル `cell` の領域の範囲。
    pub fn cell_regions(&self, cell: usize) -> std::ops::Range<usize> {
        self.offsets[cell]..self.offsets[cell + 1]
    }

    /// 領域ごとの体積（単位高さあたりの面積, cm²）の解析値。
    pub fn volumes(&self) -> Vec<f64> {
        let mut volumes: Vec<f64> = self
            .layout
            .iter()
            .flat_map(|&t| self.types[t].volumes())
            .collect();
        if self.margin > 0.0 {
            let strip = self.pitch() * self.margin;
            volumes.extend(std::iter::repeat_n(strip, 4 * self.n));
            volumes.extend(std::iter::repeat_n(self.margin * self.margin, 4));
        }
        volumes
    }

    /// 点 start から方向 direction（単位ベクトル）に長さ length 進む直線を、領域ごとの区間に分ける。
    pub fn segments(&self, start: (f64, f64), direction: (f64, f64), length: f64) -> Vec<Segment> {
        // half は格子（余白を除く）の半辺
        let (pitch, half) = (self.pitch(), 0.5 * self.n as f64 * self.pitch());
        // 格子線 x, y = -W/2 + k p（余白があれば格子の外周も）と交わる点で、セルごとの区間に分ける
        let mut cuts = vec![0.0, length];
        for (origin, d) in [(start.0, direction.0), (start.1, direction.1)] {
            if d != 0.0 {
                for k in 0..=self.n {
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
            if self.margin > 0.0 && (x.abs() > half || y.abs() > half) {
                segments.push(Segment {
                    region: self.margin_region(x, y, half),
                    length: w[1] - w[0],
                });
                continue;
            }
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

    /// 余白の中の点 (x, y) を含む余白の領域。`half` は格子の半辺
    fn margin_region(&self, x: f64, y: f64, half: f64) -> usize {
        let (n, pitch) = (self.n, self.pitch());
        let base = self.cell_region_count();
        let column = (((x + half) / pitch).floor().max(0.0) as usize).min(n - 1);
        let row = (((half - y) / pitch).floor().max(0.0) as usize).min(n - 1);
        let corner = |k: usize| base + 4 * n + k;
        match (x < -half, x > half, y > half, y < -half) {
            (true, _, true, _) => corner(0),
            (_, true, true, _) => corner(1),
            (true, _, _, true) => corner(2),
            (_, true, _, true) => corner(3),
            (_, _, true, _) => base + column,
            (_, _, _, true) => base + n + column,
            (true, _, _, _) => base + 2 * n + row,
            _ => base + 3 * n + row,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frames_need_room() {
        let cell = || PinCell::with_sectors(1.26, vec![0.4, 0.6], 8).unwrap();
        assert!(cell().with_frame(0.0).is_err());
        assert!(cell().with_frame(0.63).is_err());
        // 枠の内側（半辺 0.43）の外接円 0.608 より大きい円は置けない
        assert!(cell().with_frame(0.2).is_err());
        assert!(PinCell::grid(1.26, 3).unwrap().with_frame(0.01).is_err());
    }

    #[test]
    fn frame_takes_the_cell_edge() {
        let (pitch, t) = (1.26, 0.02);
        let cell = PinCell::with_sectors(pitch, vec![0.4, 0.6], 8)
            .unwrap()
            .with_frame(t)
            .unwrap();
        assert_eq!(cell.regions(), 3 * 8 + 1);
        let volumes = cell.volumes();
        let frame = pitch * pitch - (pitch - 2.0 * t).powi(2);
        assert!((volumes[24] - frame).abs() < 1e-12, "{volumes:?}");
        assert!((volumes.iter().sum::<f64>() - pitch * pitch).abs() < 1e-12);
        assert_eq!(cell.region_at(0.62, 0.0), 24);
        assert_eq!(cell.region_at(0.0, -0.625), 24);
        // 枠の内側の四隅は、最も外側のリング
        assert_eq!(cell.region_at(0.6, 0.6) / 8, 2);
    }

    #[test]
    fn pin_cell_matches_a_one_cell_lattice() {
        let cell = PinCell::with_sectors(1.26, vec![0.3, 0.4096], 4).unwrap();
        assert_eq!(
            Lattice::pin_cell(cell.clone()),
            Lattice::new(vec![cell], vec![vec![0]]).unwrap()
        );
    }

    #[test]
    fn validates_radii() {
        assert!(PinCell::new(1.26, vec![0.2, 0.4]).is_ok());
        assert!(PinCell::new(1.26, vec![0.4, 0.2]).is_err());
        // 内接円（0.63）より外の円はセルの外周で切り取る。外接円（0.891）以上は空になる
        assert!(PinCell::new(1.26, vec![0.54, 0.7, 0.85]).is_ok());
        assert!(PinCell::new(1.26, vec![0.54, 0.9]).is_err());
        assert!(PinCell::new(1.26, vec![0.0]).is_err());
        assert!(PinCell::new(0.0, vec![]).is_err());
    }

    #[test]
    fn rings_beyond_the_inscribed_circle_are_clipped() {
        // 半辺 h = 1 のセルの、内接円（1）と外接円（√2）の間の半径 1.2 の円は、セルの外周で切り取られる
        let cell = PinCell::new(2.0, vec![0.5, 1.2]).unwrap();
        let v = cell.volumes();
        // 半径 r > h の円と正方形の共通部分は 4 [h² tan θc + r² (π/4 - θc)]（cos θc = h/r）
        let (h, r) = (1.0_f64, 1.2_f64);
        let theta = (h / r).acos();
        let clipped = 4.0 * (h * h * theta.tan() + r * r * (0.25 * PI - theta));
        assert!((v[1] - (clipped - PI * 0.25)).abs() < 1e-14, "{}", v[1]);
        assert!((v.iter().sum::<f64>() - 4.0).abs() < 1e-14);
        // 角の近くの点は最も外の領域、辺の中点の近くは切り取られたリング
        assert_eq!(cell.region_at(0.95, 0.95), 2);
        assert_eq!(cell.region_at(0.99, 0.0), 1);
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

    /// 2×2 の格子の周りに幅 0.25 の余白
    fn with_margin() -> Lattice {
        let lattice = two_by_two();
        Lattice::with_margin(lattice.types.clone(), vec![vec![0, 1], vec![1, 0]], 0.25).unwrap()
    }

    #[test]
    fn margin_regions_and_volumes() {
        let lattice = with_margin();
        assert_eq!(lattice.width(), 4.5);
        assert_eq!(lattice.margin(), 0.25);
        // セルの 6 領域の後に、上の帯（左から）、下の帯、左の帯（上から）、右の帯、角（左上、右上、左下、右下）
        assert_eq!(lattice.regions(), 6 + 4 * 2 + 4);
        let v = lattice.volumes();
        assert_eq!(v.len(), 18);
        for (region, volume) in v.iter().enumerate().skip(6) {
            let expected = if region < 14 { 2.0 * 0.25 } else { 0.25 * 0.25 };
            assert!((volume - expected).abs() < 1e-14, "{region}: {volume}");
        }
        let total: f64 = v.iter().sum();
        assert!((total - 4.5 * 4.5).abs() < 1e-12);
        // 余白の領域が接するセル（角は角のセル）。セルの領域はそのセル
        let cells: Vec<usize> = (0..18).map(|r| lattice.region_cell(r)).collect();
        assert_eq!(
            cells,
            vec![0, 0, 1, 2, 3, 3, 0, 1, 2, 3, 0, 2, 1, 3, 0, 1, 2, 3]
        );
    }

    #[test]
    fn segments_cross_the_margin() {
        let lattice = with_margin();
        // 上の行の中心を左から右へ: 左の帯, ピンの外側, ピンの内側, ピンの外側, 水, 右の帯
        let s = lattice.segments((-2.25, 1.0), (1.0, 0.0), 4.5);
        let regions: Vec<usize> = s.iter().map(|s| s.region).collect();
        assert_eq!(regions, vec![10, 1, 0, 1, 2, 12]);
        let lengths: Vec<f64> = s.iter().map(|s| s.length).collect();
        for (a, b) in lengths.iter().zip([0.25, 0.5, 1.0, 0.5, 2.0, 0.25]) {
            assert!((a - b).abs() < 1e-14, "{lengths:?}");
        }
        // 上の余白の中を左から右へ: 左上の角, 上の帯, 上の帯, 右上の角
        let s = lattice.segments((-2.25, 2.1), (1.0, 0.0), 4.5);
        let regions: Vec<usize> = s.iter().map(|s| s.region).collect();
        assert_eq!(regions, vec![14, 6, 7, 15]);
    }

    #[test]
    fn margin_must_not_be_negative() {
        let types = two_by_two().types.clone();
        for margin in [-0.1, f64::NAN, f64::INFINITY] {
            assert!(
                Lattice::with_margin(types.clone(), vec![vec![0, 1], vec![1, 0]], margin).is_err()
            );
        }
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

    #[test]
    fn sector_cell_regions_and_segments() {
        // 1, 4, 8 以外は扇形ごとの外側の面積が等しくならないので受け付けない
        assert!(PinCell::with_sectors(2.0, vec![0.5], 3).is_err());
        assert!(PinCell::with_sectors(2.0, vec![0.5], 0).is_err());
        let cell = PinCell::with_sectors(2.0, vec![0.5], 4).unwrap();
        assert_eq!(cell.regions(), 8);
        // 内側の円の第1象限、第2象限、…、外側の第1象限、…
        assert_eq!(cell.region_at(0.1, 0.1), 0);
        assert_eq!(cell.region_at(-0.1, 0.1), 1);
        assert_eq!(cell.region_at(-0.1, -0.1), 2);
        assert_eq!(cell.region_at(0.1, -0.1), 3);
        assert_eq!(cell.region_at(0.9, 0.9), 4);
        assert_eq!(cell.region_at(0.9, -0.9), 7);
        let v = cell.volumes();
        let disk = PI * 0.25;
        for k in 0..4 {
            assert!((v[k] - disk / 4.0).abs() < 1e-15);
            assert!((v[4 + k] - (4.0 - disk) / 4.0).abs() < 1e-14);
        }
        // y = 0.25 を左から右へ: 外側の第2象限、内側の第2、第1象限、外側の第1象限
        let s = cell.segments((-1.0, 0.25), (1.0, 0.0), 2.0);
        let regions: Vec<usize> = s.iter().map(|s| s.region).collect();
        assert_eq!(regions, vec![5, 1, 0, 4]);
        let half_chord = (0.25_f64 - 0.0625).sqrt();
        let lengths: Vec<f64> = s.iter().map(|s| s.length).collect();
        for (a, b) in
            lengths
                .iter()
                .zip([1.0 - half_chord, half_chord, half_chord, 1.0 - half_chord])
        {
            assert!((a - b).abs() < 1e-14, "{lengths:?}");
        }
        // 8 個の扇形は 45° ごと
        let cell = PinCell::with_sectors(2.0, vec![], 8).unwrap();
        assert_eq!(cell.region_at(0.9, 0.1), 0);
        assert_eq!(cell.region_at(0.1, 0.9), 1);
        assert_eq!(cell.volumes(), vec![0.5; 8]);
        // 1 個なら扇形のないセルと同じ
        assert_eq!(
            PinCell::with_sectors(2.0, vec![0.5], 1).unwrap(),
            PinCell::new(2.0, vec![0.5]).unwrap()
        );
    }

    #[test]
    fn line_through_the_center_is_cut_at_sector_boundaries() {
        // 中心を通る直線では、どの扇形の境界との交点も中心になる
        let cell = PinCell::with_sectors(2.0, vec![0.5], 4).unwrap();
        let s = cell.segments((-0.6, -0.8), (0.6, 0.8), 2.0);
        let inner: Vec<(usize, f64)> = s
            .iter()
            .filter(|s| s.region < 4)
            .map(|s| (s.region, s.length))
            .collect();
        assert_eq!(inner.len(), 2, "{s:?}");
        assert_eq!((inner[0].0, inner[1].0), (2, 0));
        assert!(
            inner.iter().all(|(_, l)| (l - 0.5).abs() < 1e-14),
            "{inner:?}"
        );
    }
}
