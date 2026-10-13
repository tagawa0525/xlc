//! 格子計算の形状: 同心円で区切った正方形のセル（ピンセル）と均質な長方形のセルを、列ごとの幅と
//! 行ごとの高さを持つ長方形の格子に並べたもの。

use std::f64::consts::PI;

/// 格子の大きさに対する、同じ点とみなす距離の比。座標の丸め（格子の大きさの ε の程度）より十分大きく、列と行の幅
/// （格子の作成で、これより広いことを確かめる）やトラックの間隔より十分小さい
const COINCIDENT: f64 = 1.0e-9;

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
    /// 外周の正方形の枠の厚さ (cm)
    frame: Option<f64>,
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
    #[error("column widths and row heights must be positive and finite")]
    Sizes,
    #[error("a column or row of {size} cm is not wider than the same-point distance {limit} cm")]
    ThinCell { size: f64, limit: f64 },
    #[error(
        "{what} of {size} cm in cell type {cell_type} is not wider than the same-point distance {limit} cm"
    )]
    ThinRegion {
        what: &'static str,
        cell_type: usize,
        size: f64,
        limit: f64,
    },
    #[error("{rows} rows of cells for {heights} row heights")]
    Rows { rows: usize, heights: usize },
    #[error("row {row} has {cells} cells for {columns} column widths")]
    Columns {
        row: usize,
        cells: usize,
        columns: usize,
    },
    #[error("the pin cell at row {row}, column {column} needs a square slot of its pitch")]
    PinSlot { row: usize, column: usize },
    #[error("frame thickness must be positive, below pitch/2 and outside the circles, got {0}")]
    Frame(f64),
    #[error("a frame needs a cell with rings, not a grid of squares")]
    FrameOnGrid,
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
            frame: None,
        })
    }

    pub fn pitch(&self) -> f64 {
        self.pitch
    }

    pub fn radii(&self) -> &[f64] {
        &self.radii
    }

    /// 外周の枠の (セルの中の領域の番号, 面積 (cm²))
    pub fn frame(&self) -> Option<(usize, f64)> {
        self.frame.map(|t| {
            (
                self.regions() - 1,
                self.pitch * self.pitch - (self.pitch - 2.0 * t).powi(2),
            )
        })
    }

    /// 各リングを分ける扇形の数
    pub fn sectors(&self) -> usize {
        self.sectors
    }

    /// 1辺の小さな正方形の数（同心円のあるセルでは 1）
    pub fn grid_size(&self) -> usize {
        self.grid
    }

    /// 同心円ではなく小さな正方形に分けたセルか
    pub fn is_grid(&self) -> bool {
        self.grid > 1
    }

    pub fn regions(&self) -> usize {
        if self.grid > 1 {
            self.grid * self.grid
        } else {
            (self.radii.len() + 1) * self.sectors + usize::from(self.frame.is_some())
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
        if self.grid > 1 {
            return Err(GeometryError::FrameOnGrid);
        }
        let inner = 0.5 * self.pitch - thickness;
        let fits = thickness > 0.0
            && inner > 0.0
            && self
                .radii
                .last()
                .is_none_or(|&r| r < inner * std::f64::consts::SQRT_2);
        if !fits {
            return Err(GeometryError::Frame(thickness));
        }
        Ok(Self {
            frame: Some(thickness),
            ..self
        })
    }

    /// 同心円のセルの、面を共有する領域の組（[`Cell::adjacent`]）。同じ扇形の内外のリング、扇形の境界の半直線
    /// （x 軸から等角度）がリングを通る隣の扇形、枠の内側の正方形の辺に届くリングと枠が接する。
    fn adjacent(&self) -> Vec<(usize, usize)> {
        let (rings, sectors) = (self.radii.len() + 1, self.sectors);
        let a = self.inner_half();
        let region = |ring: usize, sector: usize| ring * sectors + sector;
        let inner = |ring: usize| if ring == 0 { 0.0 } else { self.radii[ring - 1] };
        let mut pairs = Vec::new();
        for ring in 0..rings {
            for s in 0..sectors {
                if ring + 1 < rings {
                    pairs.push((region(ring, s), region(ring + 1, s)));
                }
                // 扇形 s と s + 1 の境界の半直線は角度 2π(s + 1)/sectors で、正方形の辺までの長さは軸の上なら a、
                // 対角線の上なら a√2
                if sectors > 1 {
                    let next = (s + 1) % sectors;
                    let on_axis = (s + 1) * 4 % sectors == 0;
                    let reach = if on_axis {
                        a
                    } else {
                        a * std::f64::consts::SQRT_2
                    };
                    if inner(ring) < reach {
                        pairs.push((region(ring, s.min(next)), region(ring, s.max(next))));
                    }
                }
                // どの扇形も軸の向き（辺までの長さ a）を含む（扇形は 1, 4, 8 個だけ）ので、外径が a を超える
                // リングは辺に届く
                let reaches_edge = ring + 1 == rings || self.radii[ring] > a;
                if self.frame.is_some() && reaches_edge {
                    pairs.push((region(ring, s), rings * sectors));
                }
            }
        }
        pairs
    }

    /// セルの中の区切りの幅（名前と幅 cm）: 格子の正方形、最も内のリングの半径とリングの幅、枠の厚さ
    fn widths(&self) -> Vec<(&'static str, f64)> {
        if self.is_grid() {
            return vec![("a grid square", self.pitch / self.grid as f64)];
        }
        let mut widths: Vec<(&'static str, f64)> = self
            .radii
            .iter()
            .scan(0.0, |inner, &r| {
                let width = r - *inner;
                *inner = r;
                Some(("a ring", width))
            })
            .collect();
        // 円が枠の内側の正方形やセルの辺に接しても、間の領域は四隅の塊で細い帯ではない（細いのは接点のまわりだけで、
        // 接線に近い弦と同じく、まとめて失う体積は同じ点とみなす距離の 2 乗の程度）ので確かめない
        if let Some(thickness) = self.frame {
            widths.push(("a frame", thickness));
        }
        widths
    }

    /// 枠の内側の正方形の半辺
    fn inner_half(&self) -> f64 {
        0.5 * self.pitch - self.frame.unwrap_or(0.0)
    }

    /// 点 (x, y) が枠の中か
    fn in_frame(&self, x: f64, y: f64) -> bool {
        let a = self.inner_half();
        self.frame.is_some() && (x.abs() > a || y.abs() > a)
    }

    /// 点 (x, y) を含む領域。
    pub fn region_at(&self, x: f64, y: f64) -> usize {
        if self.grid > 1 {
            let (m, h) = (self.grid, self.pitch / self.grid as f64);
            let index = |v: f64| ((v / h).floor().max(0.0) as usize).min(m - 1);
            return index(0.5 * self.pitch - y) * m + index(x + 0.5 * self.pitch);
        }
        if self.in_frame(x, y) {
            return (self.radii.len() + 1) * self.sectors;
        }
        let ring = self.radii.partition_point(|&ri| ri <= x.hypot(y));
        let angle = y.atan2(x).rem_euclid(2.0 * PI);
        let sector = ((angle / (2.0 * PI) * self.sectors as f64) as usize).min(self.sectors - 1);
        ring * self.sectors + sector
    }

    /// セルの中心を原点とする点 start から方向 direction に長さ length 進む直線を、
    /// 同心円で区間に分ける。直線に沿って `coincident` 以内の区切りは同じ点とみてまとめる（[`Lattice::coincident`]）。
    pub(crate) fn segments(
        &self,
        start: (f64, f64),
        direction: (f64, f64),
        length: f64,
        coincident: f64,
    ) -> Vec<Segment> {
        // 区切りの t（範囲の外のものは `cuts` が捨てる）。|start + t d|² = r² の解は −b ± √(r² − p²)
        let b = start.0 * direction.0 + start.1 * direction.1;
        let mut crossings = Vec::new();
        // 小さな正方形の境界の線 x, y = -p/2 + k p/m
        let h = self.pitch / self.grid as f64;
        for (origin, d) in [(start.0, direction.0), (start.1, direction.1)] {
            if d != 0.0 {
                for k in 1..self.grid {
                    crossings.push((-0.5 * self.pitch + k as f64 * h - origin) / d);
                }
            }
        }
        // 枠の内側の正方形の辺 x, y = ±a
        if self.frame.is_some() {
            let a = self.inner_half();
            for (origin, d) in [(start.0, direction.0), (start.1, direction.1)] {
                if d != 0.0 {
                    crossings.extend([-a, a].map(|edge| (edge - origin) / d));
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
                    if u.0 * point.0 + u.1 * point.1 >= 0.0 && !self.in_frame(point.0, point.1) {
                        crossings.push(t);
                    }
                }
            }
        }
        // 中心から直線までの距離。|start + t d|² = r² の判別式 b² − (|start|² − r²) は、接線に近い直線では b² と
        // |start|² の打ち消し合いで丸めの誤差が大きいので、同じ値 r² − p² を (r − p)(r + p) で求める
        let p = (start.0 * direction.1 - start.1 * direction.0).abs();
        for &r in &self.radii {
            let disc = (r - p) * (r + p);
            if disc > 0.0 {
                let root = disc.sqrt();
                // 枠の中の交点（枠の内側の正方形の外接円の外）では区切らない
                crossings.extend([-b - root, -b + root].into_iter().filter(|&t| {
                    !self.in_frame(start.0 + t * direction.0, start.1 + t * direction.1)
                }));
            }
        }
        cuts(crossings, length, coincident)
            .windows(2)
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
        let a = self.inner_half();
        let disk = |r: f64| clipped_disk_area(r, a);
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
        volumes.push(4.0 * a * a - disk(inner));
        let mut volumes: Vec<f64> = volumes
            .iter()
            .flat_map(|v| std::iter::repeat_n(v / self.sectors as f64, self.sectors))
            .collect();
        volumes.extend(self.frame().map(|(_, area)| area));
        volumes
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

/// 格子のマスに置くセル: ピンセル（幅と高さがそのピッチのマスだけ）か、均質な長方形を `nx` × `ny` の小さな
/// 長方形に分けたもの（大きさはマスで決まる）。
#[derive(Debug, Clone, PartialEq)]
pub enum Cell {
    Pin(PinCell),
    Block { nx: usize, ny: usize },
}

impl Cell {
    /// 均質な長方形を分けないもの
    pub const BLOCK: Self = Self::Block { nx: 1, ny: 1 };

    pub fn regions(&self) -> usize {
        match self {
            Self::Pin(pin) => pin.regions(),
            Self::Block { nx, ny } => nx * ny,
        }
    }

    /// セルの中で面を共有する領域の組（セルの中の番号、小さい順）
    pub fn adjacent(&self) -> Vec<(usize, usize)> {
        match self {
            Self::Pin(pin) if pin.is_grid() => squares(pin.grid, pin.grid),
            Self::Pin(pin) => pin.adjacent(),
            Self::Block { nx, ny } => squares(*nx, *ny),
        }
    }

    /// ピンセルなら、そのピンセル
    pub fn pin(&self) -> Option<&PinCell> {
        match self {
            Self::Pin(pin) => Some(pin),
            Self::Block { .. } => None,
        }
    }
}

/// 上の行から、各行は左から番号を付けた `nx` × `ny` の長方形のうち、辺を共有する組
fn squares(nx: usize, ny: usize) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for row in 0..ny {
        for column in 0..nx {
            let k = row * nx + column;
            if column + 1 < nx {
                pairs.push((k, k + 1));
            }
            if row + 1 < ny {
                pairs.push((k, k + nx));
            }
        }
    }
    pairs
}

/// 列ごとの幅と行ごとの高さを持つ長方形の格子にセルを並べた形状（中心が原点）。
///
/// マス（セル）は上の行から順に、各行は左から順に番号を付ける。領域はセルの順に、セルの中では
/// `PinCell` の領域の順（均質な長方形では上の行から、各行は左から）に通し番号を付ける。
/// 集合体の周りの余白（水のギャップ）は、薄い行と列の均質なセルで表す（`with_margin`）。
#[derive(Debug, Clone, PartialEq)]
pub struct Lattice {
    /// 列の幅（左から, cm）
    widths: Vec<f64>,
    /// 行の高さ（上から, cm）
    heights: Vec<f64>,
    /// 列の境界の x（左から）と、行の境界の y（上から）。両端を含む
    xs: Vec<f64>,
    ys: Vec<f64>,
    types: Vec<Cell>,
    /// [セル] のセルの種類（`types` の番号）
    layout: Vec<usize>,
    /// [セル] のセルの最初の領域の番号（最後に領域の総数）
    offsets: Vec<usize>,
}

impl Lattice {
    /// ピンセル `types`（すべて同じピッチ）を正方格子に並べた形状。`layout` は上の行から順の各行のセルの種類。
    pub fn new(types: Vec<PinCell>, layout: Vec<Vec<usize>>) -> Result<Self, GeometryError> {
        Self::with_margin(types, layout, 0.0)
    }

    /// 正方格子（`new` を参照）の周りに幅 `margin` (cm) の余白を付けた形状。余白は上下左右の薄い行と列の
    /// 均質なセル（ピンセルの列・行ごとと四隅）になる。
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
        let border = types.len();
        let mut cells: Vec<Cell> = types.into_iter().map(Cell::Pin).collect();
        if margin == 0.0 {
            return Self::rectilinear(vec![pitch; n], vec![pitch; n], cells, layout);
        }
        cells.push(Cell::BLOCK);
        let mut sizes = vec![margin];
        sizes.extend(std::iter::repeat_n(pitch, n));
        sizes.push(margin);
        let edge = vec![border; n + 2];
        let mut rows = vec![edge.clone()];
        rows.extend(layout.into_iter().map(|row| {
            std::iter::once(border)
                .chain(row)
                .chain(std::iter::once(border))
                .collect()
        }));
        rows.push(edge);
        Self::rectilinear(sizes.clone(), sizes, cells, rows)
    }

    /// 列の幅 `widths`（左から）と行の高さ `heights`（上から）の格子に、セルの種類 `types` を `layout`
    /// （上の行から順の各行のセルの種類）のとおりに並べた形状。ピンセルは、幅と高さがそのピッチに等しいマスに置く。
    pub fn rectilinear(
        widths: Vec<f64>,
        heights: Vec<f64>,
        types: Vec<Cell>,
        layout: Vec<Vec<usize>>,
    ) -> Result<Self, GeometryError> {
        if layout.is_empty() || layout.iter().all(Vec::is_empty) {
            return Err(GeometryError::NoCells);
        }
        let positive = |v: &[f64]| !v.is_empty() && v.iter().all(|&x| x > 0.0 && x.is_finite());
        if !positive(&widths) || !positive(&heights) {
            return Err(GeometryError::Sizes);
        }
        // 同じ点とみなす距離（`Lattice::coincident`）以下の列や行は、両側の線を横切る点がまとまり、区間を持てない
        let limit = COINCIDENT * widths.iter().sum::<f64>().max(heights.iter().sum::<f64>());
        if let Some(&size) = widths.iter().chain(&heights).find(|&&x| x <= limit) {
            return Err(GeometryError::ThinCell { size, limit });
        }
        if types.is_empty() {
            return Err(GeometryError::NoCellTypes);
        }
        if layout.len() != heights.len() {
            return Err(GeometryError::Rows {
                rows: layout.len(),
                heights: heights.len(),
            });
        }
        if let Some(row) = layout.iter().position(|row| row.len() != widths.len()) {
            return Err(GeometryError::Columns {
                row: row + 1,
                cells: layout[row].len(),
                columns: widths.len(),
            });
        }
        if let Some(&bad) = layout.iter().flatten().find(|&&t| t >= types.len()) {
            return Err(GeometryError::CellType {
                index: bad,
                types: types.len(),
            });
        }
        for (r, row) in layout.iter().enumerate() {
            for (c, &t) in row.iter().enumerate() {
                match &types[t] {
                    Cell::Pin(pin) if widths[c] != pin.pitch() || heights[r] != pin.pitch() => {
                        return Err(GeometryError::PinSlot {
                            row: r + 1,
                            column: c + 1,
                        });
                    }
                    Cell::Block { nx, ny } if *nx == 0 || *ny == 0 => {
                        return Err(GeometryError::EmptyGrid);
                    }
                    Cell::Block { nx, ny } => {
                        let parts = [
                            ("a block column", widths[c] / *nx as f64),
                            ("a block row", heights[r] / *ny as f64),
                        ];
                        thin_region(t, &parts, limit)?;
                    }
                    Cell::Pin(pin) => thin_region(t, &pin.widths(), limit)?,
                }
            }
        }
        let layout: Vec<usize> = layout.into_iter().flatten().collect();
        let mut offsets = Vec::with_capacity(layout.len() + 1);
        let mut next = 0;
        for &t in &layout {
            offsets.push(next);
            next += types[t].regions();
        }
        offsets.push(next);
        let xs = edges(&widths, 1.0);
        let ys = edges(&heights, -1.0);
        Ok(Self {
            widths,
            heights,
            xs,
            ys,
            types,
            layout,
            offsets,
        })
    }

    /// 1つのピンセルだけの形状（`new` と同じものを、検査なしで直接作る）。
    pub fn pin_cell(cell: PinCell) -> Self {
        let (pitch, regions) = (cell.pitch(), cell.regions());
        Self {
            widths: vec![pitch],
            heights: vec![pitch],
            xs: vec![-0.5 * pitch, 0.5 * pitch],
            ys: vec![0.5 * pitch, -0.5 * pitch],
            types: vec![Cell::Pin(cell)],
            layout: vec![0],
            offsets: vec![0, regions],
        }
    }

    /// 同じ点とみなす距離 (cm): 格子の大きさの `COINCIDENT` 倍。トラックの端どうしをつなぐときと、直線が線を横切る
    /// 点をまとめるときに使う
    pub fn coincident(&self) -> f64 {
        COINCIDENT * self.width().max(self.height())
    }

    /// 列の幅（左から, cm）
    /// 列の境界の x（左から、両端を含む）。区間を作るときの線と同じもの
    pub(crate) fn column_lines(&self) -> &[f64] {
        &self.xs
    }

    /// 行の境界の y（上から、両端を含む）
    pub(crate) fn row_lines(&self) -> &[f64] {
        &self.ys
    }

    pub fn widths(&self) -> &[f64] {
        &self.widths
    }

    /// 行の高さ（上から, cm）
    pub fn heights(&self) -> &[f64] {
        &self.heights
    }

    /// 列の数
    pub fn columns(&self) -> usize {
        self.widths.len()
    }

    /// 行の数
    pub fn rows(&self) -> usize {
        self.heights.len()
    }

    /// セルの数
    pub fn cells(&self) -> usize {
        self.layout.len()
    }

    /// 横の長さ (cm)
    pub fn width(&self) -> f64 {
        self.widths.iter().sum()
    }

    /// 縦の長さ (cm)
    pub fn height(&self) -> f64 {
        self.heights.iter().sum()
    }

    pub fn regions(&self) -> usize {
        self.offsets[self.layout.len()]
    }

    /// 領域の属するセル
    pub fn region_cell(&self, region: usize) -> usize {
        self.offsets.partition_point(|&o| o <= region) - 1
    }

    /// セル `cell` の種類
    pub fn cell_type(&self, cell: usize) -> &Cell {
        &self.types[self.layout[cell]]
    }

    /// セル `cell` の領域の範囲。
    pub fn cell_regions(&self, cell: usize) -> std::ops::Range<usize> {
        self.offsets[cell]..self.offsets[cell + 1]
    }

    /// ピンセルのある行（上から）と列（左から）の番号。余白や反射体のように、ピッチに合わないマスだけの
    /// 行と列を除く（ピン出力の表の範囲）
    pub fn pin_rows_and_columns(&self) -> (Vec<usize>, Vec<usize>) {
        let columns = self.columns();
        let pin = |r: usize, c: usize| self.cell_type(r * columns + c).pin().is_some();
        let rows = (0..self.rows())
            .filter(|&r| (0..columns).any(|c| pin(r, c)))
            .collect();
        let kept = (0..columns)
            .filter(|&c| (0..self.rows()).any(|r| pin(r, c)))
            .collect();
        (rows, kept)
    }

    /// セル `cell` の幅と高さ
    pub fn cell_size(&self, cell: usize) -> (f64, f64) {
        let columns = self.columns();
        (self.widths[cell % columns], self.heights[cell / columns])
    }

    /// 領域ごとの体積（単位高さあたりの面積, cm²）の解析値。
    pub fn volumes(&self) -> Vec<f64> {
        (0..self.cells())
            .flat_map(|cell| match self.cell_type(cell) {
                Cell::Pin(pin) => pin.volumes(),
                Cell::Block { nx, ny } => {
                    let (w, h) = self.cell_size(cell);
                    vec![w * h / (nx * ny) as f64; nx * ny]
                }
            })
            .collect()
    }

    /// 点 (x, y)（格子の中心を原点とする）を含む領域。格子の外の点は、いちばん近い端のセルとする
    pub fn region_at(&self, x: f64, y: f64) -> usize {
        let (column, row) = self.cell_at(x, y);
        let cell = row * self.columns() + column;
        let (cx, cy) = self.cell_center(column, row);
        let (lx, ly) = (x - cx, y - cy);
        let local = match self.cell_type(cell) {
            Cell::Pin(pin) => pin.region_at(lx, ly),
            Cell::Block { nx, ny } => block_region(
                (self.widths[column], self.heights[row]),
                (*nx, *ny),
                (lx, ly),
            ),
        };
        self.offsets[cell] + local
    }

    /// 点 (x, y)（格子の中心を原点とする）を含むセルの（列、行（上から））。格子の外は、いちばん近い端のセル
    fn cell_at(&self, x: f64, y: f64) -> (usize, usize) {
        let column = self.xs[1..]
            .partition_point(|&e| e <= x)
            .min(self.columns() - 1);
        let row = self.ys[1..]
            .partition_point(|&e| e >= y)
            .min(self.rows() - 1);
        (column, row)
    }

    /// セル（列、行）の中心
    fn cell_center(&self, column: usize, row: usize) -> (f64, f64) {
        (
            0.5 * (self.xs[column] + self.xs[column + 1]),
            0.5 * (self.ys[row] + self.ys[row + 1]),
        )
    }

    /// 点 start から方向 direction（単位ベクトル）に長さ length 進む直線を、領域ごとの区間に分ける。
    pub fn segments(&self, start: (f64, f64), direction: (f64, f64), length: f64) -> Vec<Segment> {
        let (xs, ys) = (&self.xs, &self.ys);
        let mut crossings = Vec::new();
        for (lines, origin, d) in [(xs, start.0, direction.0), (ys, start.1, direction.1)] {
            if d != 0.0 {
                crossings.extend(lines.iter().map(|&line| (line - origin) / d));
            }
        }
        let coincident = self.coincident();
        let cuts = cuts(crossings, length, coincident);
        let columns = self.columns();
        let mut segments = Vec::new();
        for w in cuts.windows(2) {
            let mid = 0.5 * (w[0] + w[1]);
            let (x, y) = (start.0 + mid * direction.0, start.1 + mid * direction.1);
            let (column, row) = self.cell_at(x, y);
            let cell = row * columns + column;
            let center = self.cell_center(column, row);
            let local = (
                start.0 + w[0] * direction.0 - center.0,
                start.1 + w[0] * direction.1 - center.1,
            );
            let offset = self.offsets[cell];
            let pieces = match self.cell_type(cell) {
                Cell::Pin(pin) => pin.segments(local, direction, w[1] - w[0], coincident),
                Cell::Block { nx, ny } => block_segments(
                    (self.widths[column], self.heights[row]),
                    (*nx, *ny),
                    local,
                    direction,
                    w[1] - w[0],
                    coincident,
                ),
            };
            segments.extend(pieces.into_iter().map(|s| Segment {
                region: offset + s.region,
                length: s.length,
            }));
        }
        segments
    }
}

/// 長さ `length` の直線の区切り（両端 0 と `length` を含む狭義の昇順。`length` が正でなければ空）。直線に沿って
/// 隣どうしが `coincident` 以内の交点は同じ点で横切ったとみて、組の中点の 1 点にまとめる（中点は進む向きによらない。
/// 逆向きに通る像のトラックでも区間が同じになる）。端 0 と `length` を含む組は端にまとめ、両端どうしはまとめない。
/// 2 つの線（縦と横の線、円と辺、扇形の境界と円など）の交わる点をちょうど通る直線でも、それぞれを横切る t は丸めで
/// 少しずれる。まとめないと、交わる点のまわりのどれかの領域に丸めの程度の区間が残り、その領域を経由したことになる
/// （格子のセルは #0051、ピンとブロックのセルの中は #0055）。範囲の外の交点は捨てる
fn cuts(crossings: Vec<f64>, length: f64, coincident: f64) -> Vec<f64> {
    if length.is_nan() || length <= 0.0 {
        return Vec::new();
    }
    let mut points: Vec<f64> = crossings
        .into_iter()
        .filter(|&t| t > -coincident && t < length + coincident)
        .collect();
    points.sort_by(f64::total_cmp);
    let mut cuts = vec![0.0];
    let mut i = 0;
    while i < points.len() {
        // 隣どうしが coincident 以内の組 points[i..j]
        let mut j = i + 1;
        while j < points.len() && points[j] - points[j - 1] <= coincident {
            j += 1;
        }
        let (low, high) = (points[i], points[j - 1]);
        if low > coincident && high < length - coincident {
            cuts.push(0.5 * (low + high));
        }
        i = j;
    }
    cuts.push(length);
    cuts
}

/// 大きさ (w, h) を (nx, ny) に分けたブロックの、セルの中心を原点とする点 (x, y) を含む領域（行は上から）
fn block_region((w, h): (f64, f64), (nx, ny): (usize, usize), (x, y): (f64, f64)) -> usize {
    let column = (((x + 0.5 * w) / (w / nx as f64)).floor().max(0.0) as usize).min(nx - 1);
    let row = (((0.5 * h - y) / (h / ny as f64)).floor().max(0.0) as usize).min(ny - 1);
    row * nx + column
}

/// セルの種類 `cell_type` の中の区切りの幅 `parts`（名前と幅）のどれかが、同じ点とみなす距離 `limit` 以下ならエラー。
/// 区切りがまとまり、その領域が区間を持てなくなる（斜めに浅く横切るトラックでは区間が残り、体積が黙って過小になる）
fn thin_region(
    cell_type: usize,
    parts: &[(&'static str, f64)],
    limit: f64,
) -> Result<(), GeometryError> {
    match parts.iter().find(|(_, size)| *size <= limit) {
        Some(&(what, size)) => Err(GeometryError::ThinRegion {
            what,
            cell_type,
            size,
            limit,
        }),
        None => Ok(()),
    }
}

/// 大きさ `sizes` を中心を原点に並べた境界の座標（向き `sign` に。両端を含み、両端はちょうど ∓sign × 全長/2。
/// 足し合わせの丸めで端がずれると、トラックの端に極小の区間ができるため）
fn edges(sizes: &[f64], sign: f64) -> Vec<f64> {
    let half = 0.5 * sizes.iter().sum::<f64>();
    let mut edges = Vec::with_capacity(sizes.len() + 1);
    let mut at = -sign * half;
    edges.push(at);
    for &s in &sizes[..sizes.len().saturating_sub(1)] {
        at += sign * s;
        edges.push(at);
    }
    edges.push(sign * half);
    edges
}

/// 幅と高さ `(w, h)` の均質な長方形を `(nx, ny)` に分けたセルの中の区間（中心が原点。領域は上の行から、
/// 各行は左から）。
fn block_segments(
    (w, h): (f64, f64),
    (nx, ny): (usize, usize),
    start: (f64, f64),
    direction: (f64, f64),
    length: f64,
    coincident: f64,
) -> Vec<Segment> {
    let (dx, dy) = (w / nx as f64, h / ny as f64);
    let mut crossings = Vec::new();
    for (origin, d, step, count, low) in [
        (start.0, direction.0, dx, nx, -0.5 * w),
        (start.1, direction.1, dy, ny, -0.5 * h),
    ] {
        if d != 0.0 {
            crossings.extend((1..count).map(|k| (low + k as f64 * step - origin) / d));
        }
    }
    cuts(crossings, length, coincident)
        .windows(2)
        .map(|win| {
            let mid = 0.5 * (win[0] + win[1]);
            let (x, y) = (start.0 + mid * direction.0, start.1 + mid * direction.1);
            Segment {
                region: block_region((w, h), (nx, ny), (x, y)),
                length: win[1] - win[0],
            }
        })
        .collect()
}

#[cfg(test)]
mod tests;
