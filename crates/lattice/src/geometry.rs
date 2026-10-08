//! 格子計算の形状: 同心円で区切った正方形のセル（ピンセル）と均質な長方形のセルを、列ごとの幅と
//! 行ごとの高さを持つ長方形の格子に並べたもの。

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

    /// 外周の枠の厚さ (cm)
    pub fn frame(&self) -> Option<f64> {
        self.frame
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
        // 枠の内側の正方形の辺 x, y = ±a
        if self.frame.is_some() {
            let a = self.inner_half();
            for (origin, d) in [(start.0, direction.0), (start.1, direction.1)] {
                if d != 0.0 {
                    for edge in [-a, a] {
                        let t = (edge - origin) / d;
                        if t > 0.0 && t < length {
                            cuts.push(t);
                        }
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
                    if t > 0.0
                        && t < length
                        && u.0 * point.0 + u.1 * point.1 >= 0.0
                        && !self.in_frame(point.0, point.1)
                    {
                        cuts.push(t);
                    }
                }
            }
        }
        for &r in &self.radii {
            let disc = b * b - (c0 - r * r);
            if disc > 0.0 {
                let root = disc.sqrt();
                // 枠の中の交点（枠の内側の正方形の外接円の外）では区切らない
                cuts.extend([-b - root, -b + root].into_iter().filter(|&t| {
                    t > 0.0
                        && t < length
                        && !self.in_frame(start.0 + t * direction.0, start.1 + t * direction.1)
                }));
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
        if self.frame.is_some() {
            volumes.push(self.pitch * self.pitch - 4.0 * a * a);
        }
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
                    _ => {}
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

    /// 列の幅（左から, cm）
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

    /// 点 start から方向 direction（単位ベクトル）に長さ length 進む直線を、領域ごとの区間に分ける。
    pub fn segments(&self, start: (f64, f64), direction: (f64, f64), length: f64) -> Vec<Segment> {
        let (xs, ys) = (&self.xs, &self.ys);
        let mut cuts = vec![0.0, length];
        for (lines, origin, d) in [(xs, start.0, direction.0), (ys, start.1, direction.1)] {
            if d != 0.0 {
                for &line in lines.iter() {
                    let t = (line - origin) / d;
                    if t > 0.0 && t < length {
                        cuts.push(t);
                    }
                }
            }
        }
        cuts.sort_by(f64::total_cmp);
        let columns = self.columns();
        let mut segments = Vec::new();
        for w in cuts.windows(2).filter(|w| w[1] > w[0]) {
            let mid = 0.5 * (w[0] + w[1]);
            let (x, y) = (start.0 + mid * direction.0, start.1 + mid * direction.1);
            let column = xs[1..].partition_point(|&e| e <= x).min(columns - 1);
            let row = ys[1..].partition_point(|&e| e >= y).min(self.rows() - 1);
            let cell = row * columns + column;
            let center = (
                0.5 * (xs[column] + xs[column + 1]),
                0.5 * (ys[row] + ys[row + 1]),
            );
            let local = (
                start.0 + w[0] * direction.0 - center.0,
                start.1 + w[0] * direction.1 - center.1,
            );
            let offset = self.offsets[cell];
            let pieces = match self.cell_type(cell) {
                Cell::Pin(pin) => pin.segments(local, direction, w[1] - w[0]),
                Cell::Block { nx, ny } => block_segments(
                    (self.widths[column], self.heights[row]),
                    (*nx, *ny),
                    local,
                    direction,
                    w[1] - w[0],
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
) -> Vec<Segment> {
    let (dx, dy) = (w / nx as f64, h / ny as f64);
    let mut cuts = vec![0.0, length];
    for (origin, d, step, count, low) in [
        (start.0, direction.0, dx, nx, -0.5 * w),
        (start.1, direction.1, dy, ny, -0.5 * h),
    ] {
        if d != 0.0 {
            for k in 1..count {
                let t = (low + k as f64 * step - origin) / d;
                if t > 0.0 && t < length {
                    cuts.push(t);
                }
            }
        }
    }
    cuts.sort_by(f64::total_cmp);
    cuts.windows(2)
        .filter(|w| w[1] > w[0])
        .map(|win| {
            let mid = 0.5 * (win[0] + win[1]);
            let (x, y) = (start.0 + mid * direction.0, start.1 + mid * direction.1);
            let column = (((x + 0.5 * w) / dx).floor().max(0.0) as usize).min(nx - 1);
            let row = (((0.5 * h - y) / dy).floor().max(0.0) as usize).min(ny - 1);
            Segment {
                region: row * nx + column,
                length: win[1] - win[0],
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adjacent_regions_of_rings_sectors_and_frames() {
        // 半ピッチ 0.5（枠の内側 0.45）。リング: r < 0.3、0.3〜0.5、0.5 より外（4 扇形）と枠
        let pin = PinCell::with_sectors(1.0, vec![0.3, 0.5], 4)
            .unwrap()
            .with_frame(0.05)
            .unwrap();
        let mut pairs = Cell::Pin(pin).adjacent();
        pairs.sort_unstable();
        let mut expected = Vec::new();
        for s in 0..4 {
            let next = (s + 1) % 4;
            // 同じ扇形の内外のリングは接する
            for ring in 0..2 {
                expected.push((4 * ring + s, 4 * (ring + 1) + s));
            }
            // リング 0 と 1 は軸の上で隣の扇形と接する（内径 0.3 < 0.45）。リング 2（内径 0.5 > 0.45）は
            // 軸の上になく、扇形どうしは接しない
            for ring in 0..2 {
                expected.push((4 * ring + s.min(next), 4 * ring + s.max(next)));
            }
            // リング 1（外径 0.5 > 0.45）とリング 2 は枠と接する
            expected.extend([(4 + s, 12), (8 + s, 12)]);
        }
        expected.sort_unstable();
        assert_eq!(pairs, expected);
    }

    #[test]
    fn eight_sectors_meet_across_diagonals_beyond_the_edge() {
        // 半ピッチ 0.5。リング 1（内径 0.6 > 0.5）は軸の上になく、対角線（0.5√2 ≈ 0.71 まで）の上でだけ
        // 隣の扇形と接する。扇形 0 と 1 の境界は対角線（45°）、1 と 2 の境界は軸（90°）
        let pin = PinCell::with_sectors(1.0, vec![0.6], 8).unwrap();
        let pairs = Cell::Pin(pin).adjacent();
        assert!(pairs.contains(&(8, 9)), "{pairs:?}");
        assert!(!pairs.contains(&(9, 10)), "{pairs:?}");
        assert!(
            pairs.contains(&(1, 2)) && pairs.contains(&(0, 7)),
            "{pairs:?}"
        );
    }

    #[test]
    fn adjacent_squares_of_grids_and_blocks() {
        let mut pairs = Cell::Pin(PinCell::grid(1.0, 2).unwrap()).adjacent();
        pairs.sort_unstable();
        assert_eq!(pairs, vec![(0, 1), (0, 2), (1, 3), (2, 3)]);
        let mut pairs = Cell::Block { nx: 3, ny: 1 }.adjacent();
        pairs.sort_unstable();
        assert_eq!(pairs, vec![(0, 1), (1, 2)]);
    }

    #[test]
    fn frames_need_room() {
        let cell = || PinCell::with_sectors(1.26, vec![0.4, 0.6], 8).unwrap();
        assert!(cell().with_frame(0.0).is_err());
        assert!(cell().with_frame(0.63).is_err());
        // 枠の内側（半辺 0.38）の外接円 0.537 より大きい円（0.6）は置けない
        assert!(cell().with_frame(0.25).is_err());
        assert!(matches!(
            PinCell::grid(1.26, 3).unwrap().with_frame(0.01),
            Err(GeometryError::FrameOnGrid)
        ));
    }

    #[test]
    fn frame_is_not_cut_inside() {
        // 枠の中では扇形や円で区切らない（同じ領域の区間が続かない）
        let cell = PinCell::with_sectors(1.26, vec![0.4, 0.75], 8)
            .unwrap()
            .with_frame(0.05)
            .unwrap();
        let d = (0.6_f64, 0.8_f64);
        let segments = cell.segments((-0.63, -0.3), d, 2.0);
        assert!(
            segments.windows(2).all(|w| w[0].region != w[1].region),
            "{segments:?}"
        );
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
        let pin = PinCell::new(2.0, vec![0.5]).unwrap();
        let water = PinCell::new(2.0, vec![]).unwrap();
        Lattice::with_margin(vec![pin, water], vec![vec![0, 1], vec![1, 0]], 0.25).unwrap()
    }

    #[test]
    fn margin_regions_and_volumes() {
        let lattice = with_margin();
        assert_eq!(lattice.width(), 4.5);
        // 余白は周りの薄い行と列: 4 × 4 のセル。領域はセルの順（行優先）で、ピンのセルは 2 領域
        assert_eq!((lattice.columns(), lattice.rows()), (4, 4));
        assert_eq!(lattice.widths(), &[0.25, 2.0, 2.0, 0.25]);
        assert_eq!(lattice.regions(), 6 + 4 * 2 + 4);
        let v = lattice.volumes();
        let corner = 0.25 * 0.25;
        let strip = 2.0 * 0.25;
        for region in [0, 3, 14, 17] {
            assert!((v[region] - corner).abs() < 1e-14, "{region}");
        }
        for region in [1, 2, 4, 8, 9, 13, 15, 16] {
            assert!((v[region] - strip).abs() < 1e-14, "{region}");
        }
        let total: f64 = v.iter().sum();
        assert!((total - 4.5 * 4.5).abs() < 1e-12);
        let cells: Vec<usize> = (0..18).map(|r| lattice.region_cell(r)).collect();
        assert_eq!(
            cells,
            vec![0, 1, 2, 3, 4, 5, 5, 6, 7, 8, 9, 10, 10, 11, 12, 13, 14, 15]
        );
    }

    #[test]
    fn segments_cross_the_margin() {
        let lattice = with_margin();
        // 上のピンの行の中心を左から右へ: 左の余白, ピンの外側, ピンの内側, ピンの外側, 水, 右の余白
        let s = lattice.segments((-2.25, 1.0), (1.0, 0.0), 4.5);
        let regions: Vec<usize> = s.iter().map(|s| s.region).collect();
        assert_eq!(regions, vec![4, 6, 5, 6, 7, 8]);
        let lengths: Vec<f64> = s.iter().map(|s| s.length).collect();
        for (a, b) in lengths.iter().zip([0.25, 0.5, 1.0, 0.5, 2.0, 0.25]) {
            assert!((a - b).abs() < 1e-14, "{lengths:?}");
        }
        // 上の余白の中を左から右へ: 左上の角, 上の帯, 上の帯, 右上の角
        let s = lattice.segments((-2.25, 2.1), (1.0, 0.0), 4.5);
        let regions: Vec<usize> = s.iter().map(|s| s.region).collect();
        assert_eq!(regions, vec![0, 1, 2, 3]);
    }

    #[test]
    fn margin_must_not_be_negative() {
        let types = vec![PinCell::new(2.0, vec![0.5]).unwrap()];
        for margin in [-0.1, f64::NAN, f64::INFINITY] {
            assert!(Lattice::with_margin(types.clone(), vec![vec![0]], margin).is_err());
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
