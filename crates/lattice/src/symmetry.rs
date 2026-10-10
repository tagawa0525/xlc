//! 格子の対称性（正方形の回転と鏡映）による、トラック、領域、セルの像。
//!
//! 正方形の 8 つの対称操作のうち、外周の境界条件、セルの配置と大きさと種類、トラック（端の点と区間の長さ）、領域の
//! 値（物質）を変えないものを1つずつ確かめ、成り立つ操作をすべて使う（`detect`。ADR-0112）。対称性は前提とせず、
//! 真空境界や位置で違う物質のある検証の体系では、成り立つ操作だけを使う。MOC は使う操作の像（`images`）で代表の
//! トラックと領域だけを解き、自己遮蔽は同じ操作で移り合うピンと構造材の塊を組にして1回だけ解く（ADR-0114）。

use std::f64::consts::PI;

use crate::geometry::{Cell, Lattice, Segment};
use crate::tracking::{Track, Tracking};

/// 対称性のエラー。
#[derive(Debug, thiserror::Error)]
pub enum SymmetryError {
    #[error("{given} region materials for {regions} regions")]
    Regions { given: usize, regions: usize },
    #[error(
        "track {track} starts or ends at a corner of the lattice; image tracks are looked up \
         by position along the boundary, which needs ends away from the corners"
    )]
    CornerTrack { track: usize },
    #[error("symmetry operation {0} does not exist (there are 8)")]
    UnknownOperation(usize),
    #[error("the symmetry operations {0:?} are not closed under composition")]
    NotClosed(Vec<usize>),
    #[error("the lattice is not symmetric under the {operation}: {reason}")]
    Broken {
        operation: &'static str,
        reason: String,
    },
}

/// 正方形の 8 つの対称操作。格子の中心を原点とするセルの座標 (u, v)（u は右、v は下向き）に掛ける行列
/// [a, b, c, d]: (u, v) → (a u + b v, c u + d v)。恒等、回転（3つ）、左右・上下・対角線・反対角線の鏡映
const OPERATIONS: [[i64; 4]; 8] = [
    [1, 0, 0, 1],
    [0, -1, 1, 0],
    [-1, 0, 0, -1],
    [0, 1, -1, 0],
    [-1, 0, 0, 1],
    [1, 0, 0, -1],
    [0, 1, 1, 0],
    [0, -1, -1, 0],
];

/// `OPERATIONS` の順の操作の名前（エラーの表示用）
const NAMES: [&str; 8] = [
    "identity",
    "rotation by 90 degrees",
    "rotation by 180 degrees",
    "rotation by 270 degrees",
    "mirror across the vertical axis",
    "mirror across the horizontal axis",
    "mirror across the diagonal",
    "mirror across the anti-diagonal",
];

/// 解くときに使う対称操作の組み合わせ（恒等を含み、合成で閉じた操作の集合。ADR-0112）。正方形では 10 通りある。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Symmetry(usize);

/// 10 通りの対称性（`OPERATIONS` の番号のビットの集合）と、その名前（ログとエラーの表示用）。`Symmetry` はこの番号
const SYMMETRIES: [(u8, &str); 10] = [
    (0b0000_0001, "none"),
    (0b0000_0101, "rotation by 180 degrees (1/2)"),
    (0b0001_0001, "mirror across the vertical axis (1/2)"),
    (0b0010_0001, "mirror across the horizontal axis (1/2)"),
    (0b0100_0001, "mirror across the diagonal (1/2)"),
    (0b1000_0001, "mirror across the anti-diagonal (1/2)"),
    (0b0011_0101, "mirrors across both axes (1/4)"),
    (0b1100_0101, "mirrors across both diagonals (1/4)"),
    (0b0000_1111, "rotations by 90 degrees (1/4)"),
    (0b1111_1111, "all eight operations (1/8)"),
];

impl Symmetry {
    /// 恒等だけ（対称性を使わない）
    pub const NONE: Self = Self(0);
    /// 左右と上下の鏡映（1/4。軸を入れ替えないので長方形の格子にも使える）
    pub const QUARTER: Self = Self(6);
    /// 8 つの操作すべて（1/8）
    pub const OCTANT: Self = Self(9);

    /// 操作 `operations`（`OPERATIONS` の番号。恒等は含めなくてよい）を使う対称性。合成で閉じていなければエラー
    pub fn new(operations: &[usize]) -> Result<Self, SymmetryError> {
        if let Some(&unknown) = operations.iter().find(|&&i| i >= OPERATIONS.len()) {
            return Err(SymmetryError::UnknownOperation(unknown));
        }
        let mask = operations.iter().fold(1u8, |m, &i| m | (1 << i));
        SYMMETRIES
            .iter()
            .position(|&(s, _)| s == mask)
            .map(Self)
            .ok_or_else(|| {
                SymmetryError::NotClosed((0..8).filter(|i| mask & (1 << i) != 0).collect())
            })
    }

    /// 使う操作（`OPERATIONS` の番号の順。恒等が最初）
    pub fn operations(self) -> Vec<usize> {
        (0..OPERATIONS.len())
            .filter(|&i| self.contains(i))
            .collect()
    }

    /// 操作 `operation` を使うか
    fn contains(self, operation: usize) -> bool {
        SYMMETRIES[self.0].0 & (1 << operation) != 0
    }

    /// 対称性の名前
    pub fn name(self) -> &'static str {
        SYMMETRIES[self.0].1
    }

    /// 使う操作の数
    pub fn order(self) -> usize {
        SYMMETRIES[self.0].0.count_ones() as usize
    }
}

/// 格子で成り立つ対称操作（`detect`）と、使う操作の像。恒等の像を必ず含む（`detect` と `restrict` だけが作る）。
#[derive(Debug, Clone, PartialEq)]
pub struct Detected {
    symmetry: Symmetry,
    /// 使う操作の像（`Symmetry::operations` の順。恒等が最初）
    images: Vec<Image>,
    /// 格子のセルの数
    cells: usize,
    /// 確かめたトラッキングの (領域、トラック、方位角) の数（`fits`）
    shape: (usize, usize, usize),
    /// [操作] の成り立たなかった理由（成り立つ操作と使わない操作は `None`）
    rejected: Vec<Option<String>>,
}

impl Detected {
    /// トラッキング `tracking` で求めたものとして使えるか（セル、領域、トラック、方位角の数が同じ）。別の格子やトラッキングの
    /// 像は、領域とトラックの対応が成り立たない
    pub fn fits(&self, tracking: &Tracking) -> bool {
        self.cells == tracking.lattice.cells() && self.shape == shape(tracking)
    }

    /// 使う対称性
    pub fn symmetry(&self) -> Symmetry {
        self.symmetry
    }

    /// 使う操作の像（恒等が最初）
    pub fn images(&self) -> &[Image] {
        &self.images
    }

    /// [セル] の、移り合うセルのうち最も小さい番号（組の代表）
    pub fn cell_representatives(&self) -> Vec<usize> {
        (0..self.cells)
            .map(|cell| {
                self.images
                    .iter()
                    .map(|image| image.cells[cell])
                    .fold(cell, usize::min)
            })
            .collect()
    }

    /// 成り立つ操作のうち、対称性 `symmetry` の操作だけを使う（検証のための強制。ADR-0112）。`symmetry` のどれかの操作が
    /// 成り立たなければ、その操作の名前と成り立たなかった理由を示すエラー
    pub fn restrict(&self, symmetry: Symmetry) -> Result<Self, SymmetryError> {
        let held = self.symmetry.operations();
        let mut images = Vec::with_capacity(symmetry.order());
        for operation in symmetry.operations() {
            match held.iter().position(|&h| h == operation) {
                Some(i) => images.push(self.images[i].clone()),
                None => {
                    return Err(SymmetryError::Broken {
                        operation: NAMES[operation],
                        reason: self.rejected[operation]
                            .clone()
                            .unwrap_or_else(|| "it was not checked".into()),
                    });
                }
            }
        }
        Ok(Self {
            symmetry,
            images,
            cells: self.cells,
            shape: self.shape,
            rejected: self.rejected.clone(),
        })
    }
}

/// 恒等以外の 7 つの操作を1つずつ確かめ（`images` の確かめと、領域ごとの値 `values` が移った先の領域の値と同じこと）、
/// 成り立つ操作をすべて使う対称性と、その像を返す（ADR-0112）。恒等は確かめずに使う。使った対称性と、成り立たなかった
/// 操作の理由を info のログに出す。成り立つ操作が合成で閉じていなければ（確かめ方の不整合）エラー
pub fn detect<T: PartialEq + std::fmt::Display>(
    tracking: &Tracking,
    values: &[T],
) -> Result<Detected, SymmetryError> {
    let regions = tracking.lattice.regions();
    if values.len() != regions {
        return Err(SymmetryError::Regions {
            given: values.len(),
            regions,
        });
    }
    let tolerance = tracking.lattice.coincident();
    // 像のトラックの索引は一度だけ作る。作れなければ（端が外周の角にあるトラック）、恒等以外の操作はどれも確かめられない
    let starts = StartIndex::new(tracking, tolerance);
    let mut held = Vec::new();
    let mut images = vec![identity_image(tracking)];
    let mut rejected = vec![None; OPERATIONS.len()];
    for (index, &name) in NAMES.iter().enumerate().skip(1) {
        let image =
            match &starts {
                Ok(starts) => image_of(tracking, starts, index, tolerance).and_then(|image| {
                    match image
                        .regions
                        .iter()
                        .enumerate()
                        .find(|&(r, &i)| values[r] != values[i])
                    {
                        Some((r, &i)) => Err(SymmetryError::Broken {
                            operation: name,
                            reason: format!(
                                "region {r} has `{}`, but its image {i} has `{}`",
                                values[r], values[i]
                            ),
                        }),
                        None => Ok(image),
                    }
                }),
                Err(error) => Err(SymmetryError::Broken {
                    operation: name,
                    reason: error.to_string(),
                }),
            };
        match image {
            Ok(image) => {
                held.push(index);
                images.push(image);
            }
            Err(error) => {
                tracing::info!(operation = name, %error, "symmetry operation does not hold");
                // 強制したときのエラー（`restrict`）は操作の名前を付けて示すので、理由だけを残す
                rejected[index] = Some(match error {
                    SymmetryError::Broken { reason, .. } => reason,
                    other => other.to_string(),
                });
            }
        }
    }
    let symmetry = Symmetry::new(&held)?;
    tracing::info!(
        symmetry = symmetry.name(),
        operations = symmetry.order(),
        "symmetry"
    );
    Ok(Detected {
        symmetry,
        images,
        cells: tracking.lattice.cells(),
        shape: shape(tracking),
        rejected,
    })
}

/// 対称性 `symmetry` のどの操作も、形状とトラックで成り立つことを確かめた像（領域の値は確かめない。ライブラリと
/// テストで、対称性を指定して解くとき）。どれかの操作が成り立たなければエラー
pub fn require(tracking: &Tracking, symmetry: Symmetry) -> Result<Detected, SymmetryError> {
    let operations = symmetry.operations();
    let mut verified = vec![identity_image(tracking)];
    verified.extend(images(tracking, &operations[1..])?);
    Ok(Detected {
        symmetry,
        images: verified,
        cells: tracking.lattice.cells(),
        shape: shape(tracking),
        rejected: vec![None; OPERATIONS.len()],
    })
}

/// トラッキングの (領域、トラック、方位角) の数
fn shape(tracking: &Tracking) -> (usize, usize, usize) {
    (
        tracking.lattice.regions(),
        tracking.tracks.len(),
        tracking.azimuths.len(),
    )
}

/// 恒等の像（確かめなくても成り立つ）
fn identity_image(tracking: &Tracking) -> Image {
    Image {
        name: NAMES[0],
        matrix: [1.0, 0.0, 0.0, 1.0],
        cells: (0..tracking.lattice.cells()).collect(),
        regions: (0..tracking.lattice.regions()).collect(),
        tracks: (0..tracking.tracks.len()).map(|t| (t, true)).collect(),
    }
}

/// 1つの対称操作による像。
#[derive(Debug, Clone, PartialEq)]
pub struct Image {
    /// 操作の名前（エラーの表示用）
    pub name: &'static str,
    /// 格子の中心を原点とする (x, y)（y は上向き。トラックの座標）に掛ける行列 [a, b, c, d]:
    /// (x, y) → (a x + b y, c x + d y)
    pub matrix: [f64; 4],
    /// [セル] の像のセル
    pub cells: Vec<usize>,
    /// [領域] の像の領域
    pub regions: Vec<usize>,
    /// [トラック] の像のトラックと、像のトラックの順方向が元のトラックの順方向の像と同じ向きか
    pub tracks: Vec<(usize, bool)>,
}

/// 操作 `operations`（`OPERATIONS` の番号）それぞれの像。どれか1つでも、境界条件、セルの配置と大きさと種類、
/// トラック（端の点と、像の領域ごとの区間の長さ）が移り合わなければエラー。領域の物質は確かめない（使う側で
/// 確かめる）。
pub fn images(tracking: &Tracking, operations: &[usize]) -> Result<Vec<Image>, SymmetryError> {
    // 像のトラックの索引は一度だけ作る
    let tolerance = tracking.lattice.coincident();
    let starts = StartIndex::new(tracking, tolerance)?;
    operations
        .iter()
        .map(|&index| image_of(tracking, &starts, index, tolerance))
        .collect()
}

/// 操作 `index` の像（`images` を参照）
fn image_of(
    tracking: &Tracking,
    starts: &StartIndex,
    index: usize,
    tolerance: f64,
) -> Result<Image, SymmetryError> {
    let lattice = &tracking.lattice;
    let (op, name) = (&OPERATIONS[index], NAMES[index]);
    let broken = |reason: String| SymmetryError::Broken {
        operation: name,
        reason,
    };
    check_boundaries(tracking, op).map_err(broken)?;
    let (cells, regions) = region_map(lattice, op).map_err(broken)?;
    let matrix = xy_matrix(op).map(|v| v as f64);
    // トラックの像: 端の点を移し、像の向きの方位角の、同じ始点のトラックを探す
    let apply = |(x, y): (f64, f64)| (matrix[0] * x + matrix[1] * y, matrix[2] * x + matrix[3] * y);
    let azimuths = &tracking.azimuths;
    let mut tracks = Vec::with_capacity(tracking.tracks.len());
    for (index, track) in tracking.tracks.iter().enumerate() {
        let phi = azimuths[track.azimuth].phi;
        let (dx, dy) = apply((phi.cos(), phi.sin()));
        // 方位角は (0, π) で持つので、像の向きが下向きなら逆向きのトラックになる
        let forward = dy > 0.0;
        let (start, angle) = if forward {
            (apply(track.start), dy.atan2(dx))
        } else {
            (apply(track.end), (-dy).atan2(-dx))
        };
        let image = azimuths
            .iter()
            .position(|a| (a.phi - angle).abs() < 1.0e-9)
            .and_then(|a| starts.find(tracking, a, start, tolerance))
            .ok_or_else(|| broken(format!("track {index} has no image track")))?;
        let (a, b) = (
            &azimuths[track.azimuth],
            &azimuths[tracking.tracks[image].azimuth],
        );
        let close = |x: f64, y: f64| (x - y).abs() <= 1.0e-9 * x.abs().max(y.abs());
        if !(close(a.weight, b.weight) && close(a.spacing, b.spacing)) {
            return Err(broken(format!(
                "track {index} and its image {image} differ in weight or spacing"
            )));
        }
        if !segments_match(track, &tracking.tracks[image], &regions, forward, tolerance) {
            return Err(broken(format!(
                "track {index} and its image {image} differ in segments"
            )));
        }
        tracks.push((image, forward));
    }
    Ok(Image {
        name,
        matrix,
        cells,
        regions,
        tracks,
    })
}

/// 方位角ごとの、トラックの始点を外周に沿った位置の順に並べた索引（像のトラックを二分探索で探す）。
struct StartIndex {
    /// [方位角] の (外周に沿った位置, トラックの番号) を位置の順に
    sorted: Vec<Vec<(f64, usize)>>,
    /// 格子の幅と高さ
    size: (f64, f64),
}

impl StartIndex {
    /// トラックの端（始点と終点。終点は逆向きのトラックの像の始点になる）は、外周の角から `tolerance` より離れている
    /// こと。`Tracking` はトラックが外周の角を通らないように引くが、索引はそれに頼る（角の点は二つの辺の上にあり、外周に
    /// 沿った位置が二通りになる。左下の角では位置が 0 と外周の全長に分かれる）ので、確かめてエラーにする
    fn new(tracking: &Tracking, tolerance: f64) -> Result<Self, SymmetryError> {
        let size = (tracking.lattice.width(), tracking.lattice.height());
        let near_corner = |(x, y): (f64, f64)| {
            (x.abs() - 0.5 * size.0).abs() < tolerance && (y.abs() - 0.5 * size.1).abs() < tolerance
        };
        if let Some(track) = tracking
            .tracks
            .iter()
            .position(|t| near_corner(t.start) || near_corner(t.end))
        {
            return Err(SymmetryError::CornerTrack { track });
        }
        let sorted = (0..tracking.azimuths.len())
            .map(|a| {
                let mut keys: Vec<(f64, usize)> = tracking
                    .azimuth_tracks(a)
                    .map(|t| (perimeter(tracking.tracks[t].start, size), t))
                    .collect();
                keys.sort_by(|p, q| p.0.total_cmp(&q.0));
                keys
            })
            .collect();
        Ok(Self { sorted, size })
    }

    /// 方位角 `azimuth` の、始点が `point` から `tolerance` 以内のトラック。始点の間隔は許容誤差よりずっと大きいので、
    /// 外周に沿った位置で最も近い2本だけを確かめればよい
    fn find(
        &self,
        tracking: &Tracking,
        azimuth: usize,
        point: (f64, f64),
        tolerance: f64,
    ) -> Option<usize> {
        let keys = &self.sorted[azimuth];
        let at = keys.partition_point(|&(key, _)| key < perimeter(point, self.size));
        let same = |t: usize| {
            let start = tracking.tracks[t].start;
            (start.0 - point.0).abs() < tolerance && (start.1 - point.1).abs() < tolerance
        };
        keys[at.saturating_sub(1)..keys.len().min(at + 1)]
            .iter()
            .map(|&(_, t)| t)
            .find(|&t| same(t))
    }
}

/// 外周の上の点 `p`（中心が原点、幅と高さ `size`）の、左下の角から反時計回りに外周に沿って測った位置。点は、最も近い
/// 辺の上にあるとみなす。点は外周の角にないこと（`StartIndex::new` が確かめる。角では位置が二通りになる）
fn perimeter(p: (f64, f64), (w, h): (f64, f64)) -> f64 {
    let (x, y) = (p.0 + 0.5 * w, p.1 + 0.5 * h);
    // 下、右、上、左の辺までの距離
    let distances = [y, w - x, h - y, x].map(f64::abs);
    let nearest = (1..4).fold(0, |best, i| {
        if distances[i] < distances[best] {
            i
        } else {
            best
        }
    });
    match nearest {
        0 => x,
        1 => w + y,
        2 => w + h + (w - x),
        _ => 2.0 * w + h + (h - y),
    }
}

/// (u, v)（v は下向き）の操作 `op` を、(x, y) = (u, -v) の行列にしたもの
fn xy_matrix(op: &[i64; 4]) -> [i64; 4] {
    [op[0], -op[1], -op[2], op[3]]
}

/// 操作 `op` による [セル] と [領域] の像。セルの大きさと種類が移り合わなければ理由
fn region_map(lattice: &Lattice, op: &[i64; 4]) -> Result<(Vec<usize>, Vec<usize>), String> {
    let cells = (0..lattice.cells())
        .map(|cell| cell_image_checked(lattice, op, cell))
        .collect::<Result<Vec<_>, _>>()?;
    let m = xy_matrix(op);
    let mut regions = Vec::with_capacity(lattice.regions());
    for (cell, &image) in cells.iter().enumerate() {
        let offset = lattice.cell_regions(image).start;
        let local = cell_region_images(lattice.cell_type(cell), &m);
        regions.extend(local.into_iter().map(|l| offset + l));
    }
    Ok((cells, regions))
}

/// トラック `track` の区間を領域の像 `regions` で移したものが、像のトラック `image` の区間と同じか。領域ごとの長さの和と、
/// 入る領域（`forward` が偽なら像のトラックの終わりの領域）を比べる（長さの許容誤差 `tolerance`）。境界の上をかすめる
/// `tolerance` 以下の区間は、像の側にだけ現れることがあるので比べない
fn segments_match(
    track: &Track,
    image: &Track,
    regions: &[usize],
    forward: bool,
    tolerance: f64,
) -> bool {
    let long = |s: &&Segment| s.length > tolerance;
    let entry = if forward {
        image.segments.iter().find(long)
    } else {
        image.segments.iter().rev().find(long)
    };
    let entry_matches = match (track.segments.iter().find(long), entry) {
        (Some(a), Some(b)) => regions[a.region] == b.region,
        (None, None) => true,
        _ => false,
    };
    // (領域, 長さ) を領域の順に並べ、同じ領域の長さを足して比べる
    let lengths = |pairs: &mut Vec<(usize, f64)>| -> Vec<(usize, f64)> {
        pairs.sort_by_key(|&(r, _)| r);
        let mut merged: Vec<(usize, f64)> = Vec::with_capacity(pairs.len());
        for &(r, l) in pairs.iter() {
            match merged.last_mut() {
                Some((last, sum)) if *last == r => *sum += l,
                _ => merged.push((r, l)),
            }
        }
        merged.retain(|&(_, l)| l > tolerance);
        merged
    };
    let moved = lengths(
        &mut track
            .segments
            .iter()
            .map(|s| (regions[s.region], s.length))
            .collect(),
    );
    let own = lengths(
        &mut image
            .segments
            .iter()
            .map(|s| (s.region, s.length))
            .collect(),
    );
    entry_matches
        && moved.len() == own.len()
        && moved
            .iter()
            .zip(&own)
            .all(|(a, b)| a.0 == b.0 && (a.1 - b.1).abs() <= tolerance)
}

/// セルの種類 `cell` の中の領域（セルの中の番号）の、(x, y) の行列 `m` による像（像のセルの中の番号）。セルの種類が
/// `m` で変わらない（ブロックは大きさを入れ替えた種類が像のセルにある）ことは確かめてある
fn cell_region_images(cell: &Cell, m: &[i64; 4]) -> Vec<usize> {
    // 上の行から、各行は左から番号を付けた nx × ny の長方形の、像（nx と ny が入れ替わることがある）
    let squares = |nx: usize, ny: usize| -> Vec<usize> {
        let swap = swaps_axes(m);
        let (image_nx, image_ny) = if swap { (ny, nx) } else { (nx, ny) };
        (0..nx * ny)
            .map(|k| {
                let (c, r) = ((k % nx) as i64, (k / nx) as i64);
                let (x, y) = (2 * c - (nx as i64 - 1), (ny as i64 - 1) - 2 * r);
                let (x, y) = (m[0] * x + m[1] * y, m[2] * x + m[3] * y);
                let (c, r) = (x + image_nx as i64 - 1, image_ny as i64 - 1 - y);
                (r / 2) as usize * image_nx + (c / 2) as usize
            })
            .collect()
    };
    match cell {
        Cell::Block { nx, ny } => squares(*nx, *ny),
        Cell::Pin(pin) if pin.is_grid() => squares(pin.grid_size(), pin.grid_size()),
        Cell::Pin(pin) => {
            // 扇形 k は x 軸から角度 [2πk/S, 2π(k+1)/S)。S は 1, 4, 8 なので、扇形の境界は正方形の対称操作で
            // 扇形の境界に移る。扇形の中央の角の像で、像の扇形を決める。リングと枠は変わらない
            let sectors = pin.sectors();
            (0..pin.regions())
                .map(|local| {
                    if pin.frame().is_some_and(|(frame, _)| frame == local) {
                        return local;
                    }
                    let (ring, sector) = (local / sectors, local % sectors);
                    let theta = 2.0 * PI * (sector as f64 + 0.5) / sectors as f64;
                    let (x, y) = (theta.cos(), theta.sin());
                    let (x, y) = (
                        m[0] as f64 * x + m[1] as f64 * y,
                        m[2] as f64 * x + m[3] as f64 * y,
                    );
                    let angle = y.atan2(x).rem_euclid(2.0 * PI);
                    let image = ((angle / (2.0 * PI) * sectors as f64) as usize).min(sectors - 1);
                    ring * sectors + image
                })
                .collect()
        }
    }
}

/// 操作 `op` が軸を入れ替えるか（u 軸の像が v 軸の上）
fn swaps_axes(op: &[i64; 4]) -> bool {
    op[0] == 0
}

/// (u, v) に操作 `op` を掛ける
fn transform(op: &[i64; 4], (u, v): (i64, i64)) -> (i64, i64) {
    (op[0] * u + op[1] * v, op[2] * u + op[3] * v)
}

/// セル `cell` の操作 `op` による像のセル（格子の外に出れば `None`）。セルの中心の座標は2倍して整数にする:
/// U = 2c − (nc − 1)、V = 2r − (nr − 1)
fn cell_image(lattice: &Lattice, op: &[i64; 4], cell: usize) -> Option<usize> {
    let (nc, nr) = (lattice.columns(), lattice.rows());
    let (oc, or) = (nc as i64 - 1, nr as i64 - 1);
    let (u, v) = transform(
        op,
        (2 * (cell % nc) as i64 - oc, 2 * (cell / nc) as i64 - or),
    );
    let (c, r) = (u + oc, v + or);
    let inside = c >= 0 && r >= 0 && c % 2 == 0 && r % 2 == 0;
    let (c, r) = (c / 2, r / 2);
    (inside && (c as usize) < nc && (r as usize) < nr).then(|| r as usize * nc + c as usize)
}

/// セル `cell` の操作 `op` による像のセルで、大きさと種類が移り合うもの（軸を入れ替える操作では幅と高さ、ブロックの
/// nx と ny も入れ替える）。格子の外に出るか、移り合わなければ理由
fn cell_image_checked(lattice: &Lattice, op: &[i64; 4], cell: usize) -> Result<usize, String> {
    let image = cell_image(lattice, op, cell)
        .ok_or_else(|| format!("cell {cell} falls outside the lattice"))?;
    match cell_mismatch(lattice, op, cell, image) {
        Some(what) => Err(format!(
            "cell {cell} and its image {image} differ in {what}"
        )),
        None => Ok(image),
    }
}

/// セル `cell` と、操作 `op` による像のセル `image` の大きさか種類が違えば、どちらか。軸を入れ替える操作では、
/// 幅と高さ、ブロックの nx と ny を入れ替えて比べる
fn cell_mismatch(
    lattice: &Lattice,
    op: &[i64; 4],
    cell: usize,
    image: usize,
) -> Option<&'static str> {
    let swap = swaps_axes(op);
    let (w, h) = lattice.cell_size(cell);
    let size = if swap { (h, w) } else { (w, h) };
    if size != lattice.cell_size(image) {
        return Some("size");
    }
    let same = match (lattice.cell_type(cell), lattice.cell_type(image)) {
        (Cell::Block { nx, ny }, Cell::Block { nx: x, ny: y }) if swap => nx == y && ny == x,
        (a, b) => a == b,
    };
    (!same).then_some("type")
}

/// 外周の辺の境界条件が、操作 `op` で移った先の辺の境界条件と同じか。違えば理由
fn check_boundaries(tracking: &Tracking, op: &[i64; 4]) -> Result<(), String> {
    let b = tracking.boundaries();
    // 辺の外向きの向き (u, v)（v は下向き）と境界条件
    let sides = [
        ((-1, 0), b.left),
        ((1, 0), b.right),
        ((0, -1), b.top),
        ((0, 1), b.bottom),
    ];
    for &(direction, condition) in &sides {
        let image = transform(op, direction);
        let mapped = sides.iter().find(|(d, _)| *d == image).map(|(_, c)| *c);
        if mapped != Some(condition) {
            return Err(format!(
                "the side facing {direction:?} has {condition:?}, its image {mapped:?}"
            ));
        }
    }
    Ok(())
}

/// 対称操作で移り合う組（軌道）。MOC の領域と CMFD のセルは、組の代表の分だけ状態を持つ（ADR-0112）。
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Orbits {
    /// 組の代表（番号の順）
    pub representatives: Vec<usize>,
    /// [要素] の代表の、`representatives` の中の位置
    pub index: Vec<usize>,
}

impl Orbits {
    /// どの要素も自分だけの組（対称性を使わない）
    pub fn identity(count: usize) -> Self {
        Self {
            representatives: (0..count).collect(),
            index: (0..count).collect(),
        }
    }

    /// [要素] の代表 `representative`（代表は自分自身が代表）から
    pub fn new(representative: &[usize]) -> Self {
        let representatives: Vec<usize> = (0..representative.len())
            .filter(|&c| representative[c] == c)
            .collect();
        let mut position = vec![0; representative.len()];
        for (i, &c) in representatives.iter().enumerate() {
            position[c] = i;
        }
        let index: Vec<usize> = representative.iter().map(|&r| position[r]).collect();
        Self {
            representatives,
            index,
        }
    }

    /// [要素] の値 `values` の、[代表] の組ごとの和（組の要素の数や、領域の軌道の体積）
    pub fn sums(&self, values: &[f64]) -> Vec<f64> {
        let mut sums = vec![0.0; self.representatives.len()];
        for (&i, v) in self.index.iter().zip(values) {
            sums[i] += v;
        }
        sums
    }

    /// 要素 `element` が代表なら、その `representatives` の中の位置
    pub fn position(&self, element: usize) -> Option<usize> {
        let i = self.index[element];
        (self.representatives[i] == element).then_some(i)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Lattice, PinCell};
    use crate::tracking::{Boundaries, Boundary, Tracking};

    /// 3×3 の格子。`water` のセルは水、ほかはピン（燃料と水の 2 領域）。物質の番号: 燃料 0、水 1
    fn lattice(water: &[usize]) -> (Lattice, Vec<usize>) {
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let box_ = PinCell::new(1.26, vec![]).unwrap();
        let layout: Vec<Vec<usize>> = (0..3)
            .map(|r| {
                (0..3)
                    .map(|c| usize::from(water.contains(&(3 * r + c))))
                    .collect()
            })
            .collect();
        let materials = (0..9)
            .flat_map(|cell| {
                if water.contains(&cell) {
                    vec![1]
                } else {
                    vec![0, 1]
                }
            })
            .collect();
        (Lattice::new(vec![pin, box_], layout).unwrap(), materials)
    }

    fn orbits(boundaries: Boundaries, materials: Option<Vec<usize>>) -> Vec<usize> {
        let (lattice, default) = lattice(&[4]);
        let tracking = Tracking::with_boundaries(&lattice, 8, 0.2, boundaries).unwrap();
        detect(&tracking, &materials.unwrap_or(default))
            .unwrap()
            .cell_representatives()
    }

    #[test]
    fn swapping_axes_turns_blocks() {
        // 3 × 3 のピンの中央に 2 × 1 のブロックの水。90° の回転や対角線の鏡映では中央の像が 1 × 2 のブロックになり、
        // 小さな長方形の並びが移り合わないので、軸を入れ替えない操作（左右、上下、180° の回転）だけが残る:
        // 角の 4 本は同じ組、辺の上下の 2 本と左右の 2 本は別の組
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let lattice = Lattice::rectilinear(
            vec![1.26; 3],
            vec![1.26; 3],
            vec![Cell::Pin(pin), Cell::Block { nx: 2, ny: 1 }],
            vec![vec![0, 0, 0], vec![0, 1, 0], vec![0, 0, 0]],
        )
        .unwrap();
        let materials: Vec<usize> = (0..9)
            .flat_map(|cell| if cell == 4 { vec![1, 1] } else { vec![0, 1] })
            .collect();
        let tracking = Tracking::new(&lattice, 8, 0.2).unwrap();
        assert_eq!(
            detect(&tracking, &materials)
                .unwrap()
                .cell_representatives(),
            vec![0, 1, 0, 3, 4, 3, 0, 1, 0]
        );
    }

    #[test]
    fn swapped_blocks_are_images() {
        // 3 × 3 のピンの上下の辺の中央に 2 × 1、左右の辺の中央に 1 × 2 のブロックの水。軸を入れ替える操作でも
        // ブロックが移り合うので、8 つの操作がすべて成り立つ: 角の 4 本と辺の 4 つのブロックがそれぞれ同じ組
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let lattice = Lattice::rectilinear(
            vec![1.26; 3],
            vec![1.26; 3],
            vec![
                Cell::Pin(pin),
                Cell::Block { nx: 2, ny: 1 },
                Cell::Block { nx: 1, ny: 2 },
            ],
            vec![vec![0, 1, 0], vec![2, 0, 2], vec![0, 1, 0]],
        )
        .unwrap();
        let materials: Vec<usize> = (0..9)
            .flat_map(|cell| {
                if cell % 2 == 1 {
                    vec![1, 1]
                } else {
                    vec![0, 1]
                }
            })
            .collect();
        let tracking = Tracking::new(&lattice, 8, 0.2).unwrap();
        assert_eq!(
            detect(&tracking, &materials)
                .unwrap()
                .cell_representatives(),
            vec![0, 1, 0, 1, 4, 1, 0, 1, 0]
        );
    }

    #[test]
    fn reflective_lattice_has_all_symmetries() {
        // 角の 4 本と辺の 4 本がそれぞれ同じ組、中央の水は単独
        assert_eq!(
            orbits(Boundaries::REFLECTIVE, None),
            vec![0, 1, 0, 1, 4, 1, 0, 1, 0]
        );
    }

    #[test]
    fn vacuum_sides_keep_the_diagonal_mirror() {
        // 右と下が真空: 左上から右下への対角線の鏡映だけが残る
        let boundaries = Boundaries {
            right: Boundary::Vacuum,
            bottom: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        assert_eq!(orbits(boundaries, None), vec![0, 1, 2, 1, 4, 5, 2, 5, 8]);
    }

    #[test]
    fn materials_break_the_symmetry() {
        // 右下の角のピンだけ燃料の物質が違う（番号 2）: その角を動かさない対角線の鏡映だけが残る
        let (_, mut materials) = lattice(&[4]);
        let last = materials.len() - 2;
        materials[last] = 2;
        assert_eq!(
            orbits(Boundaries::REFLECTIVE, Some(materials)),
            vec![0, 1, 2, 1, 4, 5, 2, 5, 8]
        );
    }

    #[test]
    fn region_materials_must_cover_the_lattice() {
        let (lattice, materials) = lattice(&[4]);
        let tracking = Tracking::new(&lattice, 8, 0.2).unwrap();
        assert!(detect(&tracking, &materials[1..]).is_err());
    }

    /// 3 × 3 のピン（燃料 0、水 1）の中央に、扇形 8 の水のセル（扇形ごとの物質 `pattern`）を置いた格子のトラックと領域の
    /// 物質
    fn patterned(pattern: [usize; 8], boundaries: Boundaries) -> (Tracking, Vec<usize>) {
        let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
        let sectored = PinCell::with_sectors(1.26, vec![], 8).unwrap();
        let layout = vec![vec![0, 0, 0], vec![0, 1, 0], vec![0, 0, 0]];
        let lattice = Lattice::new(vec![pin, sectored], layout).unwrap();
        let values = (0..9)
            .flat_map(|cell| {
                if cell == 4 {
                    pattern.to_vec()
                } else {
                    vec![0, 1]
                }
            })
            .collect();
        let tracking = Tracking::with_boundaries(&lattice, 10, 0.05, boundaries).unwrap();
        (tracking, values)
    }

    #[test]
    fn identity_holds_even_when_images_cannot_be_found() {
        // 像のトラックは始点を外周に沿った位置で探すので、端が外周の角にあるトラックがあると、恒等以外の操作は確かめられ
        // ない（外した理由としてログに出す）。恒等は確かめなくても成り立つので、判定は止まらず、対称性なしになる
        let (mut tracking, values) = patterned([2; 8], Boundaries::REFLECTIVE);
        let half = 0.5 * tracking.lattice.width();
        tracking.tracks[0].start = (-half, -half);
        let detected = detect(&tracking, &values).unwrap();
        assert_eq!(detected.symmetry(), Symmetry::NONE);
        assert_eq!(detected.cell_representatives(), (0..9).collect::<Vec<_>>());
    }

    #[test]
    fn forced_symmetry_reports_the_failing_operation() {
        // 強制した対称性は、判定の結果から取り出す。成り立たない操作があれば、その操作の名前と理由を示すエラー
        let right = Boundaries {
            right: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let (tracking, values) = patterned([2; 8], right);
        let detected = detect(&tracking, &values).unwrap();
        // 1/4 の操作のうち、番号の順で最初に成り立たないもの（180° 回転）を、名前を重ねずに示す
        let quarter = detected
            .restrict(Symmetry::QUARTER)
            .unwrap_err()
            .to_string();
        assert!(
            quarter.contains("rotation by 180 degrees")
                && quarter.contains("Vacuum")
                && quarter.matches("not symmetric").count() == 1,
            "{quarter}"
        );
        let mirror = Symmetry::new(&[5]).unwrap();
        let kept = detected.restrict(mirror).unwrap();
        assert_eq!(kept.symmetry(), mirror);
        assert_eq!(kept.images().len(), 2);
        assert_eq!(detected.restrict(Symmetry::NONE).unwrap().images().len(), 1);
    }

    #[test]
    fn detects_every_operation_that_holds() {
        let symmetry = |ops: &[usize]| Symmetry::new(ops).unwrap();
        let reflective = Boundaries::REFLECTIVE;
        // 右と下が真空（C5G7 の 1/4 炉心の形）: 左上から右下への対角線の鏡映だけ
        let corner = Boundaries {
            right: Boundary::Vacuum,
            bottom: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        // 右だけ真空: 上下の鏡映だけ
        let right = Boundaries {
            right: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let cases = [
            ([2; 8], reflective, Symmetry::OCTANT),
            ([2; 8], corner, symmetry(&[6])),
            ([2; 8], right, symmetry(&[5])),
            // 扇形が 1 つおきに違う（風車）: 回転だけ
            ([2, 3, 2, 3, 2, 3, 2, 3], reflective, symmetry(&[1, 2, 3])),
            // 対角線をはさんで同じ: 2 つの対角線の鏡映と 180° 回転
            ([2, 2, 3, 3, 2, 2, 3, 3], reflective, symmetry(&[2, 6, 7])),
            // 半周ごとに同じ: 180° 回転だけ
            ([2, 3, 4, 5, 2, 3, 4, 5], reflective, symmetry(&[2])),
            // 扇形がどれも違う: なし
            ([2, 3, 4, 5, 6, 7, 8, 9], reflective, Symmetry::NONE),
        ];
        for (pattern, boundaries, expected) in cases {
            let (tracking, values) = patterned(pattern, boundaries);
            let detected = detect(&tracking, &values).unwrap();
            assert_eq!(
                detected.symmetry(),
                expected,
                "{pattern:?}: {} instead of {}",
                detected.symmetry().name(),
                expected.name()
            );
            // 像は使う操作ごとに 1 つ（恒等を含む）
            assert_eq!(detected.images().len(), expected.order());
        }
        // 長方形（真ん中の列が幅 1.0）では、軸を入れ替える操作は成り立たない
        let lattice = quarter_lattice(true);
        let tracking = Tracking::new(&lattice, 8, 0.05).unwrap();
        let values = vec![0; lattice.regions()];
        assert_eq!(
            detect(&tracking, &values).unwrap().symmetry(),
            Symmetry::QUARTER
        );
    }

    #[test]
    fn symmetry_names_closed_sets_only() {
        // 恒等は含めなくてよい。合成で閉じない操作の集合（90° の回転だけ）と、ない操作の番号はエラー
        assert_eq!(Symmetry::new(&[]).unwrap(), Symmetry::NONE);
        assert_eq!(Symmetry::new(&[2, 4, 5]).unwrap(), Symmetry::QUARTER);
        assert!(matches!(
            Symmetry::new(&[1]),
            Err(SymmetryError::NotClosed(_))
        ));
        assert!(matches!(
            Symmetry::new(&[8]),
            Err(SymmetryError::UnknownOperation(8))
        ));
    }

    /// 1/4 対称だが 1/8 対称でない格子。`odd` なら 3 × 3 で、真ん中の列は幅 1.0 の 2 × 2 のブロック（左右の鏡映の軸が
    /// セルの中を通る長方形）。そうでなければ 4 × 4 で、扇形 8 のピンと 2 × 2 の水を左右と上下で対称に（対角線では
    /// 非対称に）置き、周りに余白を付けたもの（軸がセルの境を通る）
    fn quarter_lattice(odd: bool) -> Lattice {
        let pin = PinCell::with_sectors(1.26, vec![0.35, 0.4096], 8).unwrap();
        if odd {
            return Lattice::rectilinear(
                vec![1.26, 1.0, 1.26],
                vec![1.26; 3],
                vec![Cell::Pin(pin), Cell::Block { nx: 2, ny: 2 }],
                vec![vec![0, 1, 0]; 3],
            )
            .unwrap();
        }
        let water = PinCell::grid(1.26, 2).unwrap();
        let layout = vec![
            vec![0, 1, 1, 0],
            vec![0, 0, 0, 0],
            vec![0, 0, 0, 0],
            vec![0, 1, 1, 0],
        ];
        Lattice::with_margin(vec![pin, water], layout, 0.04).unwrap()
    }

    #[test]
    fn quarter_images_map_tracks_onto_tracks() {
        for odd in [true, false] {
            let lattice = quarter_lattice(odd);
            let tracking = Tracking::new(&lattice, 10, 0.05).unwrap();
            let images = images(&tracking, &Symmetry::QUARTER.operations()).unwrap();
            assert_eq!(images.len(), 4);
            let exact = lattice.volumes();
            let apply =
                |m: [f64; 4], p: (f64, f64)| (m[0] * p.0 + m[1] * p.1, m[2] * p.0 + m[3] * p.1);
            let near =
                |p: (f64, f64), q: (f64, f64)| (p.0 - q.0).abs() < 1e-9 && (p.1 - q.1).abs() < 1e-9;
            for image in &images {
                // 領域は同じ体積の領域に1対1で移る
                let mut hit = vec![false; exact.len()];
                for (r, &i) in image.regions.iter().enumerate() {
                    assert!((exact[r] - exact[i]).abs() < 1e-12, "region {r} -> {i}");
                    hit[i] = true;
                }
                assert!(hit.iter().all(|&h| h));
                for (t, track) in tracking.tracks.iter().enumerate() {
                    let (u, same) = image.tracks[t];
                    let other = &tracking.tracks[u];
                    let (start, end) = if same {
                        (other.start, other.end)
                    } else {
                        (other.end, other.start)
                    };
                    assert!(near(apply(image.matrix, track.start), start), "track {t}");
                    assert!(near(apply(image.matrix, track.end), end), "track {t}");
                    // 移した領域ごとの区間の長さは、像のトラックの区間の長さと同じ
                    let mut lengths = vec![0.0_f64; exact.len()];
                    for s in &track.segments {
                        lengths[image.regions[s.region]] += s.length;
                    }
                    for s in &other.segments {
                        lengths[s.region] -= s.length;
                    }
                    assert!(lengths.iter().all(|l| l.abs() < 1e-9), "track {t}");
                }
            }
        }
    }

    #[test]
    fn segment_lengths_reject_a_wrong_region_map() {
        // 正しい領域の像ではどのトラックも区間の長さが一致し、ピンの燃料と水を入れ替えた像では一致しないトラックがある
        let (lattice, _) = lattice(&[4]);
        let tracking = Tracking::new(&lattice, 8, 0.2).unwrap();
        let images = images(&tracking, &Symmetry::QUARTER.operations()).unwrap();
        let tolerance = 1.0e-9 * lattice.width();
        for image in &images {
            for (t, &(u, forward)) in image.tracks.iter().enumerate() {
                assert!(segments_match(
                    &tracking.tracks[t],
                    &tracking.tracks[u],
                    &image.regions,
                    forward,
                    tolerance
                ));
            }
        }
        let image = &images[1];
        let mut wrong = image.regions.clone();
        wrong.swap(0, 1);
        assert!(image.tracks.iter().enumerate().any(|(t, &(u, forward))| {
            !segments_match(
                &tracking.tracks[t],
                &tracking.tracks[u],
                &wrong,
                forward,
                tolerance,
            )
        }));
    }

    #[test]
    fn region_values_must_map_onto_themselves() {
        // 入力の物質の名前で、像の領域が同じ物質か確かめる
        let (lattice, _) = lattice(&[4]);
        let tracking = Tracking::new(&lattice, 8, 0.2).unwrap();
        let names: Vec<String> = (0..lattice.cells())
            .flat_map(|cell| {
                if cell == 4 {
                    vec!["water".to_string()]
                } else {
                    vec!["fuel".to_string(), "water".to_string()]
                }
            })
            .collect();
        assert!(
            detect(&tracking, &names)
                .unwrap()
                .restrict(Symmetry::QUARTER)
                .is_ok()
        );
        // 右上の角のピンの燃料だけ別の名前
        let mut changed = names.clone();
        changed[lattice.cell_regions(2).start] = "other".into();
        let error = detect(&tracking, &changed)
            .unwrap()
            .restrict(Symmetry::QUARTER)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("`other`") && error.contains("`fuel`"),
            "{error}"
        );
    }

    #[test]
    fn track_ends_at_a_corner_are_an_error() {
        // 像のトラックは外周に沿った位置で探すので、端が外周の角にあると位置が二通りになる。トラックは角を通らない前提を
        // 確かめ、外れたら見つからないとせず理由を示す
        let (lattice, _) = lattice(&[4]);
        let mut tracking = Tracking::new(&lattice, 8, 0.2).unwrap();
        let half = (0.5 * lattice.width(), 0.5 * lattice.height());
        tracking.tracks[0].start = (-half.0, -half.1);
        let error = images(&tracking, &Symmetry::QUARTER.operations()).unwrap_err();
        assert!(
            matches!(error, SymmetryError::CornerTrack { track: 0 }),
            "{error}"
        );
    }

    #[test]
    fn quarter_needs_every_mirror() {
        // 右だけ真空
        let (symmetric, _) = lattice(&[4]);
        let boundaries = Boundaries {
            right: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let tracking = Tracking::with_boundaries(&symmetric, 8, 0.2, boundaries).unwrap();
        assert!(images(&tracking, &Symmetry::QUARTER.operations()).is_err());
        // 左の真ん中だけ水
        let (asymmetric, _) = lattice(&[3]);
        let tracking = Tracking::new(&asymmetric, 8, 0.2).unwrap();
        assert!(images(&tracking, &Symmetry::QUARTER.operations()).is_err());
        // 左右が真空で上下が鏡面反射なら、1/8 対称でなくても 1/4 対称
        let boundaries = Boundaries {
            left: Boundary::Vacuum,
            right: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let tracking = Tracking::with_boundaries(&symmetric, 8, 0.2, boundaries).unwrap();
        assert!(images(&tracking, &Symmetry::QUARTER.operations()).is_ok());
    }
}
