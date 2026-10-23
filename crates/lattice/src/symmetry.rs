//! 格子の対称性（正方形の回転と鏡映）による、トラック、領域、セルの像。
//!
//! 正方形の 8 つの対称操作のうち、外周の境界条件、セルの配置と大きさと種類、トラック（端の点と区間の長さ）、領域の
//! 値（物質）を変えないものを1つずつ確かめ、成り立つ操作をすべて使う（`detect`。ADR-0112）。対称性は前提とせず、
//! 真空境界や位置で違う物質のある検証の体系では、成り立つ操作だけを使う。MOC は使う操作の像（`images`）で代表の
//! トラックと領域だけを解き、自己遮蔽は同じ操作で移り合うピンと構造材の塊を組にして1回だけ解く（ADR-0112）。

use std::f64::consts::PI;

use rayon::prelude::*;

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
        // 成り立ったが使わない操作は、その理由を残す（もう一度絞り込んだときに「確かめていない」としない）
        let rejected = (0..OPERATIONS.len())
            .map(|operation| {
                if held.contains(&operation) && !symmetry.contains(operation) {
                    Some(format!(
                        "it holds, but the symmetry \"{}\" leaves it out",
                        symmetry.name()
                    ))
                } else {
                    self.rejected[operation].clone()
                }
            })
            .collect();
        Ok(Self {
            symmetry,
            images,
            cells: self.cells,
            shape: self.shape,
            rejected,
        })
    }

    /// 対称性を使わない（恒等だけ。確かめない）
    pub fn none(tracking: &Tracking) -> Self {
        Self {
            symmetry: Symmetry::NONE,
            images: vec![identity_image(tracking)],
            cells: tracking.lattice.cells(),
            shape: shape(tracking),
            rejected: (0..OPERATIONS.len())
                .map(|i| (i != 0).then(|| "the symmetry \"none\" does not check it".to_string()))
                .collect(),
        }
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

/// 自己遮蔽と MOC に使う対称操作（`choose`）。
#[derive(Debug, Clone, PartialEq)]
pub struct Chosen {
    shielding: Detected,
    /// MOC に使う対称性が自己遮蔽と違うとき（強制したとき）だけ
    moc: Option<Detected>,
}

impl Chosen {
    /// 自己遮蔽のピンと塊の組に使う対称操作（ADR-0112）
    pub fn shielding(&self) -> &Detected {
        &self.shielding
    }

    /// MOC に使う対称操作
    pub fn moc(&self) -> &Detected {
        self.moc.as_ref().unwrap_or(&self.shielding)
    }
}

/// 自己遮蔽と MOC に使う対称操作を選ぶ（ADR-0112）。`forced` が `None` なら `detect` の結果を両方に使う。
/// 強制した対称性が `Symmetry::NONE` なら、判定をせず、両方とも対称性なしにする（検証の基準なので判定の結果に頼らない）。
/// ほかの強制は、`detect` の結果を自己遮蔽に、その結果から取り出した強制の対称性を MOC に使い、成り立たない操作があれば
/// その名前と理由を示すエラー
pub fn choose<T: PartialEq + std::fmt::Display>(
    tracking: &Tracking,
    values: &[T],
    forced: Option<Symmetry>,
) -> Result<Chosen, SymmetryError> {
    match forced {
        Some(Symmetry::NONE) => {
            tracing::info!("symmetry \"none\": no symmetry is detected or used");
            Ok(Chosen {
                shielding: Detected::none(tracking),
                moc: None,
            })
        }
        Some(symmetry) => {
            let detected = detect(tracking, values)?;
            let moc = detected.restrict(symmetry)?;
            tracing::info!(symmetry = symmetry.name(), "forced symmetry");
            Ok(Chosen {
                shielding: detected,
                moc: Some(moc),
            })
        }
        None => Ok(Chosen {
            shielding: detect(tracking, values)?,
            moc: None,
        }),
    }
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
    // トラックごとに並列に確かめ、成り立たなければ番号の最も小さいトラックの理由を返す（エラーが分け方によらない）
    let checked: Vec<Result<(usize, bool), SymmetryError>> = tracking
        .tracks
        .par_iter()
        .enumerate()
        .map(|(index, track)| {
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
            Ok((image, forward))
        })
        .collect();
    let tracks = checked.into_iter().collect::<Result<Vec<_>, _>>()?;
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
mod tests;
