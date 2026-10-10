//! 格子の対称性（正方形の回転と鏡映）による、ピンの組と、トラックと領域の像。
//!
//! 正方形の 8 つの対称操作のうち、セルの配置と大きさ、領域の物質、外周の境界条件を変えないものを調べ、それで
//! 移り合うセルを同じ組にする。組の中のピンは、自己遮蔽の結果が同じなので1回だけ解けばよい（VERA の集合体の
//! ように 1/8 対称なら約 1/8 になる）。対称性は前提とせず、真空境界や位置で違う物質のある検証の体系では、成り立つ
//! 操作だけを使う。
//!
//! MOC で対称性を使うとき（ADR-0113）は、指定した操作のすべてで形状とトラックが移り合うことを確かめ、操作ごとの
//! セル、領域、トラックの像を求める（`images`）。

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

/// 1/4 対称の操作（`OPERATIONS` の番号）: 恒等、180° 回転、左右と上下の鏡映。軸を入れ替えないので、長方形の
/// 格子にも使える（ADR-0113）
pub const QUARTER: [usize; 4] = [0, 2, 4, 5];

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
    let lattice = &tracking.lattice;
    // 区間の端の一致を確かめる許容誤差（`Tracking` の反射のつながりと同じ）
    let tolerance = 1.0e-9 * lattice.width().max(lattice.height());
    let starts = StartIndex::new(tracking, tolerance)?;
    let mut images = Vec::with_capacity(operations.len());
    for &index in operations {
        let (op, name) = (&OPERATIONS[index], NAMES[index]);
        let broken = |reason: String| SymmetryError::Broken {
            operation: name,
            reason,
        };
        check_boundaries(tracking, op).map_err(broken)?;
        let cells = (0..lattice.cells())
            .map(|cell| cell_image_checked(lattice, op, cell))
            .collect::<Result<Vec<_>, _>>()
            .map_err(broken)?;
        // (u, v)（v は下向き）の行列を (x, y) = (u, -v) に直す
        let m = [op[0], -op[1], -op[2], op[3]];
        let matrix = m.map(|v| v as f64);
        let mut regions = Vec::with_capacity(lattice.regions());
        for (cell, &image) in cells.iter().enumerate() {
            let offset = lattice.cell_regions(image).start;
            let local = cell_region_images(lattice.cell_type(cell), &m);
            regions.extend(local.into_iter().map(|l| offset + l));
        }
        // トラックの像: 端の点を移し、像の向きの方位角の、同じ始点のトラックを探す
        let apply =
            |(x, y): (f64, f64)| (matrix[0] * x + matrix[1] * y, matrix[2] * x + matrix[3] * y);
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
        images.push(Image {
            name,
            matrix,
            cells,
            regions,
            tracks,
        });
    }
    Ok(images)
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

/// 領域ごとの値 `values`（入力の物質の名前など）が、像 `images` のどれでも移った先の領域の値と同じか。違えば、最初に
/// 違う領域とその像の値を示すエラー
pub fn check_regions<T: PartialEq + std::fmt::Display>(
    images: &[Image],
    values: &[T],
) -> Result<(), SymmetryError> {
    for image in images {
        if values.len() != image.regions.len() {
            return Err(SymmetryError::Regions {
                given: values.len(),
                regions: image.regions.len(),
            });
        }
        if let Some((r, &i)) = image
            .regions
            .iter()
            .enumerate()
            .find(|&(r, &i)| values[r] != values[i])
        {
            return Err(SymmetryError::Broken {
                operation: image.name,
                reason: format!(
                    "region {r} has `{}`, but its image {i} has `{}`",
                    values[r], values[i]
                ),
            });
        }
    }
    Ok(())
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
    let swap = swaps_axes(op);
    let (w, h) = lattice.cell_size(cell);
    let size = if swap { (h, w) } else { (w, h) };
    if size != lattice.cell_size(image) {
        return Err(format!("cell {cell} and its image {image} differ in size"));
    }
    let turned = match lattice.cell_type(cell) {
        Cell::Block { nx, ny } if swap => Cell::Block { nx: *ny, ny: *nx },
        other => other.clone(),
    };
    if &turned != lattice.cell_type(image) {
        return Err(format!("cell {cell} and its image {image} differ in type"));
    }
    Ok(image)
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

/// 対称操作で移り合う組（軌道）。MOC の領域と CMFD のセルは、組の代表の分だけ状態を持つ（ADR-0113）。
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

/// セルごとの、同じ組の代表のセル（組の中で最も小さい番号）。`region_material` は領域ごとの物質の番号。
pub fn cell_orbits(
    tracking: &Tracking,
    region_material: &[usize],
) -> Result<Vec<usize>, SymmetryError> {
    let lattice = &tracking.lattice;
    let regions = lattice.regions();
    if region_material.len() != regions {
        return Err(SymmetryError::Regions {
            given: region_material.len(),
            regions,
        });
    }
    let cells = lattice.cells();
    // セルのリングごとの物質（扇形や小さな正方形で物質が違うセルは None: 回すと領域の対応が変わる）
    let rings = |cell: usize| -> Option<Vec<usize>> {
        let materials = &region_material[lattice.cell_regions(cell)];
        let width = match lattice.cell_type(cell) {
            Cell::Pin(pin) if !pin.is_grid() => pin.sectors(),
            _ => materials.len(),
        };
        materials
            .chunks(width)
            .map(|ring| ring.iter().all(|&m| m == ring[0]).then_some(ring[0]))
            .collect()
    };
    let ring_materials: Vec<Option<Vec<usize>>> = (0..cells).map(rings).collect();
    let valid: Vec<&[i64; 4]> = OPERATIONS
        .iter()
        .filter(|op| {
            let identity = **op == OPERATIONS[0];
            // セルの種類は、軸を入れ替える操作でもそのまま比べる（ブロックの nx と ny を入れ替えない。#0053）
            let cells_match = (0..cells).all(|cell| {
                cell_image(lattice, op, cell).is_some_and(|image| {
                    let (w, h) = lattice.cell_size(cell);
                    let size = if swaps_axes(op) { (h, w) } else { (w, h) };
                    size == lattice.cell_size(image)
                        && lattice.cell_type(cell) == lattice.cell_type(image)
                        && ring_materials[cell].is_some()
                        && ring_materials[cell] == ring_materials[image]
                })
            });
            identity || (check_boundaries(tracking, op).is_ok() && cells_match)
        })
        .collect();
    Ok((0..cells)
        .map(|cell| {
            valid
                .iter()
                .filter_map(|op| cell_image(lattice, op, cell))
                .min()
                .unwrap_or(cell)
        })
        .collect())
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
        cell_orbits(&tracking, &materials.unwrap_or(default)).unwrap()
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
        assert!(cell_orbits(&tracking, &materials[1..]).is_err());
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
            let images = images(&tracking, &QUARTER).unwrap();
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
        let images = images(&tracking, &QUARTER).unwrap();
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
        let images = images(&tracking, &QUARTER).unwrap();
        let names: Vec<String> = (0..lattice.cells())
            .flat_map(|cell| {
                if cell == 4 {
                    vec!["water".to_string()]
                } else {
                    vec!["fuel".to_string(), "water".to_string()]
                }
            })
            .collect();
        assert!(check_regions(&images, &names).is_ok());
        // 右上の角のピンの燃料だけ別の名前
        let mut changed = names.clone();
        changed[lattice.cell_regions(2).start] = "other".into();
        let error = check_regions(&images, &changed).unwrap_err().to_string();
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
        let error = images(&tracking, &QUARTER).unwrap_err();
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
        assert!(images(&tracking, &QUARTER).is_err());
        // 左の真ん中だけ水
        let (asymmetric, _) = lattice(&[3]);
        let tracking = Tracking::new(&asymmetric, 8, 0.2).unwrap();
        assert!(images(&tracking, &QUARTER).is_err());
        // 左右が真空で上下が鏡面反射なら、1/8 対称でなくても 1/4 対称
        let boundaries = Boundaries {
            left: Boundary::Vacuum,
            right: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let tracking = Tracking::with_boundaries(&symmetric, 8, 0.2, boundaries).unwrap();
        assert!(images(&tracking, &QUARTER).is_ok());
    }
}
