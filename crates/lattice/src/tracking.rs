//! 特性曲線（トラック）の生成。
//!
//! セルの外周は鏡面反射とする。反射したトラックが別のトラックにちょうどつながるように、
//! 方位角を tan φ = n_x / n_y（n_x, n_y は下辺・側辺から出るトラックの本数）に補正する
//! （cyclic tracking）。方位角は (0, π) の範囲で持ち、各トラックを順方向と逆方向の両方にスイープする。

use std::f64::consts::PI;

use crate::geometry::{Lattice, Segment};

/// トラックの端から出た角度中性子束の行き先（反射後に続くトラックと、そのスイープの向き）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Link {
    pub track: usize,
    pub forward: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Track {
    pub azimuth: usize,
    pub start: (f64, f64),
    pub end: (f64, f64),
    /// 始点と終点のある外周の辺（トラックは外周の角を通らない）
    pub start_side: Side,
    pub end_side: Side,
    /// 始点から終点への順
    pub segments: Vec<Segment>,
    /// 順方向のスイープで終点から出た角度中性子束の行き先（真空境界なら `None`）
    pub forward_link: Option<Link>,
    /// 逆方向のスイープで始点から出た角度中性子束の行き先（真空境界なら `None`）
    pub backward_link: Option<Link>,
}

/// 外周の辺。`usize` にすると 0 左、1 右、2 上、3 下。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
    Left,
    Right,
    Top,
    Bottom,
}

impl Side {
    /// 縦の辺（左右）か
    pub fn is_vertical(self) -> bool {
        matches!(self, Self::Left | Self::Right)
    }
}

impl Boundaries {
    /// 辺 `side` の境界条件
    pub fn of(&self, side: Side) -> Boundary {
        match side {
            Side::Left => self.left,
            Side::Right => self.right,
            Side::Top => self.top,
            Side::Bottom => self.bottom,
        }
    }
}

/// 1つの方位角。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Azimuth {
    /// 補正後の方位角 (0 < φ < π)
    pub phi: f64,
    /// トラックの間隔 (cm)
    pub spacing: f64,
    /// 求積の重み（(0, π) の全方位角で和が 1）
    pub weight: f64,
}

/// 外周の辺の境界条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Boundary {
    /// 鏡面反射
    Reflective,
    /// 真空（出た中性子は戻らない）
    Vacuum,
    /// 出た中性子のうち割合 α（0 ≤ α ≤ 1）が鏡面反射で戻る
    Albedo(f64),
}

/// 外周の4辺の境界条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Boundaries {
    pub left: Boundary,
    pub right: Boundary,
    pub top: Boundary,
    pub bottom: Boundary,
}

impl Boundaries {
    /// 4辺とも鏡面反射
    pub const REFLECTIVE: Self = Self {
        left: Boundary::Reflective,
        right: Boundary::Reflective,
        top: Boundary::Reflective,
        bottom: Boundary::Reflective,
    };
}

/// トラッキングの条件の誤り。
#[derive(Debug, thiserror::Error)]
pub enum TrackingError {
    #[error("the number of azimuthal angles must be positive and even, got {0}")]
    Azimuthal(usize),
    #[error("an albedo must be between 0 and 1, got {0}")]
    Albedo(f64),
    #[error("track spacing must be positive, got {0}")]
    Spacing(f64),
    #[error(
        "track spacing {spacing} is too coarse for {azimuthal} azimuthal angles \
         (corrected angles coincide)"
    )]
    CoarseSpacing { spacing: f64, azimuthal: usize },
    #[error("no reflected track at {0:?}")]
    NoReflectedTrack((f64, f64)),
    #[error("two tracks reflect into track {track} ({})", if *forward { "forward" } else { "backward" })]
    DuplicateLink { track: usize, forward: bool },
}

#[derive(Debug, Clone, PartialEq)]
pub struct Tracking {
    /// トラックを引いた形状（区間の領域の番号はこの形状の番号）
    pub lattice: Lattice,
    pub azimuths: Vec<Azimuth>,
    /// 方位角の順に並ぶ
    pub tracks: Vec<Track>,
    /// [方位角] の最初のトラックの番号（最後にトラックの数）
    first: Vec<usize>,
    /// [トラック][向き (0: 順, 1: 逆)] のスイープの始めに入ってくる角度中性子束を出すスイープ
    /// （反射のつながりの逆。真空境界なら `None`）
    sources: Vec<Option<Link>>,
    /// 外周の辺の境界条件
    boundaries: Boundaries,
}

impl Tracking {
    /// `n_azimuthal` は (0, π) の方位角の数（偶数）、`spacing` は希望するトラック間隔 (cm)。
    pub fn new(lattice: &Lattice, n_azimuthal: usize, spacing: f64) -> Result<Self, TrackingError> {
        Self::with_boundaries(lattice, n_azimuthal, spacing, Boundaries::REFLECTIVE)
    }

    /// 外周の4辺の境界条件を指定して作る（`new` は4辺とも鏡面反射）。
    pub fn with_boundaries(
        lattice: &Lattice,
        n_azimuthal: usize,
        spacing: f64,
        boundaries: Boundaries,
    ) -> Result<Self, TrackingError> {
        if n_azimuthal == 0 || !n_azimuthal.is_multiple_of(2) {
            return Err(TrackingError::Azimuthal(n_azimuthal));
        }
        if !(spacing > 0.0 && spacing.is_finite()) {
            return Err(TrackingError::Spacing(spacing));
        }
        for side in [Side::Left, Side::Right, Side::Top, Side::Bottom] {
            if let Boundary::Albedo(albedo) = boundaries.of(side)
                && !(0.0..=1.0).contains(&albedo)
            {
                return Err(TrackingError::Albedo(albedo));
            }
        }
        // α = 1 は鏡面反射、α = 0 は真空と同じなので、そちらに揃える（対称性の判定で同じものとして扱うため）
        let canonical = |b: Boundary| match b {
            Boundary::Albedo(1.0) => Boundary::Reflective,
            Boundary::Albedo(0.0) => Boundary::Vacuum,
            other => other,
        };
        let boundaries = Boundaries {
            left: canonical(boundaries.left),
            right: canonical(boundaries.right),
            top: canonical(boundaries.top),
            bottom: canonical(boundaries.bottom),
        };
        let n = n_azimuthal;
        let (width, height) = (lattice.width(), lattice.height());
        // 方位角を (0, π/2) で補正し、π - φ を対にする
        let mut azimuths = vec![
            Azimuth {
                phi: 0.0,
                spacing: 0.0,
                weight: 0.0,
            };
            n
        ];
        let mut counts = vec![(0, 0); n];
        for a in 0..n / 2 {
            let desired = PI / n as f64 * (a as f64 + 0.5);
            let nx = (width / spacing * desired.sin()).floor() as usize + 1;
            let ny = (height / spacing * desired.cos()).floor() as usize + 1;
            // 下辺に nx 本、側辺に ny 本が等間隔に入り、反射で端点がつながる角: tan φ = (H/ny)/(W/nx)
            let (dx, dy) = (width / nx as f64, height / ny as f64);
            let phi = dy.atan2(dx);
            // 補正後の間隔（下辺の間隔の、トラックに垂直な成分）は希望の間隔より小さい
            let delta = dx * phi.sin();
            azimuths[a] = Azimuth {
                phi,
                spacing: delta,
                weight: 0.0,
            };
            azimuths[n - 1 - a] = Azimuth {
                phi: PI - phi,
                spacing: delta,
                weight: 0.0,
            };
            counts[a] = (nx, ny);
            counts[n - 1 - a] = (nx, ny);
        }
        if azimuths.windows(2).any(|w| w[0].phi >= w[1].phi) {
            return Err(TrackingError::CoarseSpacing {
                spacing,
                azimuthal: n_azimuthal,
            });
        }
        // 重みは、隣の方位角との中点で区切った角度の幅
        for a in 0..n {
            let lower = if a == 0 {
                0.0
            } else {
                0.5 * (azimuths[a - 1].phi + azimuths[a].phi)
            };
            let upper = if a == n - 1 {
                PI
            } else {
                0.5 * (azimuths[a].phi + azimuths[a + 1].phi)
            };
            azimuths[a].weight = (upper - lower) / PI;
        }

        let (half, half_height) = (0.5 * width, 0.5 * height);
        let mut tracks = Vec::new();
        let mut first = vec![0; n + 1];
        for a in 0..n {
            first[a] = tracks.len();
            let (nx, ny) = counts[a];
            let (dx, dy) = (azimuths[a].phi.cos(), azimuths[a].phi.sin());
            let bottom =
                (0..nx).map(|i| (-half + width / nx as f64 * (i as f64 + 0.5), -half_height));
            // 右上向きなら左辺、左上向きなら右辺から入る
            let side_x = if dx > 0.0 { -half } else { half };
            let side =
                (0..ny).map(|j| (side_x, -half_height + height / ny as f64 * (j as f64 + 0.5)));
            let entry = if dx > 0.0 { Side::Left } else { Side::Right };
            let starts = bottom
                .map(|p| (p, Side::Bottom))
                .chain(side.map(|p| (p, entry)));
            for (start, start_side) in starts {
                let to_side = if dx > 0.0 {
                    (half - start.0) / dx
                } else {
                    (-half - start.0) / dx
                };
                let to_top = (half_height - start.1) / dy;
                let t = to_side.min(to_top);
                let end = (start.0 + t * dx, start.1 + t * dy);
                // 右上向きなら右辺、左上向きなら左辺か、上辺で出る
                let end_side = if to_side < to_top {
                    if dx > 0.0 { Side::Right } else { Side::Left }
                } else {
                    Side::Top
                };
                tracks.push(Track {
                    azimuth: a,
                    start,
                    end,
                    start_side,
                    end_side,
                    segments: lattice.segments(start, (dx, dy), t),
                    forward_link: None,
                    backward_link: None,
                });
            }
        }
        first[n] = tracks.len();

        // 反射でつながるトラックを、補角の方位角のトラックから端点の一致で探す
        let tolerance = lattice.coincident();
        let same = |p: (f64, f64), q: (f64, f64)| {
            (p.0 - q.0).abs() < tolerance && (p.1 - q.1).abs() < tolerance
        };
        for k in 0..tracks.len() {
            let a = tracks[k].azimuth;
            let complement = first[n - 1 - a]..first[n - a];
            let find = |point: (f64, f64), forward: bool| -> Result<Link, TrackingError> {
                complement
                    .clone()
                    .find(|&j| {
                        let candidate = if forward {
                            tracks[j].start
                        } else {
                            tracks[j].end
                        };
                        same(candidate, point)
                    })
                    .map(|track| Link { track, forward })
                    .ok_or(TrackingError::NoReflectedTrack(point))
            };
            // 側辺で反射すると補角の順方向、上辺・下辺で反射すると逆方向に続く
            let (end, end_side) = (tracks[k].end, tracks[k].end_side);
            let forward_link = match boundaries.of(end_side) {
                Boundary::Vacuum => None,
                Boundary::Reflective | Boundary::Albedo(_) => {
                    Some(find(end, end_side.is_vertical())?)
                }
            };
            let (start, start_side) = (tracks[k].start, tracks[k].start_side);
            let backward_link = match boundaries.of(start_side) {
                Boundary::Vacuum => None,
                Boundary::Reflective | Boundary::Albedo(_) => {
                    Some(find(start, start_side == Side::Bottom)?)
                }
            };
            tracks[k].forward_link = forward_link;
            tracks[k].backward_link = backward_link;
        }
        let mut sources = vec![None; 2 * tracks.len()];
        for (index, track) in tracks.iter().enumerate() {
            for (forward, link) in [(true, track.forward_link), (false, track.backward_link)] {
                if let Some(link) = link {
                    let target = &mut sources[2 * link.track + usize::from(!link.forward)];
                    if target.is_some() {
                        return Err(TrackingError::DuplicateLink {
                            track: link.track,
                            forward: link.forward,
                        });
                    }
                    *target = Some(Link {
                        track: index,
                        forward,
                    });
                }
            }
        }
        Ok(Self {
            lattice: lattice.clone(),
            azimuths,
            tracks,
            first,
            sources,
            boundaries,
        })
    }

    /// 外周の辺の境界条件（作るときに確かめたもの。あとから変えられない）
    pub fn boundaries(&self) -> Boundaries {
        self.boundaries
    }

    /// 方位角 `azimuth` のトラックの番号の範囲。
    pub fn azimuth_tracks(&self, azimuth: usize) -> std::ops::Range<usize> {
        self.first[azimuth]..self.first[azimuth + 1]
    }

    /// トラック `track` を向き `forward` にスイープするとき、始めに入ってくる角度中性子束を出すスイープ
    /// （真空境界から始まるなら `None`）。
    pub fn incoming_source(&self, track: usize, forward: bool) -> Option<Link> {
        self.sources[2 * track + usize::from(!forward)]
    }

    /// トラックから推定した領域ごとの体積。MOC の領域平均の計算と整合する。
    pub fn volumes(&self, regions: usize) -> Vec<f64> {
        let mut volumes = vec![0.0; regions];
        for track in &self.tracks {
            let azimuth = &self.azimuths[track.azimuth];
            for s in &track.segments {
                volumes[s.region] += azimuth.weight * azimuth.spacing * s.length;
            }
        }
        volumes
    }
}

#[cfg(test)]
mod tests;
