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
    /// 始点から終点への順
    pub segments: Vec<Segment>,
    /// 順方向のスイープで終点から出た角度中性子束の行き先（真空境界なら `None`）
    pub forward_link: Option<Link>,
    /// 逆方向のスイープで始点から出た角度中性子束の行き先（真空境界なら `None`）
    pub backward_link: Option<Link>,
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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// 鏡面反射
    Reflective,
    /// 真空（出た中性子は戻らない）
    Vacuum,
}

/// 外周の4辺の境界条件。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
    pub boundaries: Boundaries,
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
            for start in bottom.chain(side) {
                let to_side = if dx > 0.0 {
                    (half - start.0) / dx
                } else {
                    (-half - start.0) / dx
                };
                let to_top = (half_height - start.1) / dy;
                let t = to_side.min(to_top);
                let end = (start.0 + t * dx, start.1 + t * dy);
                tracks.push(Track {
                    azimuth: a,
                    start,
                    end,
                    segments: lattice.segments(start, (dx, dy), t),
                    forward_link: None,
                    backward_link: None,
                });
            }
        }
        first[n] = tracks.len();

        // 反射でつながるトラックを、補角の方位角のトラックから端点の一致で探す
        let tolerance = 1.0e-9 * width.max(height);
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
            // 点のある辺の境界条件（トラックは角を通らない）
            let boundary = |p: (f64, f64)| {
                if (p.0 + half).abs() < tolerance {
                    boundaries.left
                } else if (p.0 - half).abs() < tolerance {
                    boundaries.right
                } else if (p.1 - half_height).abs() < tolerance {
                    boundaries.top
                } else {
                    boundaries.bottom
                }
            };
            // 側辺で反射すると補角の順方向、上辺・下辺で反射すると逆方向に続く
            let end = tracks[k].end;
            let forward_link = match boundary(end) {
                Boundary::Vacuum => None,
                Boundary::Reflective => Some(find(end, (end.0.abs() - half).abs() < tolerance)?),
            };
            let start = tracks[k].start;
            let backward_link = match boundary(start) {
                Boundary::Vacuum => None,
                Boundary::Reflective => {
                    Some(find(start, (start.1 + half_height).abs() < tolerance)?)
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
mod tests {
    use super::*;

    use crate::geometry::{Cell, PinCell};

    fn cell() -> Lattice {
        Lattice::pin_cell(PinCell::new(1.26, vec![0.2, 0.3, 0.4096, 0.55]).unwrap())
    }

    #[test]
    fn rejects_bad_settings() {
        assert!(Tracking::new(&cell(), 3, 0.05).is_err());
        assert!(Tracking::new(&cell(), 0, 0.05).is_err());
        assert!(Tracking::new(&cell(), 8, 0.0).is_err());
    }

    #[test]
    fn coarse_spacing_that_merges_angles_is_error() {
        // 間隔が粗いと、補正後に隣の方位角と同じ角度になり、方向が実際より少なくなる
        assert!(Tracking::new(&cell(), 54, 0.6).is_err());
    }

    #[test]
    fn corrected_angles_are_symmetric() {
        let t = Tracking::new(&cell(), 16, 0.05).unwrap();
        assert_eq!(t.azimuths.len(), 16);
        for a in 0..16 {
            let (phi, complement) = (t.azimuths[a].phi, t.azimuths[15 - a].phi);
            assert!((phi + complement - std::f64::consts::PI).abs() < 1e-12);
            assert!(t.azimuths[a].spacing <= 0.05);
        }
        let total: f64 = t.azimuths.iter().map(|a| a.weight).sum();
        assert!((total - 1.0).abs() < 1e-12);
    }

    #[test]
    fn segments_cover_each_track() {
        let t = Tracking::new(&cell(), 8, 0.05).unwrap();
        for track in &t.tracks {
            let (dx, dy) = (track.end.0 - track.start.0, track.end.1 - track.start.1);
            let length: f64 = track.segments.iter().map(|s| s.length).sum();
            assert!((length - dx.hypot(dy)).abs() < 1e-12);
        }
    }

    #[test]
    fn track_volumes_match_analytic() {
        let cell = cell();
        let t = Tracking::new(&cell, 32, 0.005).unwrap();
        for (estimated, exact) in t.volumes(cell.regions()).iter().zip(cell.volumes()) {
            assert!(
                (estimated - exact).abs() < 1e-3 * exact,
                "{estimated} vs {exact}"
            );
        }
    }

    #[test]
    fn links_join_tracks_at_reflected_points() {
        let t = Tracking::new(&cell(), 8, 0.05).unwrap();
        for track in &t.tracks {
            let link = track.forward_link.unwrap();
            let next = &t.tracks[link.track];
            let joined = if link.forward { next.start } else { next.end };
            assert!((joined.0 - track.end.0).abs() < 1e-9 && (joined.1 - track.end.1).abs() < 1e-9);
            let link = track.backward_link.unwrap();
            let prev = &t.tracks[link.track];
            let joined = if link.forward { prev.start } else { prev.end };
            assert!(
                (joined.0 - track.start.0).abs() < 1e-9 && (joined.1 - track.start.1).abs() < 1e-9
            );
            // 反射先は補角の方位角
            assert_eq!(next.azimuth, t.azimuths.len() - 1 - track.azimuth);
        }
    }

    #[test]
    fn azimuth_ranges_hold_their_tracks() {
        let t = Tracking::new(&cell(), 8, 0.05).unwrap();
        let mut next = 0;
        for a in 0..t.azimuths.len() {
            let range = t.azimuth_tracks(a);
            assert_eq!(range.start, next);
            assert!(!range.is_empty());
            assert!(
                t.tracks[range.clone()]
                    .iter()
                    .all(|track| track.azimuth == a)
            );
            next = range.end;
        }
        assert_eq!(next, t.tracks.len());
    }

    #[test]
    fn incoming_sources_invert_the_links() {
        let boundaries = Boundaries {
            right: Boundary::Vacuum,
            bottom: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let t = Tracking::with_boundaries(&two_by_two(), 8, 0.05, boundaries).unwrap();
        let mut fed = 0;
        for (index, track) in t.tracks.iter().enumerate() {
            for (forward, link) in [(true, track.forward_link), (false, track.backward_link)] {
                if let Some(link) = link {
                    let source = Link {
                        track: index,
                        forward,
                    };
                    assert_eq!(t.incoming_source(link.track, link.forward), Some(source));
                    fed += 1;
                }
            }
        }
        // 真空境界から始まるスイープには出どころがない
        let empty = (0..t.tracks.len())
            .flat_map(|i| [(i, true), (i, false)])
            .filter(|&(i, forward)| t.incoming_source(i, forward).is_none())
            .count();
        assert_eq!(fed + empty, 2 * t.tracks.len());
        assert!(empty > 0);
    }

    /// ピン 2 本の右に、幅 0.5 と 1.5 の均質な列を足した 2 行 4 列の長方形の格子
    fn rectangle() -> Lattice {
        let pin = PinCell::new(1.26, vec![0.4]).unwrap();
        Lattice::rectilinear(
            vec![1.26, 1.26, 0.5, 1.5],
            vec![1.26, 1.26],
            vec![Cell::Pin(pin), Cell::BLOCK, Cell::Block { nx: 2, ny: 3 }],
            vec![vec![0, 0, 1, 2], vec![0, 0, 1, 2]],
        )
        .unwrap()
    }

    #[test]
    fn rectangle_track_volumes_match_analytic() {
        let lattice = rectangle();
        let t = Tracking::new(&lattice, 16, 0.002).unwrap();
        for (estimated, exact) in t.volumes(lattice.regions()).iter().zip(lattice.volumes()) {
            assert!(
                (estimated - exact).abs() < 1e-3 * exact,
                "{estimated} vs {exact}"
            );
        }
    }

    #[test]
    fn rectangle_tracks_close_by_reflection() {
        let t = Tracking::new(&rectangle(), 8, 0.05).unwrap();
        for first in [0, t.tracks.len() / 2] {
            let start = Link {
                track: first,
                forward: true,
            };
            let mut link = start;
            let closed = (1..=10 * t.tracks.len()).any(|_| {
                let track = &t.tracks[link.track];
                link = if link.forward {
                    track.forward_link
                } else {
                    track.backward_link
                }
                .unwrap();
                link == start
            });
            assert!(closed, "track {first}");
        }
    }

    #[test]
    fn following_links_returns_to_start() {
        // 鏡面反射の軌跡は有限回で元のトラックに戻る（cyclic）
        let t = Tracking::new(&cell(), 8, 0.05).unwrap();
        let start = Link {
            track: 0,
            forward: true,
        };
        let mut link = start;
        for step in 1..=10 * t.tracks.len() {
            let track = &t.tracks[link.track];
            link = if link.forward {
                track.forward_link
            } else {
                track.backward_link
            }
            .unwrap();
            if link == start {
                assert!(step > 1);
                return;
            }
        }
        panic!("the trajectory does not close");
    }

    fn two_by_two() -> Lattice {
        let pin = PinCell::new(1.26, vec![0.3, 0.4096]).unwrap();
        let water = PinCell::new(1.26, vec![]).unwrap();
        Lattice::new(vec![pin, water], vec![vec![0, 1], vec![0, 0]]).unwrap()
    }

    #[test]
    fn lattice_track_volumes_match_analytic() {
        let lattice = two_by_two();
        let t = Tracking::new(&lattice, 32, 0.005).unwrap();
        for (estimated, exact) in t.volumes(lattice.regions()).iter().zip(lattice.volumes()) {
            assert!(
                (estimated - exact).abs() < 1e-3 * exact,
                "{estimated} vs {exact}"
            );
        }
    }

    #[test]
    fn lattice_tracks_reflect_at_the_outer_boundary() {
        let t = Tracking::new(&two_by_two(), 8, 0.05).unwrap();
        let half = 1.26;
        for track in &t.tracks {
            for p in [track.start, track.end] {
                let on_edge = (p.0.abs() - half).abs() < 1e-9 || (p.1.abs() - half).abs() < 1e-9;
                assert!(on_edge, "{p:?}");
            }
            let link = track.forward_link.unwrap();
            let next = &t.tracks[link.track];
            let joined = if link.forward { next.start } else { next.end };
            assert!((joined.0 - track.end.0).abs() < 1e-9 && (joined.1 - track.end.1).abs() < 1e-9);
        }
    }

    #[test]
    fn vacuum_sides_have_no_links() {
        // C5G7 と同じく右と下が真空、左と上が鏡面反射
        let boundaries = Boundaries {
            right: Boundary::Vacuum,
            bottom: Boundary::Vacuum,
            ..Boundaries::REFLECTIVE
        };
        let t = Tracking::with_boundaries(&two_by_two(), 8, 0.05, boundaries).unwrap();
        let half = 1.26;
        let on_vacuum = |p: (f64, f64)| (p.0 - half).abs() < 1e-9 || (p.1 + half).abs() < 1e-9;
        let (mut vacuum, mut reflective) = (0, 0);
        for track in &t.tracks {
            for (point, link) in [
                (track.end, track.forward_link),
                (track.start, track.backward_link),
            ] {
                assert_eq!(link.is_none(), on_vacuum(point), "{point:?}");
                if link.is_none() {
                    vacuum += 1;
                } else {
                    reflective += 1;
                }
            }
        }
        assert!(vacuum > 0 && reflective > 0);
    }

    #[test]
    fn clipped_ring_track_volumes_match_analytic() {
        // 内接円（0.63）の外の 0.7 と 0.8 の円は、セルの外周で切り取られる
        let cell = PinCell::with_sectors(1.26, vec![0.54, 0.7, 0.8], 8).unwrap();
        let lattice = Lattice::pin_cell(cell);
        let t = Tracking::new(&lattice, 32, 0.002).unwrap();
        for (estimated, exact) in t.volumes(lattice.regions()).iter().zip(lattice.volumes()) {
            assert!(
                (estimated - exact).abs() < 1e-3 * exact,
                "{estimated} vs {exact}"
            );
        }
    }

    #[test]
    fn framed_cell_track_volumes_match_analytic() {
        let cell = PinCell::with_sectors(1.26, vec![0.3, 0.54, 0.58], 8)
            .unwrap()
            .with_frame(0.02)
            .unwrap();
        let lattice = Lattice::pin_cell(cell);
        let t = Tracking::new(&lattice, 32, 0.001).unwrap();
        for (estimated, exact) in t.volumes(lattice.regions()).iter().zip(lattice.volumes()) {
            assert!(
                (estimated - exact).abs() < 1e-3 * exact,
                "{estimated} vs {exact}"
            );
        }
    }

    #[test]
    fn sector_cell_track_volumes_match_analytic() {
        let cell = PinCell::with_sectors(1.26, vec![0.3, 0.54, 0.58], 8).unwrap();
        let lattice = Lattice::pin_cell(cell);
        // 最も薄いリング（幅 0.04 cm）を 40 本のトラックが横切る間隔
        let t = Tracking::new(&lattice, 32, 0.001).unwrap();
        for (estimated, exact) in t.volumes(lattice.regions()).iter().zip(lattice.volumes()) {
            assert!(
                (estimated - exact).abs() < 1e-3 * exact,
                "{estimated} vs {exact}"
            );
        }
    }
}
