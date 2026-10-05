//! 特性曲線（トラック）の生成。
//!
//! セルの外周は鏡面反射とする。反射したトラックが別のトラックにちょうどつながるように、
//! 方位角を tan φ = n_x / n_y（n_x, n_y は下辺・側辺から出るトラックの本数）に補正する
//! （cyclic tracking）。方位角は (0, π) の範囲で持ち、各トラックを順方向と逆方向の両方にスイープする。

use crate::geometry::PinCell;

/// トラック上の1つの領域を通る区間。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Segment {
    pub region: usize,
    /// 平面に射影した長さ (cm)
    pub length: f64,
}

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
    /// 順方向のスイープで終点から出た角度中性子束の行き先
    pub forward_link: Link,
    /// 逆方向のスイープで始点から出た角度中性子束の行き先
    pub backward_link: Link,
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

#[derive(Debug, Clone, PartialEq)]
pub struct Tracking {
    pub azimuths: Vec<Azimuth>,
    pub tracks: Vec<Track>,
}

impl Tracking {
    /// `n_azimuthal` は (0, π) の方位角の数（偶数）、`spacing` は希望するトラック間隔 (cm)。
    pub fn new(cell: &PinCell, n_azimuthal: usize, spacing: f64) -> Result<Self, String> {
        todo!()
    }

    /// トラックから推定した領域ごとの体積。MOC の領域平均の計算と整合する。
    pub fn volumes(&self, regions: usize) -> Vec<f64> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cell() -> PinCell {
        PinCell::new(1.26, vec![0.2, 0.3, 0.4096, 0.55]).unwrap()
    }

    #[test]
    fn rejects_bad_settings() {
        assert!(Tracking::new(&cell(), 3, 0.05).is_err());
        assert!(Tracking::new(&cell(), 0, 0.05).is_err());
        assert!(Tracking::new(&cell(), 8, 0.0).is_err());
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
            let next = &t.tracks[track.forward_link.track];
            let joined = if track.forward_link.forward {
                next.start
            } else {
                next.end
            };
            assert!((joined.0 - track.end.0).abs() < 1e-9 && (joined.1 - track.end.1).abs() < 1e-9);
            let prev = &t.tracks[track.backward_link.track];
            let joined = if track.backward_link.forward {
                prev.start
            } else {
                prev.end
            };
            assert!(
                (joined.0 - track.start.0).abs() < 1e-9 && (joined.1 - track.start.1).abs() < 1e-9
            );
            // 反射先は補角の方位角
            assert_eq!(next.azimuth, t.azimuths.len() - 1 - track.azimuth);
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
            };
            if link == start {
                assert!(step > 1);
                return;
            }
        }
        panic!("the trajectory does not close");
    }
}
