use super::*;

use crate::geometry::{Cell, PinCell};

fn cell() -> Lattice {
    Lattice::pin_cell(PinCell::new(1.26, vec![0.2, 0.3, 0.4096, 0.55]).unwrap())
}

#[test]
fn nearly_axial_and_empty_lines_keep_their_cells() {
    // 2 × 2 の格子（ピッチ 2）。真下にほぼ平行な直線（x 成分は cos(π/2) の丸め）は、終わりの近くの横の線でも
    // セルを移り、長さ 0 の直線は区間を持たない
    let water = PinCell::new(2.0, vec![]).unwrap();
    let lattice = Lattice::new(vec![water], vec![vec![0; 2]; 2]).unwrap();
    let d = (std::f64::consts::FRAC_PI_2.cos(), -1.0);
    let s = lattice.segments((-1.0, 2.0), d, 3.0);
    let cells: Vec<usize> = s.iter().map(|s| lattice.region_cell(s.region)).collect();
    assert_eq!(cells, vec![0, 2], "{s:?}");
    assert!((s[0].length - 2.0).abs() < 1e-12 && (s[1].length - 1.0).abs() < 1e-12);
    assert!(lattice.segments((-1.0, 2.0), d, 0.0).is_empty());
}

#[test]
fn tracks_ending_on_a_ring_leave_no_round_off_segment() {
    // 2 × 2 の格子（ピッチ 1.26）の右上と左下が、半径 0.65 の円（扇形 8）が四隅に分かれたセル。方位角 16、間隔
    // 0.02 cm では、外周の辺と円の交わる点（局所 (0.16, 0.63)、0.16² + 0.63² = 0.65²）でちょうど終わるトラックが
    // ある。円を横切る点と終点が丸めでずれても、終点の手前に丸めの程度の区間を残さない（#0055）
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let corner = PinCell::with_sectors(1.26, vec![0.65], 8).unwrap();
    let lattice = Lattice::new(vec![pin, corner], vec![vec![0, 1], vec![1, 0]]).unwrap();
    let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
    for (index, track) in tracking.tracks.iter().enumerate() {
        let shortest = track
            .segments
            .iter()
            .map(|s| s.length)
            .fold(f64::INFINITY, f64::min);
        assert!(shortest > 1e-12, "track {index}: {:?}", track.segments);
    }
}

#[test]
fn tracks_leave_no_round_off_segments_in_blocks() {
    // 2 × 2 の、3 × 3 に分けたブロック（ピッチ 1.26）。45° を含む方位角では、ブロックの内側の格子の角を
    // ちょうど通るトラックがある。角のまわりの小さな長方形に丸めの程度の区間を残さない（#0055）
    let lattice = Lattice::rectilinear(
        vec![1.26; 2],
        vec![1.26; 2],
        vec![Cell::Block { nx: 3, ny: 3 }],
        vec![vec![0; 2]; 2],
    )
    .unwrap();
    for n_azimuthal in [6, 18, 30] {
        for spacing in [0.05, 0.03, 0.02] {
            let tracking = Tracking::new(&lattice, n_azimuthal, spacing).unwrap();
            for (index, track) in tracking.tracks.iter().enumerate() {
                let shortest = track
                    .segments
                    .iter()
                    .map(|s| s.length)
                    .fold(f64::INFINITY, f64::min);
                assert!(
                    shortest > 1e-12,
                    "azimuths {n_azimuthal}, spacing {spacing}, track {index}: {shortest}"
                );
            }
        }
    }
}

#[test]
fn tracks_leave_no_round_off_segments_at_cell_corners() {
    // 2 × 2、16 × 16、17 × 17 の均質なセル（区切りなし）の格子（余白なしとあり）。45° を含む方位角では、多くの
    // トラックがセルの角や、外周の辺とセルの境界の交わる点をちょうど通る。境界の線を足し合わせで作る丸めや、
    // 始点と終点の丸めがあっても、角のまわりのセルに丸めの程度の区間を残さない
    let water = PinCell::new(1.26, vec![]).unwrap();
    for n in [2usize, 16, 17] {
        for margin in [0.0, 0.04] {
            let lattice = if margin > 0.0 {
                Lattice::with_margin(vec![water.clone()], vec![vec![0; n]; n], margin).unwrap()
            } else {
                Lattice::new(vec![water.clone()], vec![vec![0; n]; n]).unwrap()
            };
            for n_azimuthal in [6, 18, 30] {
                for spacing in [0.05, 0.03] {
                    let tracking = Tracking::new(&lattice, n_azimuthal, spacing).unwrap();
                    for (index, track) in tracking.tracks.iter().enumerate() {
                        let shortest = track
                            .segments
                            .iter()
                            .map(|s| s.length)
                            .fold(f64::INFINITY, f64::min);
                        assert!(
                            shortest > 1e-12,
                            "n {n} margin {margin} azimuths {n_azimuthal}, spacing {spacing}, track {index}: {shortest}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn albedo_must_be_a_fraction() {
    for albedo in [-0.1, 1.1, f64::NAN] {
        let boundaries = Boundaries {
            top: Boundary::Albedo(albedo),
            ..Boundaries::REFLECTIVE
        };
        assert!(
            Tracking::with_boundaries(&cell(), 8, 0.05, boundaries).is_err(),
            "{albedo}"
        );
    }
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
        assert!((joined.0 - track.start.0).abs() < 1e-9 && (joined.1 - track.start.1).abs() < 1e-9);
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
