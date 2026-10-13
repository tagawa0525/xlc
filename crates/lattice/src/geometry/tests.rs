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
    let segments = cell.segments((-0.63, -0.3), d, 2.0, 0.0);
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

#[test]
fn lines_through_a_cell_corner_skip_its_neighbours() {
    // 3 × 3 の水のセル（ピッチ 1.26。境界の座標が2進で割り切れない）の内側の角を、いろいろな向きと始点で通す。
    // 角をちょうど通る直線は、角のまわりの4つのセルのうち、向かい合う2つだけを通る
    let water = PinCell::new(1.26, vec![]).unwrap();
    let lattice = Lattice::new(vec![water], vec![vec![0; 3]; 3]).unwrap();
    let cell = |r: usize| (r % 3, r / 3);
    for (cx, cy, column, row) in [(-0.63, 0.63, 1, 1), (0.63, -0.63, 2, 2), (0.63, 0.63, 2, 1)] {
        for angle in [0.1, 0.3, std::f64::consts::FRAC_PI_4, 1.1, 1.3] {
            let d = (f64::cos(angle), -f64::sin(angle));
            for back in [0.11, 0.37, 0.5, 0.61] {
                let start = (cx - back * d.0, cy - back * d.1);
                let s = lattice.segments(start, d, 2.0 * back);
                let cells: Vec<(usize, usize)> = s
                    .iter()
                    .map(|s| cell(lattice.region_cell(s.region)))
                    .collect();
                // 角の左上のセル（column - 1, row - 1）から右下のセル（column, row）へ直接移る
                assert_eq!(
                    cells,
                    vec![(column - 1, row - 1), (column, row)],
                    "corner ({cx}, {cy}), angle {angle}, back {back}: {s:?}"
                );
            }
        }
    }
}

/// ピッチ 1.26 cm の格子の同じ点とみなす距離（格子の大きさの 1e-9）
const COINCIDENT_1_26: f64 = COINCIDENT * 1.26;

/// 点 `point` をちょうど通るよう、そこから `back` だけ戻った点から、向き `angle` に `2 back` 進む直線を
/// `pin` で区切る（始点の座標の丸めで、通る点の近くの区切りがずれる）
fn through(pin: &PinCell, point: (f64, f64), angle: f64, back: f64) -> Vec<Segment> {
    let d = (angle.cos(), angle.sin());
    let start = (point.0 - back * d.0, point.1 - back * d.1);
    pin.segments(start, d, 2.0 * back, COINCIDENT_1_26)
}

#[test]
fn sector_rays_meet_once_at_the_centre() {
    // 扇形 8 のセルの中心を通る直線。8 本の扇形の境界は中心で交わり、丸めで少しずつずれても 1 点にまとまるので、
    // 中心のまわりに丸めの程度の区間を残さない
    let pin = PinCell::with_sectors(1.26, vec![0.3], 8).unwrap();
    for angle in [0.1, 0.3, 1.1, 2.0] {
        for back in [0.37, 0.41, 0.5] {
            let s = through(&pin, (0.0, 0.0), angle, back);
            assert!(s.iter().all(|s| s.length > 1e-12), "{angle} {back}: {s:?}");
        }
    }
}

#[test]
fn grid_corners_are_crossed_once() {
    // 3 × 3 の格子のセルの内側の角（±0.21, ±0.21）をちょうど通る直線は、角のまわりに丸めの程度の区間を残さない
    let pin = PinCell::grid(1.26, 3).unwrap();
    for angle in [0.1, 0.3, 1.1] {
        for back in [0.11, 0.37] {
            let s = through(&pin, (0.21, -0.21), angle, back);
            assert!(s.iter().all(|s| s.length > 1e-12), "{angle} {back}: {s:?}");
        }
    }
}

#[test]
fn frame_corners_are_crossed_once() {
    // 厚さ 0.05 cm の枠の付いたピンの、枠の内側の正方形の角（0.58, -0.58）をちょうど通る直線は、角のまわりに
    // 丸めの程度の区間を残さない
    let pin = PinCell::new(1.26, vec![0.4])
        .unwrap()
        .with_frame(0.05)
        .unwrap();
    for angle in [0.3, 0.7, 1.1] {
        for back in [0.03, 0.04] {
            let s = through(&pin, (0.58, -0.58), angle, back);
            assert!(s.iter().all(|s| s.length > 1e-12), "{angle} {back}: {s:?}");
        }
    }
}

#[test]
fn merged_cuts_do_not_depend_on_the_direction() {
    // 1 行 17 列のピン（ピッチ 1.26、半径 0.4096。同じ点とみなす距離は 2.1e-8 cm。行の中心はちょうど 0）の円を、
    // 中心からの距離が半径より 1 ulp だけ短い横の直線が横切る。円を横切る 2 点は 1.3e-8 cm しか離れておらず 1 点に
    // まとまる。逆向きに通っても、区間の長さは丸めの程度で同じ（まとめた点が進む向きによらない）
    let pin = PinCell::new(1.26, vec![0.4096]).unwrap();
    let lattice = Lattice::rectilinear(
        vec![1.26; 17],
        vec![1.26],
        vec![Cell::Pin(pin)],
        vec![vec![0; 17]],
    )
    .unwrap();
    let half = 0.5 * lattice.width();
    let y = 0.4096_f64.next_down();
    let forward = lattice.segments((-half, y), (1.0, 0.0), 2.0 * half);
    let backward = lattice.segments((half, y), (-1.0, 0.0), 2.0 * half);
    assert_eq!(forward.len(), backward.len(), "{forward:?} {backward:?}");
    for (a, b) in forward.iter().zip(backward.iter().rev()) {
        assert_eq!(a.region, b.region);
        assert!(
            (a.length - b.length).abs() < 1e-12,
            "{} vs {}",
            a.length,
            b.length
        );
    }
}

#[test]
fn regions_thinner_than_the_same_point_distance_are_errors() {
    // セルの中の区切りの幅（リングの幅、枠の厚さ、格子の幅、ブロックの小さな長方形）が同じ点とみなす
    // 距離（2 × 2、ピッチ 1.26 の格子で 2.5e-9 cm）以下なら、格子を作るときに止める（#0055）。斜めに浅く横切る
    // トラックでは区間が残り、体積が黙って過小になる
    let thin = |pin: PinCell| {
        matches!(
            Lattice::new(vec![pin], vec![vec![0; 2]; 2]),
            Err(GeometryError::ThinRegion { .. })
        )
    };
    assert!(thin(PinCell::new(1.26, vec![0.4, 0.4 + 1e-12]).unwrap()));
    assert!(thin(PinCell::new(1.26, vec![1e-12]).unwrap()));
    assert!(thin(
        PinCell::new(1.26, vec![0.4])
            .unwrap()
            .with_frame(1e-12)
            .unwrap()
    ));
    // 枠の内側の正方形に接する円は、間の四隅の領域が細い帯にならないので止めない
    assert!(!thin(
        PinCell::new(1.26, vec![0.58 - 1e-12])
            .unwrap()
            .with_frame(0.05)
            .unwrap()
    ));
    assert!(thin(PinCell::grid(1.26, 1_000_000_000).unwrap()));
    assert!(!thin(PinCell::new(1.26, vec![0.4, 0.4 + 1e-6]).unwrap()));
    let block = Lattice::rectilinear(
        vec![1.26; 2],
        vec![1.26; 2],
        vec![Cell::Block {
            nx: 1_000_000_000,
            ny: 1,
        }],
        vec![vec![0; 2]; 2],
    );
    assert!(matches!(block, Err(GeometryError::ThinRegion { .. })));
}

#[test]
fn cells_thinner_than_the_same_point_distance_are_errors() {
    // 幅が同じ点とみなす距離（格子の大きさの 1e-9）以下の列は、両側の線を横切る点が1つにまとまり、区間を持てない
    let pin = PinCell::new(1.26, vec![0.4]).unwrap();
    assert!(matches!(
        Lattice::with_margin(vec![pin.clone()], vec![vec![0; 2]; 2], 1e-12),
        Err(GeometryError::ThinCell { .. })
    ));
    assert!(Lattice::with_margin(vec![pin], vec![vec![0; 2]; 2], 1e-6).is_ok());
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
    let s = cell.segments((-1.0, -0.5), (1.0, 0.0), 2.0, 0.0);
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
    let s = cell.segments((-1.0, 0.25), (1.0, 0.0), 2.0, 0.0);
    let regions: Vec<usize> = s.iter().map(|s| s.region).collect();
    assert_eq!(regions, vec![5, 1, 0, 4]);
    let half_chord = (0.25_f64 - 0.0625).sqrt();
    let lengths: Vec<f64> = s.iter().map(|s| s.length).collect();
    for (a, b) in lengths
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
    let s = cell.segments((-0.6, -0.8), (0.6, 0.8), 2.0, 0.0);
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
