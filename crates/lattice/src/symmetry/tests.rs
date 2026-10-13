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
fn none_is_chosen_without_detection() {
    // "none" は検証の基準なので、判定をせず（判定が止まる入力でも）、自己遮蔽も MOC も対称性なしで解く。値の数が
    // 合わない入力は、判定するときだけエラーになる
    let (tracking, values) = patterned([2; 8], Boundaries::REFLECTIVE);
    let short = &values[1..];
    let chosen = choose(&tracking, short, Some(Symmetry::NONE)).unwrap();
    assert_eq!(chosen.shielding().symmetry(), Symmetry::NONE);
    assert_eq!(chosen.moc().symmetry(), Symmetry::NONE);
    assert!(choose(&tracking, short, None).is_err());
    // 判定するときは、自己遮蔽に成り立つ操作をすべて、MOC に強制した対称性を渡す
    let chosen = choose(&tracking, &values, Some(Symmetry::QUARTER)).unwrap();
    assert_eq!(chosen.shielding().symmetry(), Symmetry::OCTANT);
    assert_eq!(chosen.moc().symmetry(), Symmetry::QUARTER);
    let chosen = choose(&tracking, &values, None).unwrap();
    assert_eq!(chosen.moc().symmetry(), Symmetry::OCTANT);
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
    // 成り立ったが絞り込みで外した操作を、もう一度求めると、外した理由を示す（「確かめていない」としない）
    let (tracking, values) = patterned([2; 8], Boundaries::REFLECTIVE);
    let quarter = detect(&tracking, &values)
        .unwrap()
        .restrict(Symmetry::QUARTER)
        .unwrap();
    let error = quarter.restrict(Symmetry::OCTANT).unwrap_err().to_string();
    assert!(error.contains("leaves it out"), "{error}");
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
        let apply = |m: [f64; 4], p: (f64, f64)| (m[0] * p.0 + m[1] * p.1, m[2] * p.0 + m[3] * p.1);
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
