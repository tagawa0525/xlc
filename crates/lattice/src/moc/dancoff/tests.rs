use super::*;
use crate::geometry::{Lattice, PinCell};

fn options() -> DancoffOptions {
    DancoffOptions::default()
}

/// ピンと水のセルを `layout`（0: ピン、1: 水）に並べ、燃料以外の群の全断面積を `total` とした Dancoff 係数
fn factors(pitch: f64, layout: Vec<Vec<usize>>, total: &[f64]) -> Vec<Option<Vec<f64>>> {
    let pin = PinCell::new(pitch, vec![0.4096]).unwrap();
    let water = PinCell::new(pitch, vec![]).unwrap();
    let lattice = Lattice::new(vec![pin, water], layout).unwrap();
    let tracking = Tracking::new(&lattice, 16, 0.02).unwrap();
    let regions = lattice.volumes().len();
    // ピンのセルは (燃料, 水) の2領域、水のセルは1領域
    let mut fuel = Vec::with_capacity(regions);
    for cell in 0..lattice.cells() {
        let range = lattice.cell_regions(cell);
        fuel.push(range.len() == 2);
        fuel.extend(std::iter::repeat_n(false, range.len() - 1));
    }
    let totals = vec![total.to_vec(); regions];
    dancoff(&tracking, &totals, &fuel, &options()).unwrap()
}

#[test]
fn distant_pins_do_not_shadow_each_other() {
    // 隣のピンまで減速材の光学的厚さが 20 程度あるので、影はほぼない
    let c = factors(20.0, vec![vec![0]], &[1.0, 2.0]);
    for value in c[0].as_ref().unwrap() {
        assert!(value.abs() < 1e-6, "{value}");
    }
}

#[test]
fn transparent_moderator_shadows_completely() {
    let c = factors(1.26, vec![vec![0]], &[1.0e-6]);
    let value = c[0].as_ref().unwrap()[0];
    assert!(value > 0.999, "{value}");
}

#[test]
fn repeated_pins_match_the_single_cell_lattice() {
    let total = [0.7, 1.5];
    let single = factors(1.26, vec![vec![0]], &total);
    let four = factors(1.26, vec![vec![0, 0], vec![0, 0]], &total);
    let reference = single[0].as_ref().unwrap();
    assert!(reference[0] > reference[1], "{reference:?}");
    assert!(
        reference.iter().all(|&c| c > 0.0 && c < 1.0),
        "{reference:?}"
    );
    for cell in &four {
        for (c, r) in cell.as_ref().unwrap().iter().zip(reference) {
            assert!((c - r).abs() < 2e-3, "{c} vs {r}");
        }
    }
}

#[test]
fn void_gap_matches_a_nearly_transparent_gap() {
    // 燃料・ギャップ・被覆管・水のピン
    let cell = PinCell::new(1.26, vec![0.4096, 0.418, 0.475]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 16, 0.02).unwrap();
    let fuel = [true, false, false, false];
    let factor = |gap: f64| {
        let totals = vec![vec![], vec![gap], vec![0.3], vec![1.2]];
        dancoff(&tracking, &totals, &fuel, &options()).unwrap()[0]
            .clone()
            .unwrap()[0]
    };
    let (void, thin) = (factor(0.0), factor(1.0e-9));
    assert!((void - thin).abs() < 1e-7, "{void} vs {thin}");
}

#[test]
fn fuel_rows_are_not_read() {
    // 燃料の行は使わないので、空でもよい
    let cell = PinCell::new(1.26, vec![0.4096]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.05).unwrap();
    let c = dancoff(&tracking, &[vec![], vec![1.2]], &[true, false], &options()).unwrap();
    assert_eq!(c[0].as_ref().map(Vec::len), Some(1));
}

#[test]
fn isolated_annulus_has_its_own_chord() {
    // 光学的に厚い水の中の孤立した円環（内径 0.5、外径 0.55）: C ≈ 0、l̄ = 4V/S = 2(r₂ − r₁)。
    // 水が薄いと、内側の水は中性子が入れず空になり、外側も遠くの円環の影を受ける
    let cell = PinCell::new(4.0, vec![0.5, 0.55]).unwrap();
    let lattice = Lattice::pin_cell(cell);
    let tracking = Tracking::new(&lattice, 32, 0.005).unwrap();
    let totals = vec![vec![5.0]; 3];
    let lumps = lump_dancoff(
        &tracking,
        &totals,
        &[false, true, false],
        &[None, Some(0), None],
        1,
        &options(),
    )
    .unwrap();
    assert!(lumps[0].dancoff[0] < 0.05, "{lumps:?}");
    assert!((lumps[0].chord - 0.1).abs() < 1e-3 * 0.1, "{lumps:?}");
}

#[test]
fn connected_frames_see_only_their_open_faces() {
    // 2×2 のセルの外周の枠（厚さ 0.05）は隣の枠とつながる。セルごとの塊の l̄ = 4V/S は、露出した内側の面
    // S = 4(p − 2t) だけで決まる
    let (p, t) = (1.26, 0.05);
    let cell = PinCell::new(p, vec![]).unwrap().with_frame(t).unwrap();
    let lattice = Lattice::new(vec![cell], vec![vec![0, 0], vec![0, 0]]).unwrap();
    let tracking = Tracking::new(&lattice, 32, 0.005).unwrap();
    let totals = vec![vec![0.3]; lattice.regions()];
    // セルごとに（水, 枠）の2領域
    let lump: Vec<Option<usize>> = (0..8).map(|r| (r % 2 == 1).then_some(r / 2)).collect();
    let black: Vec<bool> = lump.iter().map(Option::is_some).collect();
    let lumps = lump_dancoff(&tracking, &totals, &black, &lump, 4, &options()).unwrap();
    let expected = (p * p - (p - 2.0 * t).powi(2)) / (p - 2.0 * t);
    for l in &lumps {
        assert!((l.chord - expected).abs() < 1e-3 * expected, "{lumps:?}");
        assert!((0.0..=1.0).contains(&l.dancoff[0]), "{lumps:?}");
        assert!(
            (l.dancoff[0] - lumps[0].dancoff[0]).abs() < 1e-6,
            "{lumps:?}"
        );
    }
}

#[test]
fn black_regions_outside_lumps_shadow_them() {
    // 4つの枠をすべて黒体にし、1つだけを塊として数えると、4つとも塊にしたときと同じ結果
    let (p, t) = (1.26, 0.05);
    let cell = PinCell::new(p, vec![]).unwrap().with_frame(t).unwrap();
    let lattice = Lattice::new(vec![cell], vec![vec![0, 0], vec![0, 0]]).unwrap();
    let tracking = Tracking::new(&lattice, 16, 0.01).unwrap();
    let totals = vec![vec![0.3]; lattice.regions()];
    let black: Vec<bool> = (0..8).map(|r| r % 2 == 1).collect();
    let all: Vec<Option<usize>> = (0..8).map(|r| (r % 2 == 1).then_some(r / 2)).collect();
    let one: Vec<Option<usize>> = (0..8).map(|r| (r == 1).then_some(0)).collect();
    let expected = lump_dancoff(&tracking, &totals, &black, &all, 4, &options()).unwrap();
    let lumps = lump_dancoff(&tracking, &totals, &black, &one, 1, &options()).unwrap();
    assert_eq!(lumps.len(), 1);
    assert!(
        (lumps[0].chord - expected[0].chord).abs() < 1e-12,
        "{lumps:?}"
    );
    assert!(
        (lumps[0].dancoff[0] - expected[0].dancoff[0]).abs() < 1e-12,
        "{lumps:?} vs {expected:?}"
    );
}

#[test]
fn rejects_lumps_that_are_not_black() {
    let cell = PinCell::new(4.0, vec![0.5, 0.55]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.05).unwrap();
    let totals = vec![vec![5.0]; 3];
    let result = lump_dancoff(
        &tracking,
        &totals,
        &[false, false, false],
        &[None, Some(0), None],
        1,
        &options(),
    );
    assert!(
        matches!(result, Err(DancoffError::NotBlack(1))),
        "{result:?}"
    );
}

#[test]
fn rejects_lump_numbers_out_of_range() {
    let cell = PinCell::new(4.0, vec![0.5, 0.55]).unwrap();
    let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.05).unwrap();
    let totals = vec![vec![5.0]; 3];
    let result = lump_dancoff(
        &tracking,
        &totals,
        &[false, true, false],
        &[None, Some(1), None],
        1,
        &options(),
    );
    assert!(
        matches!(result, Err(DancoffError::Lump { lump: 1, lumps: 1 })),
        "{result:?}"
    );
}

#[test]
fn rejects_fuel_clipped_by_the_cell() {
    // 燃料の半径がセルの半ピッチより大きいと、隣のセルから直接燃料に入る
    let pin = PinCell::new(1.0, vec![0.6]).unwrap();
    let water = PinCell::new(1.0, vec![]).unwrap();
    let lattice = Lattice::new(vec![pin, water], vec![vec![0, 1], vec![1, 1]]).unwrap();
    let tracking = Tracking::new(&lattice, 8, 0.05).unwrap();
    let mut fuel = vec![false; lattice.volumes().len()];
    fuel[0] = true;
    let totals = vec![vec![1.0]; fuel.len()];
    let result = dancoff(&tracking, &totals, &fuel, &options());
    assert!(
        matches!(result, Err(DancoffError::FuelOnBoundary(0))),
        "{result:?}"
    );
}

#[test]
fn cells_without_fuel_have_no_factors() {
    let c = factors(1.26, vec![vec![0, 1], vec![1, 0]], &[1.0]);
    assert!(c[0].is_some() && c[3].is_some());
    assert!(c[1].is_none() && c[2].is_none());
}
