use super::*;

/// 燃料 3 リング・真空のギャップ・被覆管・水のセル
const RADII: [f64; 6] = [0.2, 0.3, 0.4096, 0.418, 0.475, 0.711];
/// [リング] の物質: 燃料 0、被覆管 1、水 2
const CLASSES: [Option<usize>; 6] = [Some(0), Some(0), Some(0), None, Some(1), Some(2)];
/// 物質ごとの全断面積の範囲 (1/cm)。燃料は共鳴の山と谷、被覆管と水は狭い
const RANGES: [[f64; 2]; 3] = [[0.3, 3000.0], [0.25, 2.0], [1.2, 1.6]];

fn cell() -> Cylinder {
    Cylinder::new(RADII.to_vec()).unwrap()
}

fn totals(sigma: &[f64]) -> Vec<f64> {
    CLASSES
        .iter()
        .map(|m| m.map_or(0.0, |m| sigma[m]))
        .collect()
}

/// 表の範囲の中の、格子点から外れた全断面積の組
fn samples() -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    for f in [0.31, 0.7, 3.3, 47.0, 610.0, 2900.0] {
        for c in [0.26, 0.9, 1.97] {
            for w in [1.21, 1.43, 1.59] {
                out.push([f, c, w]);
            }
        }
    }
    out
}

/// 内挿した飛程と、格子点ごとに求めた飛程の相対差の最大（リングの組と試料について）
fn worst(table: &TrackTable) -> f64 {
    let cell = cell();
    let mut worst: f64 = 0.0;
    for s in samples() {
        let exact = cell.tracks(&totals(&s)).unwrap();
        let got = table.tracks(&s);
        for i in 0..RADII.len() {
            for j in 0..RADII.len() {
                let d = (got[j][i] - exact[i][j]).abs() / exact[i][j].abs();
                worst = worst.max(d);
            }
        }
    }
    worst
}

#[test]
fn interpolated_tracks_converge_at_fourth_order() {
    // 4 点の Lagrange 補間なので、格子点の間隔を半分にすると差は約 1/16 になる
    let table = |density| TrackTable::new(&cell(), CLASSES.to_vec(), &RANGES, density).unwrap();
    let (coarse, fine) = (worst(&table(8)), worst(&table(16)));
    assert!(fine < 1e-5, "density 16: {fine}");
    assert!(coarse / fine > 10.0, "density 8: {coarse}, 16: {fine}");
}

#[test]
fn interpolated_tracks_keep_the_neutrons() {
    // 真空でない列の衝突の和は、どこかで必ず衝突するので体積になる（補間の重みの和は 1）
    let cell = cell();
    let volumes = cell.volumes();
    let table = TrackTable::new(&cell, CLASSES.to_vec(), &RANGES, DEFAULT_DENSITY).unwrap();
    for s in samples() {
        let t = table.tracks(&s);
        let sigma = totals(&s);
        for (i, v) in volumes.iter().enumerate() {
            let collided: f64 = (0..RADII.len()).map(|j| t[j][i] * sigma[j]).sum();
            assert!(
                (collided - v).abs() <= 1e-12 * v,
                "{s:?} ring {i}: {collided} vs {v}"
            );
        }
    }
}

#[test]
fn tracks_at_the_nodes_are_the_cell_tracks() {
    // 格子点の上では、内挿の重みが 1 つの点だけ 1 になり、格子点ごとに求めた値そのものになる
    let cell = cell();
    let table = TrackTable::new(&cell, CLASSES.to_vec(), &RANGES, DEFAULT_DENSITY).unwrap();
    for s in [[0.3, 0.25, 1.2], [3000.0, 2.0, 1.6]] {
        let exact = cell.tracks(&totals(&s)).unwrap();
        let got = table.tracks(&s);
        for i in 0..RADII.len() {
            for j in 0..RADII.len() {
                let d = (got[j][i] - exact[i][j]).abs() / exact[i][j].abs();
                assert!(
                    d <= 1e-12,
                    "{s:?} ({i}, {j}): {} vs {}",
                    got[j][i],
                    exact[i][j]
                );
            }
        }
    }
}

#[test]
fn constant_materials_need_no_axis() {
    // 全断面積が変わらない物質は 1 点で済み、その点では格子点ごとに求めた値になる
    let cell = cell();
    let ranges = [[0.3, 3000.0], [0.9, 0.9], [1.4, 1.4]];
    let table = TrackTable::new(&cell, CLASSES.to_vec(), &ranges, DEFAULT_DENSITY).unwrap();
    let s = [3000.0, 0.9, 1.4];
    let exact = cell.tracks(&totals(&s)).unwrap();
    let got = table.tracks(&s);
    assert!((got[5][0] - exact[0][5]).abs() <= 1e-12 * exact[0][5]);
}

#[test]
fn rejects_bad_ranges_and_densities() {
    let cell = cell();
    let new =
        |ranges: &[[f64; 2]], density| TrackTable::new(&cell, CLASSES.to_vec(), ranges, density);
    assert!(matches!(new(&RANGES, 0), Err(TableError::Density)));
    let zero = [[0.0, 3000.0], [0.25, 2.0], [1.2, 1.6]];
    assert!(matches!(
        new(&zero, 16),
        Err(TableError::Range { material: 0, .. })
    ));
    let infinite = [[0.3, 3000.0], [0.25, f64::INFINITY], [1.2, 1.6]];
    assert!(matches!(
        new(&infinite, 16),
        Err(TableError::Range { material: 1, .. })
    ));
    let reversed = [[0.3, 3000.0], [0.25, 2.0], [1.6, 1.2]];
    assert!(matches!(
        new(&reversed, 16),
        Err(TableError::Range { material: 2, .. })
    ));
    let short = TrackTable::new(&cell, CLASSES[..5].to_vec(), &RANGES, 16);
    assert!(matches!(
        short,
        Err(TableError::Rings { given: 5, rings: 6 })
    ));
}
