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

/// 物質ごとに 1 つの軸の配置
fn layout() -> Layout {
    Layout::new(&cell(), CLASSES.to_vec(), vec![0, 1, 2]).unwrap()
}

fn table(ranges: &[[f64; 2]], density: usize) -> TrackTable {
    TrackTable::new(&cell(), layout(), ranges, density).unwrap()
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

/// 内挿した飛程 T[j][i] と格子点ごとに求めた飛程 W[i][j] の差。真空でない列は衝突確率の差 |ΔW_ij| Σ_j / V_i、
/// 真空の列はその列の最大の飛程との比（黒い燃料の奥のリングからの飛程は 1e-20 の程度になり、相対差では比べられない）
fn difference(sigma: &[f64], got: &[Vec<f64>], exact: &[Vec<f64>]) -> f64 {
    let volumes = cell().volumes();
    let totals = totals(sigma);
    let mut worst: f64 = 0.0;
    for (i, v) in volumes.iter().enumerate() {
        for (j, &t) in totals.iter().enumerate() {
            let d = (got[j][i] - exact[i][j]).abs();
            let scale = if t > 0.0 {
                v / t
            } else {
                exact.iter().map(|row| row[j]).fold(0.0, f64::max)
            };
            worst = worst.max(d / scale);
        }
    }
    worst
}

/// 試料の全断面積の組での [`difference`] の最大
fn worst(table: &TrackTable) -> f64 {
    let cell = cell();
    samples()
        .iter()
        .map(|s| {
            difference(
                s,
                &table.tracks(s).unwrap().unwrap(),
                &cell.tracks(&totals(s)).unwrap(),
            )
        })
        .fold(0.0, f64::max)
}

#[test]
fn interpolated_tracks_converge_at_fourth_order() {
    // 4 点の Lagrange 補間なので、格子点の間隔を半分にすると差は約 1/16 になる
    let (coarse, fine) = (worst(&table(&RANGES, 8)), worst(&table(&RANGES, 16)));
    assert!(fine < 1e-5, "density 16: {fine}");
    assert!(coarse / fine > 10.0, "density 8: {coarse}, 16: {fine}");
}

#[test]
fn interpolated_tracks_keep_the_neutrons() {
    // 真空でない列の衝突の和は、どこかで必ず衝突するので体積になる（補間の重みの和は 1）
    let cell = cell();
    let volumes = cell.volumes();
    let table = table(&RANGES, DEFAULT_DENSITY);
    for s in samples() {
        let t = table.tracks(&s).unwrap().unwrap();
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
    let table = table(&RANGES, DEFAULT_DENSITY);
    for s in [[0.3, 0.25, 1.2], [3000.0, 2.0, 1.6]] {
        let d = difference(
            &s,
            &table.tracks(&s).unwrap().unwrap(),
            &cell.tracks(&totals(&s)).unwrap(),
        );
        assert!(d <= 1e-12, "{s:?}: {d}");
    }
}

#[test]
fn constant_materials_need_no_axis() {
    // 全断面積が変わらない物質は 1 点で済み、その点では格子点ごとに求めた値になる
    let cell = cell();
    let ranges = [[0.3, 3000.0], [0.9, 0.9], [1.4, 1.4]];
    let table = table(&ranges, DEFAULT_DENSITY);
    let s = [3000.0, 0.9, 1.4];
    let d = difference(
        &s,
        &table.tracks(&s).unwrap().unwrap(),
        &cell.tracks(&totals(&s)).unwrap(),
    );
    assert!(d <= 1e-12, "{d}");
}

#[test]
fn rejects_bad_ranges_and_densities() {
    let cell = cell();
    let new = |ranges: &[[f64; 2]], density| TrackTable::new(&cell, layout(), ranges, density);
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
    let short = Layout::new(&cell, CLASSES[..5].to_vec(), vec![0, 1, 2]);
    assert!(matches!(
        short,
        Err(TableError::Rings { given: 5, rings: 6 })
    ));
    // 物質に軸のないもの、物質のない軸は、配置にできない
    let unknown = Layout::new(&cell, CLASSES.to_vec(), vec![0, 1]);
    assert!(matches!(unknown, Err(TableError::Axes)));
    let gap = Layout::new(&cell, CLASSES.to_vec(), vec![0, 2, 2]);
    assert!(matches!(gap, Err(TableError::Axes)));
}

#[test]
fn tracks_reject_sigmas_off_the_table() {
    // 物質の数が違う全断面積と、表の範囲の外の全断面積は、内挿せずに止める
    let table = table(&RANGES, DEFAULT_DENSITY);
    assert!(matches!(
        table.tracks(&[1.0, 1.0]),
        Err(TableError::Materials {
            given: 2,
            materials: 3
        })
    ));
    assert!(matches!(
        table.tracks(&[3001.0, 1.0, 1.4]),
        Err(TableError::Outside { material: 0, .. })
    ));
}

#[test]
fn rejects_tables_too_large_to_count() {
    // 格子点の組の数が usize に収まらない細かさは止める
    let result = TrackTable::new(&cell(), layout(), &RANGES, usize::MAX / 4);
    assert!(
        matches!(result, Err(TableError::Size)),
        "{:?}",
        result.err()
    );
}

#[test]
fn node_count_is_known_before_building() {
    // 作る前に数えた格子点の組の数は、作った表の数と同じ
    for density in [4, 8, 16] {
        let table = table(&RANGES, density);
        assert_eq!(TrackTable::count(&RANGES, density).unwrap(), table.nodes());
    }
}

/// 燃料の 3 リングがそれぞれ別の物質（物質 0, 1, 2）で 1 つの軸を共有し、被覆管（物質 3）と水（物質 4）がそれぞれの軸を持つ配置
const RING_CLASSES: [Option<usize>; 6] = [Some(0), Some(1), Some(2), None, Some(3), Some(4)];
const RING_AXES: [usize; 5] = [0, 0, 0, 1, 2];

fn ring_layout() -> Layout {
    Layout::new(&cell(), RING_CLASSES.to_vec(), RING_AXES.to_vec()).unwrap()
}

fn ring_totals(sigma: &[f64; 5]) -> Vec<f64> {
    RING_CLASSES
        .iter()
        .map(|m| m.map_or(0.0, |m| sigma[m]))
        .collect()
}

/// 燃料のリングの ln Σ の、体積の重みの平均が ln `center` になる全断面積の組。リングの ln Σ の差は外側ほど大きく、内側から
/// −a、0、+a（重みの平均を引いて中心に合わせる）
fn ring_sigma(center: f64, a: f64, clad: f64, water: f64) -> [f64; 5] {
    let volumes = cell().volumes();
    let shape = [-a, 0.0, a];
    let mean: f64 = shape.iter().zip(&volumes).map(|(s, v)| s * v).sum::<f64>()
        / volumes[..3].iter().sum::<f64>();
    let fuel = shape.map(|s| center * (s - mean).exp());
    [fuel[0], fuel[1], fuel[2], clad, water]
}

/// [`ring_totals`] のセルでの [`difference`]
fn ring_difference(sigma: &[f64; 5], got: &[Vec<f64>]) -> f64 {
    let cell = cell();
    let totals = ring_totals(sigma);
    let exact = cell.tracks(&totals).unwrap();
    let volumes = cell.volumes();
    let mut worst: f64 = 0.0;
    for (i, v) in volumes.iter().enumerate() {
        for (j, &t) in totals.iter().enumerate() {
            let d = (got[j][i] - exact[i][j]).abs();
            let scale = if t > 0.0 {
                v / t
            } else {
                exact.iter().map(|row| row[j]).fold(0.0, f64::max)
            };
            worst = worst.max(d / scale);
        }
    }
    worst
}

#[test]
fn ring_materials_of_one_axis_stay_within_the_table_accuracy() {
    // ADR-0182: リングの物質の ln Σ の中心からの差が上限の中なら、1 次の展開の誤差は表の内挿（10 倍あたり 16 点で 1e-5 未満）の水準
    let layout = ring_layout();
    let ranges = [[0.3, 3000.0], [0.25, 2.0], [1.2, 1.6]];
    let table = TrackTable::new(&cell(), ring_layout(), &ranges, DEFAULT_DENSITY).unwrap();
    for [f, c, w] in samples() {
        let sigma = ring_sigma(f, MAX_DEVIATION / 2.0, c, w);
        let (_, deviation) = layout.centers(&sigma);
        assert!(deviation <= MAX_DEVIATION, "{sigma:?}: {deviation}");
        let got = table.tracks(&sigma).unwrap().unwrap();
        let d = ring_difference(&sigma, &got);
        assert!(d < 1e-5, "{sigma:?}: {d}");
    }
}

#[test]
fn ring_expansion_error_falls_as_the_square_of_the_deviation() {
    // 燃料の軸の範囲を狭くして内挿の誤差を消すと、残る 1 次の展開の誤差は、差を半分にすると約 1/4 になる
    for center in [0.31, 3.3, 47.0, 610.0] {
        let ranges = [[center, center * 1.001], [0.9, 0.9], [1.4, 1.4]];
        let table = TrackTable::new(&cell(), ring_layout(), &ranges, DEFAULT_DENSITY).unwrap();
        let error = |a: f64| {
            let sigma = ring_sigma(center * 1.0005, a, 0.9, 1.4);
            ring_difference(&sigma, &table.tracks(&sigma).unwrap().unwrap())
        };
        let (wide, narrow) = (error(MAX_DEVIATION / 2.0), error(MAX_DEVIATION / 4.0));
        assert!(wide < 1e-5, "{center}: {wide}");
        assert!(
            wide / narrow > 3.0 && wide / narrow < 5.0,
            "{center}: {wide} vs {narrow}"
        );
    }
}

#[test]
fn ring_deviations_beyond_the_limit_are_left_to_the_exact_tracks() {
    // 中心からの差が上限を超える組は、展開せずに格子点ごとに求めさせる
    let ranges = [[0.3, 3000.0], [0.25, 2.0], [1.2, 1.6]];
    let table = TrackTable::new(&cell(), ring_layout(), &ranges, DEFAULT_DENSITY).unwrap();
    let sigma = ring_sigma(3.3, 2.0 * MAX_DEVIATION, 0.9, 1.4);
    let (_, deviation) = ring_layout().centers(&sigma);
    assert!(deviation > MAX_DEVIATION, "{deviation}");
    assert!(table.tracks(&sigma).unwrap().is_none());
}

#[test]
fn expanded_tracks_keep_the_neutrons() {
    // 微分も、真空でない列の衝突の和（体積）を変えないので、展開した飛程も中性子を保つ
    let cell = cell();
    let volumes = cell.volumes();
    let ranges = [[0.3, 3000.0], [0.25, 2.0], [1.2, 1.6]];
    let table = TrackTable::new(&cell, ring_layout(), &ranges, DEFAULT_DENSITY).unwrap();
    for [f, c, w] in samples() {
        let sigma = ring_sigma(f, MAX_DEVIATION / 2.0, c, w);
        let t = table.tracks(&sigma).unwrap().unwrap();
        let totals = ring_totals(&sigma);
        for (i, v) in volumes.iter().enumerate() {
            let collided: f64 = (0..RADII.len()).map(|j| t[j][i] * totals[j]).sum();
            assert!(
                (collided - v).abs() <= 1e-12 * v,
                "{sigma:?} ring {i}: {collided} vs {v}"
            );
        }
    }
}

#[test]
fn shared_axes_cost_two_evaluations_per_material() {
    // 格子点の組ごとの衝突確率の数: 表の値 1 つと、軸を共有する物質ごとの中心差分の 2 つ
    assert_eq!(layout().evaluations(), 1);
    assert_eq!(ring_layout().evaluations(), 1 + 2 * 3);
}
