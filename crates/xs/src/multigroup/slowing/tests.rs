use super::*;
use crate::endf::Tab1;
use crate::multigroup::SlowingDown;

const RANGE: (f64, f64) = (1.0, 1.0e3);

fn constant(value: f64) -> Tab1 {
    Tab1::linear(vec![1.0e-5, 2.0e7], vec![value, value])
}

/// 範囲 RANGE の減速計算（ポテンシャル散乱 `potential` b）
fn slowing(awr: f64, potential: f64) -> SlowingDown {
    SlowingDown {
        awr,
        potential,
        range: RANGE,
    }
}

#[test]
fn pure_scatterer_keeps_the_weight_spectrum() {
    let (total, elastic) = (constant(10.0), constant(10.0));
    let shape = slowing_down_shape(
        &total,
        &elastic,
        &slowing(236.0, 10.0),
        &Weight::lwr(),
        50.0,
    )
    .unwrap();
    for y in &shape.y {
        assert!((y - 1.0).abs() < 1e-12, "{y}");
    }
}

#[test]
fn source_above_the_range_is_potential_scattering() {
    // 範囲の上の吸収体の散乱は、NR 近似どおりポテンシャル散乱 σp が重み関数のスペクトルで散乱するとする
    // （NJOY GROUPR の flux calculator と同じ）。範囲の上に散乱の共鳴があっても、範囲の中の純粋な
    // 散乱体（σs = σp）の中性子束は重み関数のまま。共鳴は上端のすぐ上（最初の窓の中）に置く
    let above = 1.005 * RANGE.1;
    let x = vec![1.0e-5, above, above, 2.0e7];
    let step = Tab1::linear(x, vec![10.0, 10.0, 30.0, 30.0]);
    let shape =
        slowing_down_shape(&step, &step, &slowing(236.0, 10.0), &Weight::lwr(), 50.0).unwrap();
    for y in &shape.y {
        assert!((y - 1.0).abs() < 1e-12, "{y}");
    }
}

#[test]
fn shape_ends_exactly_at_the_range() {
    // 端のエネルギーを計算し直すと丸めでずれ、範囲の端で表を引けなくなる
    let (total, elastic) = (constant(10.0), constant(8.0));
    let range = (0.1, 1.83e5);
    let slowing = SlowingDown {
        awr: 56.44629,
        potential: 8.0,
        range,
    };
    let shape = slowing_down_shape(&total, &elastic, &slowing, &Weight::lwr(), 50.0).unwrap();
    assert_eq!(shape.x.first(), Some(&range.0));
    assert_eq!(shape.x.last(), Some(&range.1));
}

#[test]
fn energies_with_equal_lethargy_give_finite_flux() {
    // 別の値なのに、レサジー ln(E_high / E) が丸めで同じになる2つのエネルギーを探す
    let lethargy = |e: f64| (RANGE.1 / e).ln();
    let mut e = 1.0_f64;
    while lethargy(e.next_up()) != lethargy(e) {
        e = e.next_up();
    }
    let x = vec![1.0e-5, e, e.next_up(), 2.0e7];
    let total = Tab1::linear(x.clone(), vec![10.0, 10.0, 12.0, 12.0]);
    let elastic = Tab1::linear(x, vec![8.0, 8.0, 9.0, 9.0]);
    let shape =
        slowing_down_shape(&total, &elastic, &slowing(236.0, 9.0), &Weight::lwr(), 50.0).unwrap();
    assert!(shape.y.iter().all(|y| y.is_finite()), "{:?}", shape.y);
}

#[test]
fn rejects_fewer_than_two_points_per_window() {
    let (total, elastic) = (constant(10.0), constant(8.0));
    for points in [0, 1] {
        let result = slowing_down_shape_with(
            &total,
            &elastic,
            &slowing(236.0, 8.0),
            &Weight::lwr(),
            50.0,
            points,
        );
        assert!(
            matches!(result, Err(SlowingDownError::WindowPoints(_))),
            "{points}"
        );
    }
}

#[test]
fn constant_absorber_approaches_the_asymptotic_flux() {
    // σa = 2, σs = 8, σ0 = 50: 上端から十分下では F = σ0 / (σa + σ0)
    let (total, elastic) = (constant(10.0), constant(8.0));
    let shape =
        slowing_down_shape(&total, &elastic, &slowing(236.0, 8.0), &Weight::lwr(), 50.0).unwrap();
    let low = shape.eval(2.0).unwrap();
    assert!((low - 50.0 / 52.0).abs() < 1e-10, "{low}");
}

/// 10 eV に共鳴（幅 0.5 eV、ピークで捕獲 2000 b・散乱 1000 b）、ほかは散乱 10 b
fn resonance() -> (Tab1, Tab1) {
    let x: Vec<f64> = (0..=20_000)
        .map(|i| 5.0 + 10.0 * f64::from(i) / 20_000.0)
        .collect();
    let peak = |e: f64| 1.0 / (1.0 + ((e - 10.0) / 0.25).powi(2));
    let elastic: Vec<f64> = x.iter().map(|&e| 10.0 + 1000.0 * peak(e)).collect();
    let total = x
        .iter()
        .zip(&elastic)
        .map(|(&e, s)| s + 2000.0 * peak(e))
        .collect();
    let wrap = |y: Vec<f64>| {
        // 範囲の外は一定値で延ばす
        let mut xs = vec![1.0e-5];
        xs.extend(&x);
        xs.push(2.0e7);
        let mut ys = vec![y[0]];
        ys.extend(&y);
        ys.push(y[y.len() - 1]);
        Tab1::linear(xs, ys)
    };
    (wrap(total), wrap(elastic))
}

#[test]
fn heavy_absorber_scatters_without_losing_energy() {
    // 質量が大きく、散乱でのエネルギーの損失（窓の幅 Δ ≈ 4/A = 4e-4）が共鳴の幅（5%）よりずっと
    // 小さいと、共鳴の中でも F ≈ σ0 / (σa + σ0)（散乱断面積によらない。NR 近似の 1/(σt + σ0) とは違う）
    let (total, elastic) = resonance();
    let shape = slowing_down_shape(
        &total,
        &elastic,
        &slowing(1.0e4, 10.0),
        &Weight::lwr(),
        50.0,
    )
    .unwrap();
    let at_peak = shape.eval(10.0).unwrap();
    assert!(
        (at_peak - 50.0 / 2050.0).abs() < 1e-3 * at_peak,
        "{at_peak}"
    );
}

#[test]
fn finer_grids_do_not_change_the_flux() {
    let (total, elastic) = resonance();
    let coarse = slowing_down_shape(
        &total,
        &elastic,
        &slowing(236.0, 10.0),
        &Weight::lwr(),
        50.0,
    )
    .unwrap();
    let fine = slowing_down_shape_with(
        &total,
        &elastic,
        &slowing(236.0, 10.0),
        &Weight::lwr(),
        50.0,
        4 * WINDOW_POINTS,
    )
    .unwrap();
    for e in [9.0, 9.9, 10.0, 10.1, 11.0, 3.0] {
        let (a, b) = (coarse.eval(e).unwrap(), fine.eval(e).unwrap());
        assert!((a - b).abs() < 1e-4 * b, "{e}: {a} vs {b}");
    }
}
