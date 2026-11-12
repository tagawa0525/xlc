use super::*;
use xlc_format::fewgroup::Pins;
use xlc_format::table::{Transform, product_points};

/// 各量が (ppm, 減速材の温度, √燃料の温度) の多重線形な関数（交差の項を含む）である状態
fn value(ppm: f64, tm: f64, tf: f64) -> f64 {
    1.0 + 1.0e-4 * ppm + 2.0e-3 * tm + 3.0e-2 * tf.sqrt() + 1.0e-7 * ppm * tm
}

fn few(ppm: f64, tm: f64, tf: f64) -> FewGroup {
    let v = value(ppm, tm, tf);
    FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        k_inf: v,
        total: vec![0.5 * v, v],
        absorption: vec![0.01 * v, 0.08 * v],
        nu_fission: vec![0.006 * v, 0.1 * v],
        fission: vec![0.0025 * v, 0.04 * v],
        kappa_fission: None,
        chi: vec![1.0, 0.0],
        diffusion: vec![1.4 * v, 0.4 * v],
        scatter: vec![vec![0.47 * v, 0.02 * v], vec![0.0, 0.92 * v]],
        adf: vec![vec![v, 2.0 * v]; 6],
        vacuum_ratio: vec![
            None,
            Some(vec![0.5 * v / 3.0, 0.4 * v / 3.0]),
            None,
            None,
            None,
            None,
        ],
        normalized_adf: vec![false; 6],
        height: None,
        // 形状関数と角の不連続因子も内挿し、セルの境界はそのまま（ADR-0096）
        pins: Some(Pins {
            dx: vec![1.0, 1.0],
            dy: vec![2.0],
            form: vec![vec![vec![v, 2.0 * v]], vec![vec![3.0 * v, 0.0]]],
            corners: vec![vec![1.1 * v, 1.2 * v]; 4],
            normalized_corners: vec![false; 4],
            halves: vec![false; 4],
        }),
    }
}

const PPM: [f64; 3] = [900.0, 1300.0, 1700.0];
const TM: [f64; 3] = [550.0, 565.0, 580.0];
const TF: [f64; 2] = [550.0, 900.0];

fn table() -> Table {
    let mut states = Vec::new();
    for ppm in PPM {
        for tm in TM {
            for tf in TF {
                states.push(few(ppm, tm, tf));
            }
        }
    }
    let mut table = Table {
        points: Vec::new(),
        split: 1,
        heavy_metal: None,
        micro: None,
        histories: Vec::new(),
        kinetics: None,
        cooling: None,
        axes: vec![
            Axis::boron(PPM.to_vec(), 0.199),
            Axis::moderator(TM.to_vec(), vec![0.76972, 0.74271, 0.71190]),
            Axis::fuel(TF.to_vec()),
        ],
        states,
    };
    table.points = product_points(&table.axes);
    table
}

fn state(boron: f64, moderator: f64, fuel: f64) -> State {
    State {
        boron,
        b10: 0.199,
        moderator,
        fuel,
        burnup: None,
    }
}

fn close(a: f64, b: f64) -> bool {
    ((a - b) / b).abs() < 1e-13
}

#[test]
fn grid_points_return_the_lattice_states() {
    let t = table();
    let got = interpolate(&t, Some(&state(1300.0, 580.0, 900.0)), None)
        .unwrap()
        .remove(0);
    assert_eq!(got, few(1300.0, 580.0, 900.0));
}

#[test]
fn multilinear_quantities_are_reproduced_between_points() {
    let t = table();
    for (ppm, tm, tf) in [
        (1000.0, 560.0, 700.0),
        (1291.0, 570.0, 565.0),
        (1650.0, 551.0, 899.0),
    ] {
        let got = interpolate(&t, Some(&state(ppm, tm, tf)), None)
            .unwrap()
            .remove(0);
        let want = few(ppm, tm, tf);
        // すべての量を比べる
        let flat = |f: &FewGroup| {
            let mut v = vec![f.k_inf];
            for x in [
                &f.total,
                &f.absorption,
                &f.nu_fission,
                &f.fission,
                &f.chi,
                &f.diffusion,
            ] {
                v.extend(x);
            }
            v.extend(f.scatter.iter().flatten());
            v.extend(f.adf.iter().flatten());
            v.extend(f.vacuum_ratio.iter().flatten().flatten());
            let pins = f.pins.as_ref().unwrap();
            v.extend(pins.form.iter().flatten().flatten());
            v.extend(pins.corners.iter().flatten());
            v
        };
        let (g, w) = (flat(&got), flat(&want));
        assert_eq!(g.len(), w.len());
        for (g, w) in g.iter().zip(&w) {
            assert!(*w == 0.0 && *g == 0.0 || close(*g, *w), "{g} vs {w}");
        }
        assert_eq!(got.vacuum_ratio[0], None);
        assert_eq!(got.bounds, want.bounds);
        assert_eq!(got.normalized_adf, want.normalized_adf);
        let (gp, wp) = (got.pins.unwrap(), want.pins.unwrap());
        assert_eq!((gp.dx, gp.dy), (wp.dx, wp.dy));
    }
}

#[test]
fn fuel_temperature_is_interpolated_in_its_square_root() {
    // 燃料の温度だけを変えると、量は √T について線形に変わる
    let t = table();
    let got = interpolate(&t, Some(&state(1300.0, 565.0, 700.0)), None)
        .unwrap()
        .remove(0);
    let w = (700f64.sqrt() - 550f64.sqrt()) / (900f64.sqrt() - 550f64.sqrt());
    let want = (1.0 - w) * few(1300.0, 565.0, 550.0).k_inf + w * few(1300.0, 565.0, 900.0).k_inf;
    assert!(close(got.k_inf, want));
}

/// Problem 9 の燃料の温度の格子点（ADR-0239）
const TF_P9: [f64; 7] = [560.0, 600.0, 750.0, 900.0, 1050.0, 1200.0, 1350.0];

/// 燃料の温度の軸だけが点 `points` を持ち（ホウ素 1300 ppm、減速材 565 K の 1 点）、燃料の温度 T の状態が table() の状態の
/// `g(T)` 倍のテーブル
fn fuel_table_of(points: &[f64], g: impl Fn(f64) -> f64) -> Table {
    let mut table = Table {
        points: Vec::new(),
        split: 1,
        heavy_metal: None,
        micro: None,
        histories: Vec::new(),
        kinetics: None,
        cooling: None,
        axes: vec![
            Axis::boron(vec![1300.0], 0.199),
            Axis::moderator(vec![565.0], vec![0.74271]),
            Axis::fuel(points.to_vec()),
        ],
        states: points
            .iter()
            .map(|&tf| curved(g(tf), 1300.0, 565.0, 565.0))
            .collect(),
    };
    table.points = product_points(&table.axes);
    table
}

/// 燃料の温度の依存のモデル a + b√T + cT（高速群の共鳴吸収の √T と熱群の T。#0274）の倍率
fn doppler_model(tf: f64) -> f64 {
    0.9 + 4.0e-3 * tf.sqrt() + 2.0e-5 * tf
}

#[test]
fn fuel_temperatures_reproduce_quadratics_in_the_square_root_between_points() {
    // √T の 3 次 Hermite は √T の 2 次式を厳密に再現する（ADR-0242）。端の区間と内側の区間の点
    let t = fuel_table_of(&TF_P9, doppler_model);
    let nominal = few(1300.0, 565.0, 565.0).k_inf;
    for tf in [565.0, 580.0, 700.0, 820.0, 1111.0, 1300.0, 1349.0] {
        let got = interpolate(&t, Some(&state(1300.0, 565.0, tf)), None)
            .unwrap()
            .remove(0)
            .k_inf;
        let want = nominal * doppler_model(tf);
        assert!(close(got, want), "fuel {tf}: {got} vs {want}");
    }
}

#[test]
fn fuel_grid_points_return_the_lattice_states() {
    let t = fuel_table_of(&TF_P9, doppler_model);
    for (i, &tf) in TF_P9.iter().enumerate() {
        let got = interpolate(&t, Some(&state(1300.0, 565.0, tf)), None)
            .unwrap()
            .remove(0);
        assert_eq!(got, t.states[i], "fuel {tf}");
    }
}

/// 燃料の温度の点 `points` の内挿で、状態の燃料の温度 `tf` が点ごとの値に掛ける重み（点の状態だけを 1.1 倍にしたときの k∞ の増え方を
/// 0.1 で割る。2 倍にすると真空の境界の比が 1 を超えて検査に止められる）
fn fuel_weights(points: &[f64], tf: f64) -> Vec<f64> {
    let base = interpolate(
        &fuel_table_of(points, |_| 1.0),
        Some(&state(1300.0, 565.0, tf)),
        None,
    )
    .unwrap()
    .remove(0)
    .k_inf;
    (0..points.len())
        .map(|j| {
            let t = fuel_table_of(points, |x| if x == points[j] { 1.1 } else { 1.0 });
            (interpolate(&t, Some(&state(1300.0, 565.0, tf)), None)
                .unwrap()
                .remove(0)
                .k_inf
                / base
                - 1.0)
                / 0.1
        })
        .collect()
}

#[test]
fn fuel_weights_match_the_python_comparison_of_the_p9_grid() {
    // examples/depletion/compare_grid_interpolation.py の `hermite_weights` を √T の点に当てた重み。#0274 の 1 点を除いた誤差
    // （Slurm 752 の 2A のテーブルで 0.38 pcm。ADR-0242）は、7 点から内側の 1 点を除いた点で、除いた点の温度を引いたこの重み
    let without = |i: usize| -> Vec<f64> {
        TF_P9
            .iter()
            .enumerate()
            .filter(|&(k, _)| k != i)
            .map(|(_, &x)| x)
            .collect()
    };
    let cases = [
        (
            without(1),
            600.0,
            vec![
                (560.0, 0.6749953466231184),
                (750.0, 0.47001054765426514),
                (900.0, -0.1450058942773837),
            ],
        ),
        (
            without(2),
            750.0,
            vec![
                (560.0, -0.6818865068470399),
                (600.0, 1.1684138362216125),
                (900.0, 0.7222544137021589),
                (1050.0, -0.20878174307673159),
            ],
        ),
        (
            without(3),
            900.0,
            vec![
                (600.0, -0.131631116500441),
                (750.0, 0.5966216729620588),
                (1050.0, 0.7366594471155862),
                (1200.0, -0.20165000357720414),
            ],
        ),
        (
            without(4),
            1050.0,
            vec![
                (750.0, -0.1367118337025198),
                (900.0, 0.6067707898894349),
                (1200.0, 0.7265349651417771),
                (1350.0, -0.19659392132869194),
            ],
        ),
        (
            without(5),
            1200.0,
            vec![
                (900.0, -0.2900874895928406),
                (1050.0, 0.9351312343892599),
                (1350.0, 0.35495625520358065),
            ],
        ),
        (
            TF_P9.to_vec(),
            580.0,
            vec![
                (560.0, 0.43990153479800564),
                (600.0, 0.5761247225891927),
                (750.0, -0.01602625738719826),
            ],
        ),
        (
            TF_P9.to_vec(),
            1300.0,
            vec![
                (1050.0, -0.10008586475668797),
                (1200.0, 0.5335050628467114),
                (1350.0, 0.5665808019099766),
            ],
        ),
        (
            TF_P9.to_vec(),
            820.0,
            vec![
                (600.0, -0.05590689276888743),
                (750.0, 0.5775774841412364),
                (900.0, 0.5458990433575228),
                (1050.0, -0.06756963472987183),
            ],
        ),
    ];
    for (points, tf, want) in cases {
        let got = fuel_weights(&points, tf);
        for (j, &x) in points.iter().enumerate() {
            let w = want.iter().find(|&&(p, _)| p == x).map_or(0.0, |&(_, w)| w);
            assert!(
                (got[j] - w).abs() < 1e-13,
                "points {points:?} at {tf}: weight of {x} is {} vs {w}",
                got[j]
            );
        }
    }
}

#[test]
fn negative_fuel_interpolations_name_the_fuel_axis() {
    // 熱群の吸収が 900 K まで 0 で、そこから増える折れ線。900 K の傾きが正になり、750〜900 K の Hermite は 0 を下回る
    let mut t = fuel_table_of(&TF_P9, |_| 1.0);
    for (s, &tf) in t.states.iter_mut().zip(&TF_P9) {
        s.absorption[1] = 1.0e-5 * (tf - 900.0).max(0.0);
    }
    let err = interpolate(&t, Some(&state(1300.0, 565.0, 850.0)), None).unwrap_err();
    assert!(
        err.contains("negative") && err.contains("ADR-0242"),
        "{err}"
    );
    assert!(interpolate(&t, Some(&state(1300.0, 565.0, 1100.0)), None).is_ok());
}

#[test]
fn invalid_states_are_errors() {
    let t = table();
    let err = |s: Option<&State>, table: &Table| interpolate(table, s, None).unwrap_err();
    // 範囲の外（外挿しない）
    assert!(err(Some(&state(800.0, 565.0, 600.0)), &t).contains("boron"));
    assert!(err(Some(&state(1300.0, 590.0, 600.0)), &t).contains("moderator"));
    assert!(err(Some(&state(1300.0, 565.0, 1000.0)), &t).contains("fuel"));
    assert!(err(Some(&state(f64::NAN, 565.0, 600.0)), &t).contains("boron"));
    // B-10 の割合が範囲の外（ADR-0237）
    let mut s = state(1300.0, 565.0, 600.0);
    s.b10 = 1.5;
    assert!(err(Some(&s), &t).contains("B-10"));
    // 状態があるのに軸がない、状態がないのに軸がある
    let mut two = t.clone();
    two.axes.pop();
    two.points = product_points(&two.axes);
    two.states = two.states.iter().step_by(2).cloned().collect();
    assert!(err(Some(&state(1300.0, 565.0, 600.0)), &two).contains("fuel"));
    // 軸がないことは、ほかの軸の範囲の外より先に示す
    assert!(err(Some(&state(800.0, 565.0, 600.0)), &two).contains("fuel"));
    assert!(err(None, &t).contains("state"));
}

#[test]
fn single_point_axes_need_that_value() {
    let mut t = table();
    t.axes[2] = Axis::fuel(vec![550.0]);
    t.points = product_points(&t.axes);
    t.states = t.states.iter().step_by(2).cloned().collect();
    let got = interpolate(&t, Some(&state(1300.0, 565.0, 550.0)), None)
        .unwrap()
        .remove(0);
    assert_eq!(got, few(1300.0, 565.0, 550.0));
    assert!(interpolate(&t, Some(&state(1300.0, 565.0, 551.0)), None).is_err());
}

#[test]
fn tables_without_axes_need_no_state() {
    let single = Table::single(few(1300.0, 565.0, 565.0));
    assert_eq!(
        interpolate(&single, None, None).unwrap()[0],
        few(1300.0, 565.0, 565.0)
    );
    assert!(interpolate(&single, Some(&state(1300.0, 565.0, 565.0)), None).is_err());
}

#[test]
fn states_differing_in_pins_are_errors() {
    // テーブルの検査を経ずに作ったテーブルでも、形状関数の有無やセルが状態で違えば黙って内挿しない
    let mut t = table();
    t.states[1].pins = None;
    let err = interpolate(&t, Some(&state(900.0, 550.0, 700.0)), None).unwrap_err();
    assert!(err.contains("pins"), "{err}");
    let mut t = table();
    if let Some(p) = t.states[1].pins.as_mut() {
        p.dx[0] = 1.5;
    }
    let err = interpolate(&t, Some(&state(900.0, 550.0, 700.0)), None).unwrap_err();
    assert!(err.contains("pins"), "{err}");
}

#[test]
fn split_tables_interpolate_each_sub_node() {
    // 2 × 2 の部分ノード（ADR-0108）。部分ノード q の値は few(...) の (1 + q/10) 倍で、部分ノードごとに内挿する
    let mut t = table();
    let points = t.states.len();
    t.split = 2;
    t.states = (0..points)
        .flat_map(|p| {
            let base = t.states[p].clone();
            (0..4).map(move |q| {
                let mut s = base.clone();
                s.k_inf *= 1.0 + f64::from(u8::try_from(q).unwrap()) / 10.0;
                s
            })
        })
        .collect();
    let got = interpolate(&t, Some(&state(1291.0, 570.0, 565.0)), None).unwrap();
    assert_eq!(got.len(), 4);
    let want = few(1291.0, 570.0, 565.0).k_inf;
    for (q, sub) in got.iter().enumerate() {
        let expected = want * (1.0 + q as f64 / 10.0);
        assert!(
            close(sub.k_inf, expected),
            "sub-node {q}: {} vs {expected}",
            sub.k_inf
        );
    }
}

/// table() の軸の間（ホウ素の次）に制御棒の軸（なし、AIC、B4C）を挟んだテーブル。棒の値ごとに量を
/// 1、0.5、0.25 倍にする
fn rodded_table() -> Table {
    let factors = [1.0, 0.5, 0.25];
    let mut states = Vec::new();
    for ppm in PPM {
        for factor in factors {
            for tm in TM {
                for tf in TF {
                    let mut f = few(ppm, tm, tf);
                    f.absorption.iter_mut().for_each(|a| *a *= factor);
                    states.push(f);
                }
            }
        }
    }
    let mut t = table();
    t.axes.insert(
        1,
        Axis::Rod {
            values: vec!["out".into(), "aic".into(), "b4c".into()],
        },
    );
    t.points = product_points(&t.axes);
    t.states = states;
    t
}

#[test]
fn rod_axes_select_the_named_value_and_out_without_a_rod() {
    let t = rodded_table();
    let s = state(1100.0, 560.0, 700.0);
    let plain = interpolate(&table(), Some(&s), None).unwrap();
    let out = interpolate(&t, Some(&s), None).unwrap();
    let aic = interpolate(&t, Some(&s), Some("aic")).unwrap();
    let b4c = interpolate(&t, Some(&s), Some("b4c")).unwrap();
    for g in 0..2 {
        let a = plain[0].absorption[g];
        assert!(close(out[0].absorption[g], a));
        assert!(close(aic[0].absorption[g], 0.5 * a));
        assert!(close(b4c[0].absorption[g], 0.25 * a));
    }
    // 制御棒の軸のほかの量は同じ内挿
    assert!(close(b4c[0].total[1], plain[0].total[1]));
}

#[test]
fn rod_values_the_table_cannot_select_are_errors() {
    let t = rodded_table();
    let s = state(1100.0, 560.0, 700.0);
    // 知らない値
    assert!(interpolate(&t, Some(&s), Some("hafnium")).is_err());
    // 制御棒の軸のないテーブルの層に棒がある（黙って棒のない値にしない）
    assert!(interpolate(&table(), Some(&s), Some("aic")).is_err());
    // 棒のない状態（out）を持たない軸
    let mut always = rodded_table();
    always.axes[1] = Axis::Rod {
        values: vec!["aic".into(), "b4c".into(), "plenum".into()],
    };
    assert!(interpolate(&always, Some(&s), None).is_err());
    assert!(interpolate(&always, Some(&s), Some("plenum")).is_ok());
    // 軸のないテーブルは棒を選べない
    let single = Table::single(few(1300.0, 565.0, 565.0));
    assert!(interpolate(&single, None, Some("aic")).is_err());
}

const BURNUP: [f64; 3] = [0.0, 10.0, 20.0];

/// table() の各量を (1 + bu / 40) 倍した状態（燃焼度についても多重線形）
fn burnt(bu: f64, ppm: f64, tm: f64, tf: f64) -> FewGroup {
    let f = few(ppm, tm, tf);
    let mut out = scaled(&f, 1.0 + bu / 40.0);
    // 群境界と比で決まる量は倍にしない（chi は和が 1）
    out.chi = f.chi;
    out
}

/// table() の前に燃焼度の軸を置いたテーブル（燃焼度がいちばん遅く変わる。ADR-0187）
fn burnup_table() -> Table {
    let mut states = Vec::new();
    for bu in BURNUP {
        for ppm in PPM {
            for tm in TM {
                for tf in TF {
                    states.push(burnt(bu, ppm, tm, tf));
                }
            }
        }
    }
    let mut t = table();
    t.axes.insert(
        0,
        Axis::Burnup {
            points: BURNUP.to_vec(),
        },
    );
    t.points = product_points(&t.axes);
    t.states = states;
    t
}

fn at_burnup(burnup: f64, boron: f64, moderator: f64, fuel: f64) -> State {
    State {
        burnup: Some(burnup),
        ..state(boron, moderator, fuel)
    }
}

#[test]
fn burnup_grid_points_return_the_lattice_states() {
    let t = burnup_table();
    for bu in BURNUP {
        let got = interpolate(&t, Some(&at_burnup(bu, 1300.0, 580.0, 900.0)), None)
            .unwrap()
            .remove(0);
        assert_eq!(got, burnt(bu, 1300.0, 580.0, 900.0));
    }
}

#[test]
fn burnup_is_interpolated_linearly_with_the_other_axes() {
    let t = burnup_table();
    for (bu, ppm, tm, tf) in [(4.0, 1300.0, 565.0, 900.0), (15.5, 1000.0, 560.0, 700.0)] {
        let got = interpolate(&t, Some(&at_burnup(bu, ppm, tm, tf)), None)
            .unwrap()
            .remove(0);
        let want = burnt(bu, ppm, tm, tf);
        assert!(
            close(got.k_inf, want.k_inf),
            "{} vs {}",
            got.k_inf,
            want.k_inf
        );
        assert!(close(got.absorption[1], want.absorption[1]));
        assert!(close(got.chi[0], 1.0));
    }
}

#[test]
fn burnups_outside_the_points_or_missing_are_errors() {
    let t = burnup_table();
    let err = |s: &State| interpolate(&t, Some(s), None).unwrap_err();
    assert!(err(&at_burnup(20.5, 1300.0, 565.0, 900.0)).contains("burnup"));
    assert!(err(&at_burnup(-1.0, 1300.0, 565.0, 900.0)).contains("burnup"));
    assert!(err(&at_burnup(f64::NAN, 1300.0, 565.0, 900.0)).contains("burnup"));
    // 燃焼度の軸を持つテーブルに、燃焼度のない状態
    assert!(err(&state(1300.0, 565.0, 900.0)).contains("burnup"));
}

#[test]
fn fissile_tables_without_a_burnup_axis_reject_burnups() {
    // 燃料（核分裂のある）のテーブルを、燃焼度の軸なしに燃えた状態で使わない
    let err =
        interpolate(&table(), Some(&at_burnup(10.0, 1300.0, 565.0, 900.0)), None).unwrap_err();
    assert!(err.contains("burnup"), "{err}");
    // 燃焼度 0 は新しい燃料のテーブルの状態そのもの
    let fresh = interpolate(&table(), Some(&state(1300.0, 565.0, 900.0)), None).unwrap();
    let zero = interpolate(&table(), Some(&at_burnup(0.0, 1300.0, 565.0, 900.0)), None).unwrap();
    assert_eq!(fresh, zero);
    // 核分裂のないテーブル（反射体）は燃焼しないので、どの燃焼度でも同じ値
    let mut reflector = table();
    for s in &mut reflector.states {
        s.nu_fission = vec![0.0; 2];
        s.fission = vec![0.0; 2];
    }
    let fresh = interpolate(&reflector, Some(&state(1300.0, 565.0, 900.0)), None).unwrap();
    let burnt = interpolate(
        &reflector,
        Some(&at_burnup(10.0, 1300.0, 565.0, 900.0)),
        None,
    )
    .unwrap();
    assert_eq!(fresh, burnt);
}

/// table() の各量を `g(燃焼度)` 倍した状態を、燃焼度の点 `points` で並べたテーブル（chi は倍にしない）
fn burnup_table_of(points: &[f64], g: impl Fn(f64) -> f64) -> Table {
    let mut states = Vec::new();
    for &bu in points {
        for ppm in PPM {
            for tm in TM {
                for tf in TF {
                    states.push(curved(g(bu), ppm, tm, tf));
                }
            }
        }
    }
    let mut t = table();
    t.axes.insert(
        0,
        Axis::Burnup {
            points: points.to_vec(),
        },
    );
    t.points = product_points(&t.axes);
    t.states = states;
    t
}

/// table() の状態を `g` 倍したもの（chi は和が 1 なので倍にしない）
fn curved(g: f64, ppm: f64, tm: f64, tf: f64) -> FewGroup {
    let f = few(ppm, tm, tf);
    let mut out = scaled(&f, g);
    out.chi = f.chi;
    out
}

/// 間隔の揃わない燃焼度の点。0 から最初の点の間は、Xe と Sm の平衡までの急な変化を模す
const CURVED: [f64; 6] = [0.0, 0.5, 2.0, 5.0, 10.0, 18.0];

/// 燃焼度について 2 次の倍率。燃焼度 0 の点だけ曲線の外（Xe と Sm のない新しい燃料）
fn curve(bu: f64) -> f64 {
    if bu == 0.0 {
        1.3
    } else {
        1.0 + 0.03 * bu - 8.0e-4 * bu * bu
    }
}

#[test]
fn curved_burnup_grid_points_return_the_lattice_states() {
    let t = burnup_table_of(&CURVED, curve);
    for bu in CURVED {
        let got = interpolate(&t, Some(&at_burnup(bu, 1300.0, 580.0, 900.0)), None)
            .unwrap()
            .remove(0);
        assert_eq!(got, curved(curve(bu), 1300.0, 580.0, 900.0), "burnup {bu}");
    }
}

#[test]
fn burnups_past_the_first_point_reproduce_quadratics_with_the_other_axes() {
    let t = burnup_table_of(&CURVED, curve);
    for bu in [0.5, 0.7, 1.3, 3.1, 7.5, 12.0, 17.9] {
        for (ppm, tm, tf) in [(1300.0, 565.0, 900.0), (1000.0, 560.0, 700.0)] {
            let got = interpolate(&t, Some(&at_burnup(bu, ppm, tm, tf)), None)
                .unwrap()
                .remove(0);
            let want = curved(curve(bu), ppm, tm, tf);
            assert!(
                close(got.k_inf, want.k_inf),
                "burnup {bu}: {} vs {}",
                got.k_inf,
                want.k_inf
            );
            assert!(close(got.absorption[1], want.absorption[1]));
            assert!(close(got.scatter[0][1], want.scatter[0][1]));
            let form = |f: &FewGroup| f.pins.as_ref().unwrap().form[1][0][0];
            assert!(close(form(&got), form(&want)));
            assert!(close(got.chi[0], 1.0));
        }
    }
}

#[test]
fn burnups_below_the_first_point_are_linear_from_zero_and_zero_does_not_shape_the_rest() {
    let t = burnup_table_of(&CURVED, curve);
    let at = |t: &Table, bu: f64| {
        interpolate(t, Some(&at_burnup(bu, 1300.0, 565.0, 900.0)), None)
            .unwrap()
            .remove(0)
    };
    let fresh = few(1300.0, 565.0, 900.0).k_inf;
    let w = 0.2 / 0.5;
    let want = ((1.0 - w) * curve(0.0) + w * curve(0.5)) * fresh;
    assert!(close(at(&t, 0.2).k_inf, want));
    // 燃焼度 0 の値を変えても、最初の点から先は変わらない
    let other = burnup_table_of(&CURVED, |bu| if bu == 0.0 { 2.0 } else { curve(bu) });
    for bu in [0.5, 0.6, 1.0, 2.0, 4.0, 18.0] {
        assert_eq!(at(&t, bu), at(&other, bu), "burnup {bu}");
    }
}

#[test]
fn two_point_burnup_spans_are_linear() {
    // 0 のない 2 点の軸と、0 のあとに 2 点の軸（最初の点から先が 2 点）
    for points in [vec![5.0, 15.0], vec![0.0, 10.0, 20.0]] {
        let t = burnup_table_of(&points, |bu| 1.0 + 0.02 * bu + 1.0e-3 * bu * bu);
        let last = points.len() - 1;
        let (a, b) = (points[last - 1], points[last]);
        let g = |bu: f64| 1.0 + 0.02 * bu + 1.0e-3 * bu * bu;
        for bu in [a + 0.25 * (b - a), a + 0.5 * (b - a)] {
            let got = interpolate(&t, Some(&at_burnup(bu, 1300.0, 565.0, 900.0)), None)
                .unwrap()
                .remove(0);
            let w = (bu - a) / (b - a);
            let want = ((1.0 - w) * g(a) + w * g(b)) * few(1300.0, 565.0, 900.0).k_inf;
            assert!(close(got.k_inf, want), "{points:?} at {bu}");
        }
    }
}

#[test]
fn knees_overshoot_only_next_to_them_by_two_27ths_of_the_slope_change() {
    // 10 MWd/kg で傾きが 0 から 0.05 に変わる折れ線（Gd の燃え尽きの折れ曲がりを模す）。間隔 h = 2
    let points: Vec<f64> = (1..=10).map(|i| 2.0 * f64::from(i)).collect();
    let g = |bu: f64| 1.0 + 0.05 * (bu - 10.0).max(0.0);
    let t = burnup_table_of(&points, g);
    let fresh = few(1300.0, 565.0, 900.0).k_inf;
    let ratio = |bu: f64| {
        interpolate(&t, Some(&at_burnup(bu, 1300.0, 565.0, 900.0)), None)
            .unwrap()
            .remove(0)
            .k_inf
            / fresh
    };
    let samples = |a: f64, b: f64| (0..=200).map(move |i| a + (b - a) * f64::from(i) / 200.0);
    // 折れ曲がりから 2 区間より離れたところは、折れ線そのもの（両側の 3 点が同じ直線に乗る）
    for bu in samples(2.0, 8.0).chain(samples(12.0, 20.0)) {
        assert!((ratio(bu) - g(bu)).abs() < 1e-13, "{bu}: {}", ratio(bu));
    }
    // 折れ曲がりの手前の区間で、データの範囲（1）を 2/27 × h × 傾きの変化（3 次 Hermite の t = 2/3 の値）だけ下回る
    let lowest = samples(8.0, 10.0).map(ratio).fold(f64::INFINITY, f64::min);
    let bound = 2.0 / 27.0 * 2.0 * 0.05;
    assert!(
        (1.0 - lowest - bound).abs() < 1e-6,
        "lowest {lowest}, bound {bound}"
    );
    // 折れ曲がりの後の区間は、データの範囲の中
    for bu in samples(10.0, 12.0) {
        let r = ratio(bu);
        assert!((g(10.0)..=g(12.0)).contains(&r), "{bu}: {r}");
    }
}

#[test]
fn burnup_weights_match_the_python_comparison() {
    // examples/depletion/test_compare_burnup_interpolation.py の同じ点と値（ADR-0192 の比較の Python と同じ重み）
    let values = [
        1.3,
        1.069710678118655,
        1.1374213562373094,
        1.213606797749979,
        1.296227766016838,
        1.3882640687119285,
    ];
    let g = |bu: f64| values[CURVED.iter().position(|&p| p == bu).unwrap()];
    let t = burnup_table_of(&CURVED, g);
    let fresh = few(1300.0, 565.0, 900.0).k_inf;
    for (bu, want) in [
        (0.3, 1.161826406871193),
        (1.2, 1.103766188073716),
        (3.3, 1.1769920198453223),
        (7.7, 1.2626845558655815),
        (15.0, 1.3595423657587526),
    ] {
        let got = interpolate(&t, Some(&at_burnup(bu, 1300.0, 565.0, 900.0)), None)
            .unwrap()
            .remove(0)
            .k_inf
            / fresh;
        assert!(close(got, want), "burnup {bu}: {got} vs {want}");
    }
}

#[test]
fn negative_interpolated_cross_sections_are_errors() {
    // 熱群の吸収が 10 MWd/kg まで 0 で、そこから増える折れ線。折れ曲がりの手前の区間で Hermite は 0 を下回る
    let points: Vec<f64> = (1..=10).map(|i| 2.0 * f64::from(i)).collect();
    let mut t = burnup_table_of(&points, |_| 1.0);
    let per_burnup = t.states.len() / points.len();
    for (i, s) in t.states.iter_mut().enumerate() {
        s.absorption[1] = 0.01 * (points[i / per_burnup] - 10.0).max(0.0);
    }
    let err = interpolate(&t, Some(&at_burnup(9.0, 1300.0, 565.0, 900.0)), None).unwrap_err();
    assert!(err.contains("negative"), "{err}");
    // データの範囲の中なら通る
    assert!(interpolate(&t, Some(&at_burnup(12.0, 1300.0, 565.0, 900.0)), None).is_ok());
}

#[test]
fn negative_interpolated_form_functions_name_the_burnup_points() {
    // 形状関数の 1 つが 10 MWd/kg まで 0 で、そこから増える折れ線
    let points: Vec<f64> = (1..=10).map(|i| 2.0 * f64::from(i)).collect();
    let mut t = burnup_table_of(&points, |_| 1.0);
    let per_burnup = t.states.len() / points.len();
    for (i, s) in t.states.iter_mut().enumerate() {
        s.pins.as_mut().unwrap().form[1][0][1] = 0.01 * (points[i / per_burnup] - 10.0).max(0.0);
    }
    let err = interpolate(&t, Some(&at_burnup(9.0, 1300.0, 565.0, 900.0)), None).unwrap_err();
    assert!(err.contains("form") && err.contains("ADR-0192"), "{err}");
}

/// 軸の直積の重みつきの和（軸ごとの格子点と重みの直積を、最初の軸がいちばん遅く変わる順に、軸の順に掛けた重みで足す）。
/// 直積のテーブルの内挿がこれとビットで同じことを確かめる（ADR-0097、ADR-0192、ADR-0200）
fn product_sum(table: &Table, state: &State, rod: Option<&str>) -> Vec<FewGroup> {
    let terms: Vec<Terms> = table
        .axes
        .iter()
        .map(|axis| match axis {
            Axis::Burnup { points } => burnup_terms(points, state.burnup.unwrap()).unwrap(),
            Axis::Boron { ppm, .. } => bracket("boron", ppm, state.boron, |x| x).unwrap(),
            Axis::Moderator { temperature, .. } => {
                bracket("moderator", temperature, state.moderator, |x| x).unwrap()
            }
            Axis::Fuel { temperature, .. } => {
                cubic_terms("fuel", temperature, state.fuel, f64::sqrt).unwrap()
            }
            Axis::Rod { values } => vec![(
                values
                    .iter()
                    .position(|v| v == rod.unwrap_or(ROD_OUT))
                    .unwrap(),
                1.0,
            )],
        })
        .collect();
    let sub = table.sub_nodes();
    let mut out: Vec<FewGroup> = (0..sub).map(|q| scaled(&table.states[q], 0.0)).collect();
    let mut pick = vec![0usize; terms.len()];
    'corners: loop {
        let mut index = 0;
        let mut weight = 1.0;
        for ((axis, t), &p) in table.axes.iter().zip(&terms).zip(&pick) {
            let (i, w) = t[p];
            index = index * axis.len() + i;
            weight *= w;
        }
        if weight != 0.0 {
            for (q, acc) in out.iter_mut().enumerate() {
                accumulate(acc, weight, &table.states[index * sub + q]).unwrap();
            }
        }
        let mut k = terms.len();
        loop {
            if k == 0 {
                break 'corners;
            }
            k -= 1;
            pick[k] += 1;
            if pick[k] < terms[k].len() {
                break;
            }
            pick[k] = 0;
        }
    }
    out
}

/// `table` の 3 つの状態の軸に、公称の値（`nominal` の番号の格子点）を与えたもの
fn with_nominal(mut table: Table, nominal: [usize; 3]) -> Table {
    let mut n = nominal.into_iter();
    table.axes = table
        .axes
        .into_iter()
        .map(|a| match a.state_axis() {
            Some(_) => {
                let value = a.points()[n.next().unwrap()];
                a.with_nominal(value)
            }
            None => a,
        })
        .collect();
    table
}

#[test]
fn product_tables_interpolate_bit_for_bit_like_the_axis_product_sum() {
    let states = [
        (1000.0, 560.0, 700.0),
        (1291.0, 570.0, 565.0),
        (1650.0, 551.0, 899.0),
        (1300.0, 580.0, 900.0),
    ];
    // 公称の値のない直積と、公称の値を持つ直積（すべての成分の部分格子を持つ cut-HDMR）
    for t in [table(), with_nominal(table(), [1, 1, 0])] {
        for &(ppm, tm, tf) in &states {
            let s = state(ppm, tm, tf);
            assert_eq!(
                interpolate(&t, Some(&s), None).unwrap(),
                product_sum(&t, &s, None)
            );
        }
    }
    for t in [rodded_table(), with_nominal(rodded_table(), [2, 0, 1])] {
        for &(ppm, tm, tf) in &states {
            let s = state(ppm, tm, tf);
            for rod in [None, Some("aic"), Some("b4c")] {
                assert_eq!(
                    interpolate(&t, Some(&s), rod).unwrap(),
                    product_sum(&t, &s, rod)
                );
            }
        }
    }
    let curved = burnup_table_of(&CURVED, curve);
    for t in [curved.clone(), with_nominal(curved, [1, 1, 1])] {
        for bu in [0.0, 0.3, 0.5, 1.3, 7.5, 18.0] {
            for &(ppm, tm, tf) in &states {
                let s = at_burnup(bu, ppm, tm, tf);
                assert_eq!(
                    interpolate(&t, Some(&s), None).unwrap(),
                    product_sum(&t, &s, None),
                    "burnup {bu}"
                );
            }
        }
    }
}

#[test]
fn axes_interpolate_in_the_variable_of_their_transform() {
    // 燃料の温度の軸を温度そのもので内挿するテーブル（内挿の変数は軸の属性。ADR-0200）
    let mut t = table();
    t.axes[2] = Axis::Fuel {
        temperature: TF.to_vec(),
        nominal: None,
        transform: Transform::Linear,
    };
    let got = interpolate(&t, Some(&state(1300.0, 565.0, 700.0)), None)
        .unwrap()
        .remove(0);
    let w = (700.0 - 550.0) / (900.0 - 550.0);
    let want = (1.0 - w) * few(1300.0, 565.0, 550.0).k_inf + w * few(1300.0, 565.0, 900.0).k_inf;
    assert!(close(got.k_inf, want), "{} vs {want}", got.k_inf);
    // ホウ素の軸を √ppm で
    let mut t = table();
    t.axes[0] = Axis::Boron {
        ppm: PPM.to_vec(),
        b10: 0.199,
        nominal: None,
        transform: Transform::Sqrt,
    };
    let got = interpolate(&t, Some(&state(1000.0, 565.0, 550.0)), None)
        .unwrap()
        .remove(0);
    let w = (1000f64.sqrt() - 900f64.sqrt()) / (1300f64.sqrt() - 900f64.sqrt());
    let want = (1.0 - w) * few(900.0, 565.0, 550.0).k_inf + w * few(1300.0, 565.0, 550.0).k_inf;
    assert!(close(got.k_inf, want), "{} vs {want}", got.k_inf);
}

const TF3: [f64; 3] = [550.0, 565.0, 580.0];

/// (ppm, 減速材の温度, √燃料の温度) の公称（1300、565、565）からの偏差の多項式で、2 軸の積の項を持ち、3 軸の積の項の係数が
/// `triple` の値
fn pair_value(triple: f64, ppm: f64, tm: f64, tf: f64) -> f64 {
    let (b, m, f) = (ppm - 1300.0, tm - 565.0, tf.sqrt() - 565f64.sqrt());
    1.0 + 1.0e-4 * b
        + 2.0e-3 * m
        + 3.0e-2 * f
        + 1.0e-7 * b * m
        + 2.0e-6 * b * f
        + 4.0e-5 * m * f
        + triple * b * m * f
}

/// ホウ素、減速材、燃料が 3 点ずつ（公称は中央）で、`keep` を満たす点（公称の点の番号 1 の数で選ぶ）だけを持つテーブル。状態は
/// `pair_value(triple, …)` 倍
fn subgrid_table(triple: f64, keep: impl Fn(&[usize]) -> bool) -> Table {
    let axes = vec![
        Axis::boron(PPM.to_vec(), 0.199).with_nominal(1300.0),
        Axis::moderator(TM.to_vec(), vec![0.76972, 0.74271, 0.71190]).with_nominal(565.0),
        Axis::fuel(TF3.to_vec()).with_nominal(565.0),
    ];
    let points: Vec<Vec<usize>> = product_points(&axes)
        .into_iter()
        .filter(|p| keep(p))
        .collect();
    let states = points
        .iter()
        .map(|p| {
            curved(
                pair_value(triple, PPM[p[0]], TM[p[1]], TF3[p[2]]),
                1300.0,
                565.0,
                565.0,
            )
        })
        .collect();
    let table = Table {
        axes,
        split: 1,
        points,
        states,
        heavy_metal: None,
        micro: None,
        histories: Vec::new(),
        kinetics: None,
        cooling: None,
    };
    table.validate().unwrap();
    table
}

/// 公称でない軸が 2 つ以下の点（2 軸の組すべての部分格子。ADR-0200）
fn pairs(p: &[usize]) -> bool {
    p.contains(&1)
}

#[test]
fn cut_hdmr_reproduces_pair_terms_and_drops_only_the_triple_term() {
    let template = few(1300.0, 565.0, 565.0).k_inf;
    for (ppm, tm, tf) in [
        (1000.0, 560.0, 560.0),
        (1291.0, 570.0, 565.0),
        (1650.0, 551.0, 579.0),
        (900.0, 550.0, 550.0),
    ] {
        // 2 軸の積の項までの関数は、どの状態でも再現する
        let t = subgrid_table(0.0, pairs);
        let got = interpolate(&t, Some(&state(ppm, tm, tf)), None)
            .unwrap()
            .remove(0)
            .k_inf;
        let want = pair_value(0.0, ppm, tm, tf) * template;
        assert!((got - want).abs() < 1e-13 * want, "{got} vs {want}");
        // 3 軸の積の項だけを落とす
        let triple = 1.0e-9;
        let t = subgrid_table(triple, pairs);
        let got = interpolate(&t, Some(&state(ppm, tm, tf)), None)
            .unwrap()
            .remove(0)
            .k_inf;
        let dropped = triple * (ppm - 1300.0) * (tm - 565.0) * (tf.sqrt() - 565f64.sqrt());
        let want = (pair_value(triple, ppm, tm, tf) - dropped) * template;
        assert!((got - want).abs() < 1e-13 * want, "{got} vs {want}");
    }
}

#[test]
fn cut_hdmr_points_return_the_table_states() {
    let t = subgrid_table(1.0e-9, pairs);
    assert_eq!(t.points.len(), 19);
    for (p, s) in t.points.iter().zip(&t.states) {
        let got = interpolate(&t, Some(&state(PPM[p[0]], TM[p[1]], TF3[p[2]])), None)
            .unwrap()
            .remove(0);
        assert_eq!(&got, s, "{p:?}");
    }
}

#[test]
fn rod_values_assemble_from_their_own_nominal_state() {
    // 制御棒を入れた組（公称の状態を制御棒を入れた状態にした別の組。ADR-0200）は 1 軸ずつの分岐だけを持ち、吸収を半分にした値
    let out = subgrid_table(0.0, pairs);
    let rodded = subgrid_table(0.0, |p| p.iter().filter(|&&i| i != 1).count() <= 1);
    let mut t = out.clone();
    t.axes.push(Axis::Rod {
        values: vec!["out".into(), "aic".into()],
    });
    let mut rows: Vec<(Vec<usize>, FewGroup)> = out
        .points
        .iter()
        .zip(&out.states)
        .map(|(p, s)| ([p.as_slice(), &[0]].concat(), s.clone()))
        .chain(rodded.points.iter().zip(&rodded.states).map(|(p, s)| {
            let mut s = s.clone();
            s.absorption.iter_mut().for_each(|a| *a *= 0.5);
            ([p.as_slice(), &[1]].concat(), s)
        }))
        .collect();
    rows.sort_by(|a, b| a.0.cmp(&b.0));
    (t.points, t.states) = rows.into_iter().unzip();
    t.validate().unwrap();
    // 公称の状態で、制御棒を入れた組の公称の状態を返す
    let nominal = interpolate(&t, Some(&state(1300.0, 565.0, 565.0)), Some("aic"))
        .unwrap()
        .remove(0);
    assert_eq!(nominal, t.states[t.point_index(&[1, 1, 1, 1]).unwrap()]);
    // 制御棒を入れた組の 1 軸ずつの和（2 軸の積の項を落とす）
    let (ppm, tm, tf) = (1000.0, 560.0, 560.0);
    let got = interpolate(&t, Some(&state(ppm, tm, tf)), Some("aic"))
        .unwrap()
        .remove(0);
    let one_axis = pair_value(0.0, ppm, 565.0, 565.0)
        + pair_value(0.0, 1300.0, tm, 565.0)
        + pair_value(0.0, 1300.0, 565.0, tf)
        - 2.0;
    let base = few(1300.0, 565.0, 565.0);
    let want = 0.5 * one_axis * base.absorption[1];
    assert!(
        (got.absorption[1] - want).abs() < 1e-13 * want,
        "{} vs {want}",
        got.absorption[1]
    );
    // 制御棒のない組は 2 軸の組まで
    let got = interpolate(&t, Some(&state(ppm, tm, tf)), None)
        .unwrap()
        .remove(0);
    let want = pair_value(0.0, ppm, tm, tf) * base.absorption[1];
    assert!((got.absorption[1] - want).abs() < 1e-13 * want);
}

#[test]
fn negative_cut_hdmr_assemblies_are_errors() {
    // 熱群の吸収が、2 軸の組の点で小さく、1 軸ずつの分岐と公称で大きいと、3 つの軸の端（解いていない格子点）で負になる
    let mut t = subgrid_table(0.0, pairs);
    for (p, s) in t.points.iter().zip(t.states.iter_mut()) {
        let off = p.iter().filter(|&&i| i != 1).count();
        s.absorption[1] = if off == 2 { 0.001 } else { 0.08 };
    }
    let err = interpolate(&t, Some(&state(900.0, 550.0, 550.0)), None).unwrap_err();
    assert!(
        err.contains("negative") && err.contains("ADR-0200"),
        "{err}"
    );
    // 部分格子の点では負にならない
    assert!(interpolate(&t, Some(&state(900.0, 550.0, 565.0)), None).is_ok());
}

/// 燃焼度の点 CURVED のテーブルで、ホウ素 3 点 × 減速材 1 点 × 燃料 2 点（公称は 1300 ppm、565 K、550 K）の状態を解き、公称の状態は
/// すべての燃焼度の点で、ほかの点は `thinned` の燃焼度の点だけで解いたもの。状態は few(1300, 565, 550) の
/// `shape(燃焼度) + c × quadratic(燃焼度)` 倍（c は点ごとの係数で、公称で 0）
fn thinned_table(thinned: &[f64], shape: impl Fn(f64) -> f64) -> Table {
    let axes = vec![
        Axis::Burnup {
            points: CURVED.to_vec(),
        },
        Axis::boron(PPM.to_vec(), 0.199).with_nominal(1300.0),
        Axis::moderator(vec![565.0], vec![0.74271]).with_nominal(565.0),
        Axis::fuel(TF.to_vec()).with_nominal(550.0),
    ];
    let points: Vec<Vec<usize>> = product_points(&axes)
        .into_iter()
        .filter(|p| p[1..] == [1, 0, 0] || thinned.contains(&CURVED[p[0]]))
        .collect();
    let states = points
        .iter()
        .map(|p| {
            let bu = CURVED[p[0]];
            curved(
                shape(bu) + coefficient(&p[1..]) * quadratic(bu),
                1300.0,
                565.0,
                550.0,
            )
        })
        .collect();

    Table {
        axes,
        split: 1,
        points,
        states,
        heavy_metal: None,
        micro: None,
        histories: Vec::new(),
        kinetics: None,
        cooling: None,
    }
}

/// 点（ホウ素、減速材、燃料の番号）の、燃焼度の 2 次の差の係数（公称で 0）
fn coefficient(p: &[usize]) -> f64 {
    0.1 * (p[0] as f64 - 1.0) + 0.05 * p[2] as f64 + 0.02 * (p[0] as f64 - 1.0) * p[2] as f64
}

fn quadratic(bu: f64) -> f64 {
    0.01 * bu - 2.0e-4 * bu * bu
}

/// 燃焼度の点の間で 2 次でない形（公称の状態だけが持つ曲がり）
fn wiggle(bu: f64) -> f64 {
    if bu == 0.0 {
        1.3
    } else {
        1.0 + 0.03 * bu - 8.0e-4 * bu * bu + 0.02 * (bu / 3.0).sin()
    }
}

#[test]
fn off_nominal_states_add_their_difference_from_the_nominal_state_over_their_burnups() {
    // 公称でない状態は 0、0.5、5、18 MWd/kg だけで解く（ADR-0200 の案 B2）
    let t = thinned_table(&[0.0, 0.5, 5.0, 18.0], wiggle);
    t.validate().unwrap();
    let template = few(1300.0, 565.0, 550.0).k_inf;
    for (i, ppm) in PPM.into_iter().enumerate() {
        for (k, tf) in TF.into_iter().enumerate() {
            let c = coefficient(&[i, 0, k]);
            // 公称の状態だけを解いた燃焼度の点では、公称の状態の値に 2 次の差（3 点の放物線の傾きの Hermite が再現する）を足した値
            for bu in [2.0, 10.0] {
                let got = interpolate(&t, Some(&at_burnup(bu, ppm, 565.0, tf)), None)
                    .unwrap()
                    .remove(0)
                    .k_inf;
                let want = (wiggle(bu) + c * quadratic(bu)) * template;
                assert!(
                    (got - want).abs() < 1e-13 * want,
                    "{bu} {ppm} {tf}: {got} vs {want}"
                );
            }
            // 解いた燃焼度の点では、テーブルの状態そのもの
            for bu in [0.5, 5.0, 18.0] {
                let got = interpolate(&t, Some(&at_burnup(bu, ppm, 565.0, tf)), None)
                    .unwrap()
                    .remove(0);
                let p = t.point_index(&[CURVED.iter().position(|&b| b == bu).unwrap(), i, 0, k]);
                assert_eq!(got, t.states[p.unwrap()], "{bu} {ppm} {tf}");
            }
        }
    }
    // 燃焼度の点の間でも、公称の状態の内挿に差の内挿を足したもの
    let bu = 7.5;
    let s = at_burnup(bu, 1000.0, 565.0, 700.0);
    let nominal = interpolate(&t, Some(&at_burnup(bu, 1300.0, 565.0, 550.0)), None)
        .unwrap()
        .remove(0)
        .k_inf;
    let got = interpolate(&t, Some(&s), None).unwrap().remove(0).k_inf;
    let w = (700f64.sqrt() - 550f64.sqrt()) / (900f64.sqrt() - 550f64.sqrt());
    let b = (1000.0 - 900.0) / 400.0;
    // ホウ素と燃料の係数の多重線形の内挿
    let c = (1.0 - b) * ((1.0 - w) * coefficient(&[0, 0, 0]) + w * coefficient(&[0, 0, 1]))
        + b * ((1.0 - w) * coefficient(&[1, 0, 0]) + w * coefficient(&[1, 0, 1]));
    let want = nominal + c * quadratic(bu) * template;
    assert!((got - want).abs() < 1e-13 * want, "{got} vs {want}");
}

#[test]
fn burnups_outside_the_points_of_an_off_nominal_state_are_errors() {
    // 公称でない状態を 10 MWd/kg までしか解いていないテーブルは、その先の燃焼度で外挿しない（読み込みの検査も止める）
    let t = thinned_table(&[0.0, 0.5, 5.0, 10.0], wiggle);
    assert!(t.validate().is_err());
    assert!(interpolate(&t, Some(&at_burnup(15.0, 1300.0, 565.0, 550.0)), None).is_ok());
    let err = interpolate(&t, Some(&at_burnup(15.0, 1000.0, 565.0, 550.0)), None).unwrap_err();
    assert!(err.contains("burnup"), "{err}");
}

#[test]
fn negative_differences_from_the_nominal_state_are_errors() {
    // 公称の状態の熱群の吸収が燃焼度 2 の点だけ小さく、公称でない状態の差（0、0.5、5、18 で解いた）が大きく負だと、燃焼度 2 で
    // 公称の状態に差を足した値が負になる。燃焼度の重みは凸結合でも、差の内挿は公称の状態を負の重みで引く
    let mut t = thinned_table(&[0.0, 0.5, 5.0, 18.0], |_| 1.0);
    for (p, s) in t.points.iter().zip(t.states.iter_mut()) {
        let nominal = p[1..] == [1, 0, 0];
        s.absorption[1] = match (nominal, CURVED[p[0]] == 2.0) {
            (true, true) => 1.0e-4,
            (true, false) => 0.01,
            (false, _) => 0.001,
        };
    }
    t.validate().unwrap();
    let err = interpolate(&t, Some(&at_burnup(2.0, 900.0, 565.0, 550.0)), None).unwrap_err();
    assert!(
        err.contains("negative") && err.contains("ADR-0200"),
        "{err}"
    );
}

/// 状態ごとに、少数群定数の項目から作った 1 核種の微視的断面積（捕獲の経路は吸収）と、燃焼度の点ごとの数密度 `density(燃焼度)` を
/// 足したテーブル
fn with_micro(mut t: Table, density: impl Fn(f64) -> f64) -> Table {
    let Axis::Burnup { points } = t.axes[0].clone() else {
        panic!("the test tables put the burnup axis first")
    };
    t.micro = Some(xlc_format::table::Micro {
        nuclides: vec!["U-235".into()],
        densities: points.iter().map(|&b| vec![density(b)]).collect(),
        states: t
            .states
            .iter()
            .map(|s| {
                vec![MicroNuclide {
                    total: s.total.clone(),
                    nu_fission: s.nu_fission.clone(),
                    fission: s.fission.clone(),
                    scatter: s.scatter.clone(),
                    reactions: vec![xlc_format::table::MicroReaction {
                        mt: 102,
                        channels: vec![s.absorption.clone()],
                    }],
                    thermal_fission: Some(0.9),
                    derivative: None,
                }]
            })
            .collect(),
    });
    t.validate().unwrap();
    t
}

#[test]
fn micro_cross_sections_take_the_weights_of_the_few_group_constants() {
    // 線形の燃焼度の点と、3 次 Hermite の燃焼度の点（ADR-0192）
    for (t, last) in [
        (with_micro(burnup_table(), curve), 20.0),
        (with_micro(burnup_table_of(&CURVED, curve), curve), 18.0),
    ] {
        let Axis::Burnup { points } = &t.axes[0] else {
            panic!("the test tables put the burnup axis first")
        };
        let densities = &t.micro.as_ref().unwrap().densities;
        for (bu, ppm, tm, tf) in [
            (0.0, 1300.0, 565.0, 550.0),
            (0.3, 1000.0, 560.0, 700.0),
            (3.0, 1700.0, 580.0, 900.0),
            (10.0, 1450.0, 551.0, 600.0),
            (14.0, 900.0, 575.0, 820.0),
            (last, 1300.0, 565.0, 565.0),
        ] {
            let s = at_burnup(bu, ppm, tm, tf);
            let f = interpolate(&t, Some(&s), None).unwrap().remove(0);
            let m = interpolate_micro(&t, &s, None).unwrap();
            let [n] = &m.nuclides[..] else {
                panic!("one nuclide")
            };
            assert_eq!(n.total, f.total, "{bu}");
            assert_eq!(n.nu_fission, f.nu_fission);
            assert_eq!(n.fission, f.fission);
            assert_eq!(n.scatter, f.scatter);
            assert_eq!(n.reactions[0].mt, 102);
            assert_eq!(n.reactions[0].channels, std::slice::from_ref(&f.absorption));
            assert!((n.thermal_fission.unwrap() - 0.9).abs() <= 1e-15);
            // 数密度は公称の履歴の値なので、燃焼度だけで内挿する
            let want: f64 = burnup_terms(points, bu)
                .unwrap()
                .iter()
                .map(|&(k, h)| h * densities[k][0])
                .sum();
            assert_eq!(m.densities, [want], "{bu}");
        }
    }
}

#[test]
fn density_derivatives_take_the_weights_of_the_few_group_constants() {
    // 数密度の微分（ADR-0221）は σᵢ と同じ重みで状態と燃焼度で内挿する。微分は負になりうる（負の値で内挿のエラーにしない）
    for (t, last) in [
        (with_micro(burnup_table(), curve), 20.0),
        (with_micro(burnup_table_of(&CURVED, curve), curve), 18.0),
    ] {
        let mut t = t;
        let states = t.states.clone();
        for (state, few) in t.micro.as_mut().unwrap().states.iter_mut().zip(&states) {
            let neg = |v: &[f64]| v.iter().map(|x| -2.0 * x).collect::<Vec<_>>();
            // 反応率の微分（ADR-0215）は −2 σᵢ。境界より下の核分裂は断面積（割合 × 群の和の σf）の微分
            let sigma = state[0].clone();
            state[0].derivative = Some(xlc_format::table::MicroDerivative {
                total: neg(&few.total),
                nu_fission: neg(&few.nu_fission),
                fission: neg(&few.fission),
                scatter: few.scatter.iter().map(|r| neg(r)).collect(),
                rate: xlc_format::table::RateResidual {
                    fission: neg(&sigma.fission),
                    reactions: sigma
                        .reactions
                        .iter()
                        .map(|r| xlc_format::table::MicroReaction {
                            mt: r.mt,
                            channels: r.channels.iter().map(|c| neg(c)).collect(),
                        })
                        .collect(),
                    thermal_fission: sigma
                        .thermal_fission
                        .map(|f| -2.0 * f * sigma.fission.iter().sum::<f64>()),
                },
            });
        }
        t.validate().unwrap();
        for (bu, ppm, tm, tf) in [
            (0.0, 1300.0, 565.0, 550.0),
            (3.0, 1700.0, 580.0, 900.0),
            (14.0, 900.0, 575.0, 820.0),
            (last, 1300.0, 565.0, 565.0),
        ] {
            let s = at_burnup(bu, ppm, tm, tf);
            let f = interpolate(&t, Some(&s), None).unwrap().remove(0);
            let m = interpolate_micro(&t, &s, None).unwrap();
            let d = m.nuclides[0].derivative.as_ref().unwrap();
            let close = |a: &[f64], b: &[f64]| {
                a.iter()
                    .zip(b)
                    .all(|(x, y)| (x + 2.0 * y).abs() <= 1e-14 * y.abs().max(1e-300))
            };
            assert!(
                close(&d.total, &f.total),
                "{bu}: {:?} {:?}",
                d.total,
                f.total
            );
            assert!(close(&d.nu_fission, &f.nu_fission), "{bu}");
            assert!(close(&d.fission, &f.fission), "{bu}");
            for (a, b) in d.scatter.iter().zip(&f.scatter) {
                assert!(close(a, b), "{bu}");
            }
            // 反応率の微分も σᵢ と同じ重みで内挿する
            let sigma = &m.nuclides[0];
            assert!(close(&d.rate.fission, &sigma.fission), "{bu}");
            for (a, b) in d.rate.reactions.iter().zip(&sigma.reactions) {
                for (a, b) in a.channels.iter().zip(&b.channels) {
                    assert!(close(a, b), "{bu}");
                }
            }
            let below = sigma.thermal_fission.unwrap() * sigma.fission.iter().sum::<f64>();
            assert!(close(&[d.rate.thermal_fission.unwrap()], &[below]), "{bu}");
        }
    }
    // 微分を持たないテーブルは、内挿しても持たない
    let t = with_micro(burnup_table(), curve);
    let m = interpolate_micro(&t, &at_burnup(3.0, 1700.0, 580.0, 900.0), None).unwrap();
    assert!(m.nuclides[0].derivative.is_none());
}

#[test]
fn micro_cross_sections_need_the_micro_data_and_the_rod_axis() {
    let s = at_burnup(3.0, 1300.0, 565.0, 600.0);
    assert!(interpolate_micro(&burnup_table(), &s, None).is_err());
    let t = with_micro(burnup_table(), curve);
    assert!(interpolate_micro(&t, &s, Some("aic")).is_err());
    // 範囲の外の状態は少数群定数と同じくエラー
    assert!(interpolate_micro(&t, &at_burnup(25.0, 1300.0, 565.0, 600.0), None).is_err());
}

#[test]
fn thermal_fractions_near_one_stay_in_range_between_burnup_points() {
    // 境界より下の核分裂の断面積 a(u) が 5 MWd/kg までは小さく 10 MWd/kg から大きく、上の断面積は 1e-6 で一定。割合 a / (a + 1e-6) は
    // 10 と 18 MWd/kg で 1 に近く、割合そのものを 3 次 Hermite（負の重み）で内挿すると 1 を超える。断面積で内挿した比は超えない
    let below = |u: f64| if u <= 5.0 { 1.0e-5 } else { 1.0 };
    let mut t = with_micro(burnup_table_of(&CURVED, curve), curve);
    let points = t.points.clone();
    for (p, state) in points.iter().zip(&mut t.micro.as_mut().unwrap().states) {
        let a = below(CURVED[p[0]]);
        state[0].fission = vec![1.0e-6, a];
        state[0].thermal_fission = Some(a / (a + 1.0e-6));
    }
    t.validate().unwrap();
    for bu in [11.0, 13.0, 15.0, 17.0] {
        let s = at_burnup(bu, 1300.0, 580.0, 900.0);
        let m = interpolate_micro(&t, &s, None).unwrap();
        let f = m.nuclides[0].thermal_fission.unwrap();
        assert!((0.0..=1.0).contains(&f), "{bu}: {f}");
        // 境界より下と全体の断面積の比
        let fission = &m.nuclides[0].fission;
        let want = fission[1] / (fission[0] + fission[1]);
        assert!((f - want).abs() <= 1e-15, "{bu}: {f} vs {want}");
    }
}

// 履歴の補正（ADR-0215）

#[test]
fn moderator_densities_follow_the_temperature_weights_of_the_axis_points() {
    let mut table = Table::single(few(1000.0, 565.0, 900.0));
    table.axes = vec![Axis::moderator(
        vec![560.0, 570.0, 600.0],
        vec![0.75, 0.74, 0.70],
    )];
    // 格子点では、その点の密度そのもの（重みが厳密に 1）
    assert_eq!(moderator_density(&table, 560.0).unwrap(), 0.75);
    assert_eq!(moderator_density(&table, 570.0).unwrap(), 0.74);
    assert_eq!(moderator_density(&table, 600.0).unwrap(), 0.70);
    // 点の間は、少数群定数と同じく温度について線形
    let between = moderator_density(&table, 580.0).unwrap();
    assert!(
        (between - (0.74 * 2.0 / 3.0 + 0.70 / 3.0)).abs() < 1.0e-15,
        "{between}"
    );
    // 範囲の外は外挿しない
    for t in [550.0, 610.0, f64::NAN] {
        let err = moderator_density(&table, t).unwrap_err();
        assert!(err.contains("moderator"), "{err}");
    }
    // 減速材の軸のないテーブル
    let err = moderator_density(&Table::single(few(1000.0, 565.0, 900.0)), 565.0).unwrap_err();
    assert!(err.contains("moderator"), "{err}");
}

/// 燃焼度について `slope` × 燃焼度の残り（2 群、辺 6 つ）
fn residual_at(burnup: f64, slope: f64) -> Residual {
    let v = |a: f64| a * slope * burnup;
    Residual {
        total: vec![v(1.0), v(-2.0)],
        absorption: vec![v(0.1), v(-0.3)],
        nu_fission: vec![v(-0.2), v(0.4)],
        fission: vec![v(-0.08), v(0.16)],
        chi: vec![v(0.01), v(-0.01)],
        diffusion: vec![v(5.0), v(-1.0)],
        scatter: vec![vec![v(0.9), v(-0.05)], vec![v(0.002), v(-1.7)]],
        adf: vec![vec![v(0.1), v(-0.2)]; 6],
    }
}

#[test]
fn residuals_are_interpolated_in_burnup_like_the_table_and_not_extrapolated() {
    let points = [0.0, 10.0, 20.0, 30.0];
    let rod = HistoryRod {
        rod: ROD_OUT.into(),
        burnup: points.to_vec(),
        untracked: points.iter().map(|&b| residual_at(b, 1.0e-4)).collect(),
        linearization: points.iter().map(|&b| residual_at(b, -3.0e-4)).collect(),
        densities: vec![Vec::new(); points.len()],
        micro: xlc_format::table::RateResiduals {
            untracked: Vec::new(),
            linearization: Vec::new(),
        },
    };
    // 点では、その点の残りそのもの（追わない核種の分と線形化の誤差のどちらも）
    for (k, &b) in points.iter().enumerate() {
        assert_eq!(
            interpolate_residual(&rod, &rod.untracked, b).unwrap(),
            rod.untracked[k]
        );
        assert_eq!(
            interpolate_residual(&rod, &rod.linearization, b).unwrap(),
            rod.linearization[k]
        );
    }
    // 燃焼度について線形な残りは、3 次 Hermite でも線形のまま（ADR-0192 の重み）
    let r = interpolate_residual(&rod, &rod.untracked, 17.0).unwrap();
    let expected = residual_at(17.0, 1.0e-4);
    let close = |a: f64, b: f64| (a - b).abs() <= 1.0e-15 * b.abs().max(1.0e-6);
    assert!(close(r.total[1], expected.total[1]), "{r:?}");
    assert!(close(r.scatter[1][0], expected.scatter[1][0]), "{r:?}");
    assert!(close(r.adf[5][1], expected.adf[5][1]), "{r:?}");
    assert!(close(r.diffusion[0], expected.diffusion[0]), "{r:?}");
    // 点の範囲の外は外挿しない
    for b in [-1.0, 30.5] {
        let err = interpolate_residual(&rod, &rod.untracked, b).unwrap_err();
        assert!(err.contains("burnup"), "{err}");
    }
}

#[test]
fn branch_density_differences_are_interpolated_like_the_residuals() {
    // ΔN̄ᵢ^j = N̄ᵢ^j − N̄ᵢ⁰ を R の点で作り、R と同じ重み（ADR-0192。R の点で）で引く。N̄ᵢ⁰ はテーブルの燃焼度の軸の点の値
    let axis = [0.0, 5.0, 10.0, 20.0, 30.0];
    let nominal: Vec<Vec<f64>> = axis
        .iter()
        .map(|&b| vec![1.0 - 0.01 * b, 0.02 * b])
        .collect();
    let points = [0.0, 10.0, 20.0, 30.0];
    // 燃焼度について線形な差（ΔN̄ = (2e-3 b, −1e-3 b)）
    let rod = HistoryRod {
        rod: ROD_OUT.into(),
        burnup: points.to_vec(),
        untracked: Vec::new(),
        linearization: Vec::new(),
        densities: points
            .iter()
            .map(|&b| vec![1.0 - 0.01 * b + 2.0e-3 * b, 0.02 * b - 1.0e-3 * b])
            .collect(),
        micro: xlc_format::table::RateResiduals {
            untracked: Vec::new(),
            linearization: Vec::new(),
        },
    };
    for &b in &points {
        let d = interpolate_branch_difference(&rod, &axis, &nominal, b).unwrap();
        assert!(
            (d[0] - 2.0e-3 * b).abs() < 1e-15 && (d[1] + 1.0e-3 * b).abs() < 1e-15,
            "{b}: {d:?}"
        );
    }
    let d = interpolate_branch_difference(&rod, &axis, &nominal, 17.0).unwrap();
    assert!((d[0] - 2.0e-3 * 17.0).abs() < 1e-15, "{d:?}");
    // 範囲の外は外挿しない。R の点がテーブルの燃焼度の軸の点でなければエラー
    assert!(interpolate_branch_difference(&rod, &axis, &nominal, 30.5).is_err());
    let err = interpolate_branch_difference(&rod, &[0.0, 10.0, 20.0, 25.0], &nominal[..4], 17.0)
        .unwrap_err();
    assert!(err.contains("burnup"), "{err}");
}

/// 同じ B-10 の個数密度になる、テーブルの B-10 の割合（0.199）の ppm（ADR-0237）
fn converted(ppm: f64, b10: f64) -> f64 {
    ppm * xlc_format::table::b10_density_ratio(b10, 0.199).unwrap()
}

#[test]
fn differing_b10_interpolates_at_the_ppm_of_the_same_b10_density() {
    let t = table();
    for (ppm, b10, tm, tf) in [
        (1300.0, 0.18, 580.0, 900.0),
        (1000.0, 0.18, 560.0, 700.0),
        (1000.0, 0.30, 551.0, 899.0),
    ] {
        let mut s = state(ppm, tm, tf);
        s.b10 = b10;
        let got = interpolate(&t, Some(&s), None).unwrap().remove(0);
        // 換算した ppm の状態を、テーブルの割合で引いたものとビットで同じ
        let same = interpolate(&t, Some(&state(converted(ppm, b10), tm, tf)), None)
            .unwrap()
            .remove(0);
        assert_eq!(got, same, "{ppm} ppm at {b10}");
        // 多重線形な量は、換算した ppm の値を再現する
        let want = few(converted(ppm, b10), tm, tf).k_inf;
        assert!(
            close(got.k_inf, want),
            "{ppm} ppm at {b10}: {} vs {want}",
            got.k_inf
        );
    }
}

#[test]
fn converted_ppm_outside_the_boron_points_is_an_error() {
    // 1300 ppm の 30 at% は、19.9 at% で約 1978 ppm で、テーブルの上端 1700 ppm の外
    let t = table();
    let mut s = state(1300.0, 565.0, 600.0);
    s.b10 = 0.30;
    let e = interpolate(&t, Some(&s), None).unwrap_err();
    assert!(e.contains("boron") && e.contains("b10"), "{e}");
    assert!(e.contains(&format!("{}", converted(1300.0, 0.30))), "{e}");
}

#[test]
fn table_b10_is_not_converted() {
    // B-10 のないテーブル（b10 0）は換算の比が決まらない（ADR-0237）が、状態も b10 0 なら換算せずに引き、b10 だけが違う同じテーブルを
    // その b10 で引いたものとビットで同じ
    let t = table();
    let mut none = t.clone();
    let Axis::Boron { b10, .. } = &mut none.axes[0] else {
        panic!("the first axis is boron");
    };
    *b10 = 0.0;
    let mut s = state(1291.0, 570.0, 565.0);
    let want = interpolate(&t, Some(&s), None).unwrap();
    s.b10 = 0.0;
    assert_eq!(interpolate(&none, Some(&s), None).unwrap(), want);
}

#[test]
fn zero_ppm_of_another_b10_outside_the_boron_points_names_the_conversion() {
    // 0 ppm は換算しても 0 ppm で、テーブルの下端 900 ppm の外。B-10 の割合が違えば、エラーに換算の値を添える（ADR-0237）
    let t = table();
    let mut s = state(0.0, 565.0, 600.0);
    s.b10 = 0.18;
    let e = interpolate(&t, Some(&s), None).unwrap_err();
    assert!(e.contains("ADR-0237") && e.contains("b10 0.18"), "{e}");
}

// 動特性と冷却計算のデータ（ADR-0238）

/// 状態ごとに少数群定数の項目から作った 1 族の動特性のデータ（遅発の生成は Σf、崩壊定数で割った量は全断面積、スペクトルは χ、群速度の
/// 逆数は拡散係数）と、燃焼度の点ごとに `density(燃焼度)` から作った 2 核種の冷却計算の断面積を足したテーブル
fn with_kinetics(mut t: Table, sigma: impl Fn(f64) -> f64) -> Table {
    let Axis::Burnup { points } = t.axes[0].clone() else {
        panic!("the test tables put the burnup axis first")
    };
    t.kinetics = Some(xlc_format::table::Kinetics {
        states: t
            .states
            .iter()
            .map(|s| xlc_format::table::KineticsState {
                delayed_nu_fission: vec![s.fission.clone()],
                delayed_precursors: vec![s.total.clone()],
                delayed_chi: vec![s.chi.clone()],
                inverse_velocity: s.diffusion.clone(),
            })
            .collect(),
    });
    t.cooling = Some(xlc_format::table::Cooling {
        nuclides: vec!["U-235".into(), "Cs-133".into()],
        points: points
            .iter()
            .map(|&b| {
                vec![
                    xlc_format::table::CoolingNuclide {
                        fission: 40.0 * sigma(b),
                        thermal_fission: Some(0.9 - b / 1000.0),
                        reactions: vec![xlc_format::table::CoolingReaction {
                            mt: 102,
                            channels: vec![8.0 * sigma(b), 0.5 * sigma(b)],
                        }],
                    },
                    xlc_format::table::CoolingNuclide {
                        fission: 0.0,
                        thermal_fission: None,
                        reactions: vec![xlc_format::table::CoolingReaction {
                            mt: 102,
                            channels: vec![2.6 * sigma(b)],
                        }],
                    },
                ]
            })
            .collect(),
    });
    t.validate().unwrap();
    t
}

#[test]
fn kinetics_data_at_the_grid_points_are_the_lattice_values() {
    let t = with_kinetics(burnup_table_of(&CURVED, curve), curve);
    let kinetics = t.kinetics.as_ref().unwrap();
    for (p, point) in t.points.iter().enumerate() {
        let s = at_burnup(CURVED[point[0]], PPM[point[1]], TM[point[2]], TF[point[3]]);
        assert_eq!(
            interpolate_kinetics(&t, &s, None).unwrap(),
            kinetics.states[p],
            "{point:?}"
        );
    }
}

#[test]
fn kinetics_data_take_the_weights_of_the_few_group_constants() {
    for (t, last) in [
        (with_kinetics(burnup_table(), curve), 20.0),
        (with_kinetics(burnup_table_of(&CURVED, curve), curve), 18.0),
    ] {
        for (bu, ppm, tm, tf) in [
            (0.0, 1300.0, 565.0, 550.0),
            (0.3, 1000.0, 560.0, 700.0),
            (3.0, 1700.0, 580.0, 900.0),
            (14.0, 900.0, 575.0, 820.0),
            (last, 1300.0, 565.0, 565.0),
        ] {
            let s = at_burnup(bu, ppm, tm, tf);
            let f = interpolate(&t, Some(&s), None).unwrap().remove(0);
            let k = interpolate_kinetics(&t, &s, None).unwrap();
            assert_eq!(
                k.delayed_nu_fission,
                std::slice::from_ref(&f.fission),
                "{bu}"
            );
            assert_eq!(k.delayed_precursors, std::slice::from_ref(&f.total), "{bu}");
            assert_eq!(k.delayed_chi, std::slice::from_ref(&f.chi), "{bu}");
            assert_eq!(k.inverse_velocity, f.diffusion, "{bu}");
        }
    }
}

#[test]
fn decay_constants_keep_the_delayed_production_and_the_precursors_under_the_given_flux() {
    // 2 群で、族の崩壊定数で割った量が群ごとに違う（核種の混ざり方が群で違う）。λ = Σ_G ν_dΣf φ / Σ_G (ν_dΣf/λ) φ
    let kinetics = xlc_format::table::KineticsState {
        delayed_nu_fission: vec![vec![1.0e-5, 3.0e-4], vec![2.0e-5, 6.0e-4]],
        delayed_precursors: vec![vec![1.0e-3, 2.0e-2], vec![2.0e-5, 4.0e-4]],
        delayed_chi: vec![vec![0.9, 0.1], vec![1.0, 0.0]],
        inverse_velocity: vec![5.0e-8, 2.0e-6],
    };
    let flux = [3.0, 1.0];
    let lambda = decay_constants(&kinetics, &flux).unwrap();
    let want = [
        (1.0e-5 * 3.0 + 3.0e-4) / (1.0e-3 * 3.0 + 2.0e-2),
        (2.0e-5 * 3.0 + 6.0e-4) / (2.0e-5 * 3.0 + 4.0e-4),
    ];
    for (got, want) in lambda.iter().zip(want) {
        assert!((got - want).abs() <= 1e-15 * want, "{got} vs {want}");
    }
    // 群の数が違う中性子束と、遅発の生成のない族はエラー
    assert!(decay_constants(&kinetics, &[1.0]).is_err());
    let mut silent = kinetics.clone();
    silent.delayed_nu_fission[1] = vec![0.0, 0.0];
    silent.delayed_precursors[1] = vec![0.0, 0.0];
    assert!(decay_constants(&silent, &flux).is_err());
}

#[test]
fn cooling_cross_sections_are_interpolated_in_burnup_only() {
    let t = with_kinetics(burnup_table_of(&CURVED, curve), curve);
    let cooling = t.cooling.as_ref().unwrap();
    // 燃焼度の点では、その点の値（割合は比に戻すので丸めの差だけ）
    for (b, &bu) in CURVED.iter().enumerate() {
        let got = interpolate_cooling(&t, bu).unwrap();
        for (g, w) in got.iter().zip(&cooling.points[b]) {
            assert_eq!(g.fission, w.fission, "{bu}");
            assert_eq!(g.reactions, w.reactions, "{bu}");
            match (g.thermal_fission, w.thermal_fission) {
                (Some(a), Some(b)) => assert!((a - b).abs() <= 1e-15, "{a} vs {b}"),
                (a, b) => assert_eq!(a, b),
            }
        }
    }
    // 点の間は燃焼度の 3 次 Hermite（ADR-0192）
    for bu in [0.3, 3.0, 14.0] {
        let got = interpolate_cooling(&t, bu).unwrap();
        let terms = burnup_terms(&CURVED, bu).unwrap();
        let want: f64 = terms
            .iter()
            .map(|&(k, h)| h * cooling.points[k][1].reactions[0].channels[0])
            .sum();
        assert!(
            (got[1].reactions[0].channels[0] - want).abs() <= 1e-14 * want,
            "{bu}"
        );
        let fission: f64 = terms
            .iter()
            .map(|&(k, h)| h * cooling.points[k][0].fission)
            .sum();
        let below: f64 = terms
            .iter()
            .map(|&(k, h)| {
                let p = &cooling.points[k][0];
                h * p.fission * p.thermal_fission.unwrap()
            })
            .sum();
        assert!((got[0].fission - fission).abs() <= 1e-14 * fission, "{bu}");
        let fraction = got[0].thermal_fission.unwrap();
        assert!((fraction - below / fission).abs() <= 1e-14, "{bu}");
        assert_eq!(got[1].thermal_fission, None);
    }
}

#[test]
fn kinetics_and_cooling_data_need_the_data_and_the_burnup_range() {
    let s = at_burnup(3.0, 1300.0, 565.0, 600.0);
    assert!(interpolate_kinetics(&burnup_table(), &s, None).is_err());
    assert!(interpolate_cooling(&burnup_table(), 3.0).is_err());
    let t = with_kinetics(burnup_table(), curve);
    assert!(interpolate_kinetics(&t, &s, Some("aic")).is_err());
    assert!(interpolate_kinetics(&t, &at_burnup(25.0, 1300.0, 565.0, 600.0), None).is_err());
    assert!(interpolate_cooling(&t, 25.0).is_err());
    assert!(interpolate_cooling(&t, -1.0).is_err());
}
