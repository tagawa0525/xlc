use super::*;
use xlc_format::fewgroup::Pins;

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
    Table {
        split: 1,
        heavy_metal: None,
        axes: vec![
            Axis::Boron {
                ppm: PPM.to_vec(),
                b10: 0.199,
            },
            Axis::Moderator {
                temperature: TM.to_vec(),
                density: vec![0.76972, 0.74271, 0.71190],
            },
            Axis::Fuel {
                temperature: TF.to_vec(),
            },
        ],
        states,
    }
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

#[test]
fn invalid_states_are_errors() {
    let t = table();
    let err = |s: Option<&State>, table: &Table| interpolate(table, s, None).unwrap_err();
    // 範囲の外（外挿しない）
    assert!(err(Some(&state(800.0, 565.0, 600.0)), &t).contains("boron"));
    assert!(err(Some(&state(1300.0, 590.0, 600.0)), &t).contains("moderator"));
    assert!(err(Some(&state(1300.0, 565.0, 1000.0)), &t).contains("fuel"));
    assert!(err(Some(&state(f64::NAN, 565.0, 600.0)), &t).contains("boron"));
    // B-10 の割合が違う
    let mut s = state(1300.0, 565.0, 600.0);
    s.b10 = 0.2;
    assert!(err(Some(&s), &t).contains("b10"));
    assert!(err(Some(&s), &t).contains("exactly"));
    // 状態があるのに軸がない、状態がないのに軸がある
    let mut two = t.clone();
    two.axes.pop();
    two.states = two.states.iter().step_by(2).cloned().collect();
    assert!(err(Some(&state(1300.0, 565.0, 600.0)), &two).contains("fuel"));
    // 軸がないことは、ほかの軸の範囲の外より先に示す
    assert!(err(Some(&state(800.0, 565.0, 600.0)), &two).contains("fuel"));
    assert!(err(None, &t).contains("state"));
}

#[test]
fn single_point_axes_need_that_value() {
    let mut t = table();
    t.axes[2] = Axis::Fuel {
        temperature: vec![550.0],
    };
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
