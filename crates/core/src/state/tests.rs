use super::*;

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
        adf: vec![vec![v, 2.0 * v]; 4],
        vacuum_ratio: vec![None, Some(vec![0.5 * v / 3.0, 0.4 * v / 3.0]), None, None],
        normalized_adf: vec![false; 4],
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
    }
}

fn close(a: f64, b: f64) -> bool {
    ((a - b) / b).abs() < 1e-13
}

#[test]
fn grid_points_return_the_lattice_states() {
    let t = table();
    let got = interpolate(&t, Some(&state(1300.0, 580.0, 900.0))).unwrap();
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
        let got = interpolate(&t, Some(&state(ppm, tm, tf))).unwrap();
        let want = few(ppm, tm, tf);
        assert!(
            close(got.k_inf, want.k_inf),
            "{} vs {}",
            got.k_inf,
            want.k_inf
        );
        assert!(close(got.scatter[0][1], want.scatter[0][1]));
        assert!(close(got.adf[3][1], want.adf[3][1]));
        let (Some(g), Some(w)) = (&got.vacuum_ratio[1], &want.vacuum_ratio[1]) else {
            panic!("vacuum ratio")
        };
        assert!(close(g[0], w[0]));
        assert_eq!(got.vacuum_ratio[0], None);
        assert_eq!(got.bounds, want.bounds);
        assert_eq!(got.normalized_adf, want.normalized_adf);
    }
}

#[test]
fn fuel_temperature_is_interpolated_in_its_square_root() {
    // 燃料の温度だけを変えると、量は √T について線形に変わる
    let t = table();
    let got = interpolate(&t, Some(&state(1300.0, 565.0, 700.0))).unwrap();
    let w = (700f64.sqrt() - 550f64.sqrt()) / (900f64.sqrt() - 550f64.sqrt());
    let want = (1.0 - w) * few(1300.0, 565.0, 550.0).k_inf + w * few(1300.0, 565.0, 900.0).k_inf;
    assert!(close(got.k_inf, want));
}

#[test]
fn invalid_states_are_errors() {
    let t = table();
    let err = |s: Option<&State>, table: &Table| interpolate(table, s).unwrap_err();
    // 範囲の外（外挿しない）
    assert!(err(Some(&state(800.0, 565.0, 600.0)), &t).contains("boron"));
    assert!(err(Some(&state(1300.0, 590.0, 600.0)), &t).contains("moderator"));
    assert!(err(Some(&state(1300.0, 565.0, 1000.0)), &t).contains("fuel"));
    assert!(err(Some(&state(f64::NAN, 565.0, 600.0)), &t).contains("boron"));
    // B-10 の割合が違う
    let mut s = state(1300.0, 565.0, 600.0);
    s.b10 = 0.2;
    assert!(err(Some(&s), &t).contains("b10"));
    // 状態があるのに軸がない、状態がないのに軸がある
    let mut two = t.clone();
    two.axes.pop();
    two.states = two.states.iter().step_by(2).cloned().collect();
    assert!(err(Some(&state(1300.0, 565.0, 600.0)), &two).contains("fuel"));
    assert!(err(None, &t).contains("state"));
}

#[test]
fn single_point_axes_need_that_value() {
    let mut t = table();
    t.axes[2] = Axis::Fuel {
        temperature: vec![550.0],
    };
    t.states = t.states.iter().step_by(2).cloned().collect();
    let got = interpolate(&t, Some(&state(1300.0, 565.0, 550.0))).unwrap();
    assert_eq!(got, few(1300.0, 565.0, 550.0));
    assert!(interpolate(&t, Some(&state(1300.0, 565.0, 551.0))).is_err());
}

#[test]
fn tables_without_axes_need_no_state() {
    let single = Table::single(few(1300.0, 565.0, 565.0));
    assert_eq!(
        interpolate(&single, None).unwrap(),
        few(1300.0, 565.0, 565.0)
    );
    assert!(interpolate(&single, Some(&state(1300.0, 565.0, 565.0))).is_err());
}
