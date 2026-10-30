use xlc_format::table::{
    Axis, MODERATOR_DENSITY, Micro, MicroDerivative, MicroReaction, ROD_OUT, Residual,
};

use super::*;

/// 2 群の少数群定数。`s` で値を変える（10 進で短く書けない値）
fn few(s: f64) -> FewGroup {
    let x = |v: f64| v * (1.0 + s / 7.0) / 3.0;
    FewGroup {
        bounds: vec![1.0e7, 0.625, 1.0e-5],
        k_inf: x(3.5),
        total: vec![x(0.55), x(1.3)],
        absorption: vec![x(0.01), x(0.1)],
        nu_fission: vec![x(0.007), x(0.14)],
        fission: vec![x(0.0028), x(0.057)],
        chi: vec![1.0 - s / 1000.0, s / 1000.0],
        diffusion: vec![x(1.4), x(0.4)],
        scatter: vec![vec![x(0.52), x(0.017)], vec![x(0.0015), x(1.2)]],
        adf: vec![vec![x(1.0), x(1.01)]; 6],
        vacuum_ratio: vec![None; 6],
        normalized_adf: vec![false; 6],
        height: None,
        pins: None,
    }
}

/// 核種の微視的断面積（2 群）。`derivative` なら数密度の微分も持つ
fn sigma(s: f64, derivative: bool) -> MicroNuclide {
    let x = |v: f64| v * (1.0 + s / 5.0) / 3.0;
    MicroNuclide {
        total: vec![x(11.0), x(130.0)],
        nu_fission: vec![x(1.1), x(170.0)],
        fission: vec![x(0.4), x(70.0)],
        scatter: vec![vec![x(9.0), x(0.1)], vec![x(0.01), x(8.0)]],
        reactions: vec![MicroReaction {
            mt: 102,
            channels: vec![vec![x(0.3), x(20.0)]],
        }],
        thermal_fission: Some(0.9),
        derivative: derivative.then(|| MicroDerivative {
            total: vec![x(9.0), x(110.0)],
            nu_fission: vec![x(1.0), x(-150.0)],
            fission: vec![x(0.35), x(-60.0)],
            scatter: vec![vec![x(8.5), x(0.12)], vec![x(0.02), x(7.0)]],
        }),
    }
}

const NUCLIDES: [&str; 2] = ["U-235", "Pu-240"];

/// 公称のテーブル: 燃焼度 0、10、20、ホウ素 1300 ppm、減速材 565 K と 0.743 g/cm³（公称）と 600 K と 0.661 g/cm³、燃料 565 K。
/// 公称でない状態は 0 と 20 で解く。Pu-240 は数密度の微分を持つ
fn nominal() -> Table {
    let axes = vec![
        Axis::Burnup {
            points: vec![0.0, 10.0, 20.0],
        },
        Axis::boron(vec![1300.0], 0.199).with_nominal(1300.0),
        Axis::moderator(vec![565.0, 600.0], vec![0.743, 0.661]).with_nominal(565.0),
        Axis::fuel(vec![565.0]).with_nominal(565.0),
    ];
    let points: Vec<Vec<usize>> = vec![
        vec![0, 0, 0, 0],
        vec![0, 0, 1, 0],
        vec![1, 0, 0, 0],
        vec![2, 0, 0, 0],
        vec![2, 0, 1, 0],
    ];
    let n = points.len();
    let table = Table {
        axes,
        split: 1,
        states: (0..n).map(|p| few(p as f64)).collect(),
        heavy_metal: Some(400.0 / 3.0),
        micro: Some(Micro {
            nuclides: NUCLIDES.iter().map(|&s| s.into()).collect(),
            densities: vec![
                vec![7.0e-4, 0.0],
                vec![6.0e-4, 1.0e-5],
                vec![5.0e-4, 2.0e-5],
            ],
            states: (0..n)
                .map(|p| vec![sigma(p as f64, false), sigma(2.0 * p as f64, true)])
                .collect(),
        }),
        histories: Vec::new(),
        points,
    };
    table.validate().unwrap();
    table
}

/// 履歴の分岐の自分の条件（600 K と 0.661 g/cm³）のテーブル: 燃焼度 0 と 20、状態 1 つ
fn branch() -> Table {
    let axes = vec![
        Axis::Burnup {
            points: vec![0.0, 20.0],
        },
        Axis::boron(vec![1300.0], 0.199).with_nominal(1300.0),
        Axis::moderator(vec![600.0], vec![0.661]).with_nominal(600.0),
        Axis::fuel(vec![565.0]).with_nominal(565.0),
    ];
    let table = Table {
        axes,
        split: 1,
        points: vec![vec![0, 0, 0, 0], vec![1, 0, 0, 0]],
        states: vec![few(1.5), few(4.5)],
        heavy_metal: Some(400.0 / 3.0),
        micro: Some(Micro {
            nuclides: NUCLIDES.iter().map(|&s| s.into()).collect(),
            densities: vec![vec![7.0e-4, 0.0], vec![5.2e-4, 2.6e-5]],
            states: (0..2)
                .map(|p| vec![sigma(10.0 + p as f64, false), sigma(20.0 + p as f64, false)])
                .collect(),
        }),
        histories: Vec::new(),
    };
    table.validate().unwrap();
    table
}

/// 丸めの範囲で等しい: 差が、足し引きした項の大きさ `scale` の 1e-14 倍以内（足す順が違えば、打ち消し合う項の丸めが残る）
fn near(a: f64, b: f64, scale: f64) -> bool {
    (a - b).abs() <= 1e-14 * scale
}

fn close(a: f64, b: f64) -> bool {
    near(a, b, a.abs().max(b.abs()))
}

#[test]
fn residuals_are_the_branch_constants_less_the_linearly_corrected_nominal_condition() {
    let n = nominal();
    let b = branch();
    let (h0, hb, rod) = residuals(&n, MODERATOR_DENSITY, ROD_OUT, &b).unwrap();
    assert_eq!((h0, hb), (0.743, 0.661));
    assert_eq!(rod.rod, ROD_OUT);
    // 公称のテーブルが条件の状態を解いた燃焼度の点
    assert_eq!(rod.burnup, vec![0.0, 20.0]);
    let micro = n.micro.as_ref().unwrap();
    let branch_micro = b.micro.as_ref().unwrap();
    // 公称のテーブルの条件の状態の点（燃焼度 0 と 20）と、分岐のテーブルの点
    for (k, (&p, &u)) in [1usize, 4].iter().zip(&[0usize, 2]).enumerate() {
        let f0 = &n.states[p];
        let fb = &b.states[k];
        let c = &micro.states[p];
        let d: Vec<f64> = (0..2)
            .map(|i| branch_micro.densities[k][i] - micro.densities[u][i])
            .collect();
        // U-235 は σᵢ、Pu-240 は数密度の微分（core の補正と同じ係数）
        let coefficient = |i: usize| -> MicroDerivative {
            c[i].derivative.clone().unwrap_or_else(|| MicroDerivative {
                total: c[i].total.clone(),
                nu_fission: c[i].nu_fission.clone(),
                fission: c[i].fission.clone(),
                scatter: c[i].scatter.clone(),
            })
        };
        assert!(c[1].derivative.is_some() && c[0].derivative.is_none());
        let r = &rod.residuals[k];
        for g in 0..2 {
            let sum = |pick: &dyn Fn(usize) -> f64| (0..2).map(|i| d[i] * pick(i)).sum::<f64>();
            let total = sum(&|i| coefficient(i).total[g]);
            let absorption =
                sum(&|i| coefficient(i).total[g] - coefficient(i).scatter[g].iter().sum::<f64>());
            // R = 分岐 − 公称 − 補正（丸めは 3 つの項の大きさの範囲）
            let expect = |r: f64, b: f64, n: f64, c: f64| {
                assert!(
                    near(r, b - n - c, b.abs() + n.abs() + c.abs()),
                    "{g}: {r} {b} {n} {c}"
                );
            };
            expect(r.total[g], fb.total[g], f0.total[g], total);
            expect(
                r.absorption[g],
                fb.absorption[g],
                f0.absorption[g],
                absorption,
            );
            expect(
                r.nu_fission[g],
                fb.nu_fission[g],
                f0.nu_fission[g],
                sum(&|i| coefficient(i).nu_fission[g]),
            );
            expect(
                r.fission[g],
                fb.fission[g],
                f0.fission[g],
                sum(&|i| coefficient(i).fission[g]),
            );
            for h in 0..2 {
                expect(
                    r.scatter[g][h],
                    fb.scatter[g][h],
                    f0.scatter[g][h],
                    sum(&|i| coefficient(i).scatter[g][h]),
                );
            }
            // 線形の補正を持たない項目は 2 つの定数の差
            assert_eq!(r.chi[g], fb.chi[g] - f0.chi[g]);
            assert_eq!(r.diffusion[g], fb.diffusion[g] - f0.diffusion[g]);
            for side in 0..6 {
                assert_eq!(r.adf[side][g], fb.adf[side][g] - f0.adf[side][g]);
            }
        }
    }
}

#[test]
fn the_nominal_condition_with_the_correction_and_the_residual_is_the_branch() {
    // core が履歴の分岐の点（w = 1、ノードの数密度が分岐のもの）で返す定数は、分岐の格子計算の定数（ADR-0215）
    let n = nominal();
    let b = branch();
    let (_, _, rod) = residuals(&n, MODERATOR_DENSITY, ROD_OUT, &b).unwrap();
    let micro = n.micro.as_ref().unwrap();
    for (k, (&p, &u)) in [1usize, 4].iter().zip(&[0usize, 2]).enumerate() {
        let d: Vec<f64> = (0..2)
            .map(|i| b.micro.as_ref().unwrap().densities[k][i] - micro.densities[u][i])
            .collect();
        let corrected = linear(&n.states[p], &micro.states[p], &d);
        let r: &Residual = &rod.residuals[k];
        let fb = &b.states[k];
        for g in 0..2 {
            assert!(close(corrected.total[g] + r.total[g], fb.total[g]));
            assert!(close(
                corrected.absorption[g] + r.absorption[g],
                fb.absorption[g]
            ));
            assert!(close(
                corrected.nu_fission[g] + r.nu_fission[g],
                fb.nu_fission[g]
            ));
            for h in 0..2 {
                assert!(close(
                    corrected.scatter[g][h] + r.scatter[g][h],
                    fb.scatter[g][h]
                ));
            }
        }
        // 補正は χ、拡散係数、ADF を変えない
        assert_eq!(corrected.chi, n.states[p].chi);
        assert_eq!(corrected.diffusion, n.states[p].diffusion);
        assert_eq!(corrected.adf, n.states[p].adf);
        // 補正の差が 0 なら定数はそのまま
        assert_eq!(
            linear(&n.states[p], &micro.states[p], &[0.0, 0.0]),
            n.states[p]
        );
    }
}

#[test]
fn history_corrections_are_written_to_a_valid_nominal_table() {
    let mut n = nominal();
    let correction =
        history_correction(&n, MODERATOR_DENSITY, &[(ROD_OUT.into(), branch())]).unwrap();
    assert_eq!(correction.variable, MODERATOR_DENSITY);
    assert_eq!((correction.nominal, correction.branch), (0.743, 0.661));
    assert_eq!(correction.rods.len(), 1);
    n.histories.push(correction);
    n.validate().unwrap();
}

#[test]
fn branch_conditions_differ_from_the_nominal_state_on_the_variable_axis_only() {
    let n = nominal();
    let fails = |change: &dyn Fn(&mut Table)| {
        let mut b = branch();
        change(&mut b);
        let err = residuals(&n, MODERATOR_DENSITY, ROD_OUT, &b)
            .unwrap_err()
            .to_string();
        // 折り返した文の空白が文に入らない
        assert!(!err.contains("  "), "{err}");
        err
    };
    // 公称のテーブルの格子点でない
    let err = fails(&|b| b.axes[2] = Axis::moderator(vec![610.0], vec![0.661]).with_nominal(610.0));
    assert!(err.contains("point"), "{err}");
    let err = fails(&|b| b.axes[2] = Axis::moderator(vec![600.0], vec![0.65]).with_nominal(600.0));
    assert!(err.contains("point"), "{err}");
    // 公称の状態そのもの
    let err = fails(&|b| b.axes[2] = Axis::moderator(vec![565.0], vec![0.743]).with_nominal(565.0));
    assert!(err.contains("moderator"), "{err}");
    // 変数の軸のほかが公称と違う
    let mut n2 = nominal();
    n2.axes[3] = Axis::fuel(vec![565.0, 600.0]).with_nominal(565.0);
    let mut b = branch();
    b.axes[3] = Axis::fuel(vec![600.0]).with_nominal(600.0);
    b.axes[2] = Axis::moderator(vec![565.0], vec![0.743]).with_nominal(565.0);
    let err = residuals(&n2, MODERATOR_DENSITY, ROD_OUT, &b)
        .unwrap_err()
        .to_string();
    assert!(err.contains("moderator"), "{err}");
    // 分岐のテーブルは状態 1 つ
    let err = fails(&|b| {
        b.axes[2] = Axis::moderator(vec![565.0, 600.0], vec![0.743, 0.661]).with_nominal(600.0)
    });
    assert!(err.contains("one"), "{err}");
    // 知らない変数
    let err = residuals(&n, "boron", ROD_OUT, &branch())
        .unwrap_err()
        .to_string();
    assert!(err.contains("variable"), "{err}");
}

#[test]
fn branch_burnups_are_those_solved_at_the_condition() {
    let n = nominal();
    let mut b = branch();
    b.axes[0] = Axis::Burnup {
        points: vec![0.0, 10.0, 20.0],
    };
    b.points = vec![vec![0, 0, 0, 0], vec![1, 0, 0, 0], vec![2, 0, 0, 0]];
    b.states.push(few(7.0));
    let micro = b.micro.as_mut().unwrap();
    micro.densities.insert(1, vec![6.1e-4, 1.2e-5]);
    micro.states.push(micro.states[1].clone());
    b.validate().unwrap();
    // 公称のテーブルは条件の状態を 10 で解いていない
    let err = residuals(&n, MODERATOR_DENSITY, ROD_OUT, &b)
        .unwrap_err()
        .to_string();
    assert!(err.contains("burnup"), "{err}");
}

#[test]
fn branch_tables_hold_the_same_micro_nuclides_and_heavy_metal() {
    let n = nominal();
    let mut b = branch();
    b.micro.as_mut().unwrap().nuclides[1] = "Pu-239".into();
    let err = residuals(&n, MODERATOR_DENSITY, ROD_OUT, &b)
        .unwrap_err()
        .to_string();
    assert!(err.contains("nuclides"), "{err}");
    let mut b = branch();
    b.heavy_metal = Some(401.0);
    let err = residuals(&n, MODERATOR_DENSITY, ROD_OUT, &b)
        .unwrap_err()
        .to_string();
    assert!(err.contains("heavy metal"), "{err}");
    let mut b = branch();
    b.micro = None;
    assert!(residuals(&n, MODERATOR_DENSITY, ROD_OUT, &b).is_err());
}

#[test]
fn rod_values_of_tables_without_a_rod_axis_are_out() {
    let err = residuals(&nominal(), MODERATOR_DENSITY, "aic", &branch())
        .unwrap_err()
        .to_string();
    assert!(err.contains("rod"), "{err}");
}

#[test]
fn the_input_names_the_tables_per_variable_and_rod_value() {
    let input = ResidualInput::parse(
        r#"
nominal = "n.msgpack"
output = "out.msgpack"
[[history]]
variable = "moderator-density"
branches = { out = "b.msgpack" }
"#,
    )
    .unwrap();
    assert_eq!(
        input,
        ResidualInput {
            nominal: "n.msgpack".into(),
            output: "out.msgpack".into(),
            histories: vec![HistoryEntry {
                variable: MODERATOR_DENSITY.into(),
                branches: vec![(ROD_OUT.into(), "b.msgpack".into())],
            }],
        }
    );
    for bad in [
        // 履歴がない
        "nominal = \"n\"\noutput = \"o\"\n",
        // 分岐のテーブルがない
        "nominal = \"n\"\noutput = \"o\"\n[[history]]\nvariable = \"moderator-density\"\nbranches = {}\n",
        // 知らない項目
        "nominal = \"n\"\noutput = \"o\"\nextra = 1\n[[history]]\nvariable = \"moderator-density\"\nbranches = { out = \"b\" }\n",
        // 出力が公称のテーブル（読んだファイルを上書きしない）
        "nominal = \"n\"\noutput = \"n\"\n[[history]]\nvariable = \"moderator-density\"\nbranches = { out = \"b\" }\n",
    ] {
        assert!(ResidualInput::parse(bad).is_err(), "{bad}");
    }
}

#[test]
fn run_writes_the_nominal_table_with_the_history_residuals() {
    let dir = std::env::temp_dir().join(format!("xlc-residual-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("n.msgpack"), nominal().to_msgpack().unwrap()).unwrap();
    std::fs::write(dir.join("b.msgpack"), branch().to_msgpack().unwrap()).unwrap();
    let input = dir.join("histories.toml");
    std::fs::write(
        &input,
        "nominal = \"n.msgpack\"\noutput = \"o.msgpack\"\n[[history]]\nvariable = \"moderator-density\"\n\
         branches = { out = \"b.msgpack\" }\n",
    )
    .unwrap();
    run(&input).unwrap();
    let written = Table::from_msgpack(&std::fs::read(dir.join("o.msgpack")).unwrap()).unwrap();
    let mut expected = nominal();
    expected.histories = vec![
        history_correction(&expected, MODERATOR_DENSITY, &[(ROD_OUT.into(), branch())]).unwrap(),
    ];
    assert_eq!(written, expected);
    // R を持つ公称のテーブルに重ねて足さない
    std::fs::write(dir.join("n.msgpack"), written.to_msgpack().unwrap()).unwrap();
    let err = run(&input).unwrap_err().to_string();
    assert!(err.contains("already"), "{err}");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn branch_values_on_axes_the_nominal_table_lacks_are_errors() {
    // 公称のテーブルに燃料の軸が無いと、分岐の燃料の温度を照らせず、温度の差の効きが R に混ざる
    let mut n = nominal();
    n.axes.remove(3);
    for p in &mut n.points {
        p.pop();
    }
    for axis in &mut n.axes {
        *axis = axis.clone().without_nominal();
    }
    n.axes[1] = n.axes[1].clone().with_nominal(1300.0);
    n.axes[2] = n.axes[2].clone().with_nominal(565.0);
    n.validate().unwrap();
    let mut b = branch();
    b.axes[3] = Axis::fuel(vec![900.0]).with_nominal(900.0);
    let err = residuals(&n, MODERATOR_DENSITY, ROD_OUT, &b)
        .unwrap_err()
        .to_string();
    assert!(err.contains("fuel"), "{err}");
}

#[test]
fn branch_conditions_with_the_nominal_density_are_errors() {
    // 温度だけが公称と違う格子点では Hᵇ = H⁰ で、core の重みが 0 で割る
    let mut n = nominal();
    n.axes[2] = Axis::moderator(vec![565.0, 600.0], vec![0.743, 0.743]).with_nominal(565.0);
    let mut b = branch();
    b.axes[2] = Axis::moderator(vec![600.0], vec![0.743]).with_nominal(600.0);
    let err = residuals(&n, MODERATOR_DENSITY, ROD_OUT, &b)
        .unwrap_err()
        .to_string();
    assert!(err.contains("density"), "{err}");
}
