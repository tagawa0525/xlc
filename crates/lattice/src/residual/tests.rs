use xlc_format::table::{
    Axis, MODERATOR_DENSITY, Micro, MicroDerivative, MicroNuclide, MicroReaction, ROD_OUT, Residual,
};

use super::*;

/// 公称の状態の定数 `state` に、線形の補正を足したもの
fn linear(state: &FewGroup, micro: &[MicroNuclide], delta: &[f64]) -> FewGroup {
    let mut f = state.clone();
    add_linear_correction(&mut f, micro, delta);
    f
}

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

/// 履歴の分岐の組成のうち追う核種だけを公称の組成に足したハイブリッドを、同じ条件で解いたテーブル: 分岐と同じ軸、点、追う核種、
/// 均質化した数密度で、定数と微視的断面積が違う
fn hybrid() -> Table {
    let mut table = branch();
    table.states = vec![few(1.5), few(3.5)];
    table.micro.as_mut().unwrap().states = (0..2)
        .map(|p| vec![sigma(30.0 + p as f64, false), sigma(40.0 + p as f64, false)])
        .collect();
    table.validate().unwrap();
    table
}

/// 制御棒の値 out の、分岐とハイブリッドのテーブルの組
fn pair() -> BranchTables {
    BranchTables {
        branch: branch(),
        hybrid: hybrid(),
    }
}

/// 丸めの範囲で等しい: 差が、足し引きした項の大きさ `scale` の 1e-14 倍以内（足す順が違えば、打ち消し合う項の丸めが残る）
fn near(a: f64, b: f64, scale: f64) -> bool {
    (a - b).abs() <= 1e-14 * scale
}

fn close(a: f64, b: f64) -> bool {
    near(a, b, a.abs().max(b.abs()))
}

#[test]
fn linearizations_are_the_hybrid_less_the_linearly_corrected_nominal_condition() {
    let n = nominal();
    let b = branch();
    let h = hybrid();
    let (h0, hb, rod) = residuals(&n, MODERATOR_DENSITY, ROD_OUT, &pair()).unwrap();
    assert_eq!((h0, hb), (0.743, 0.661));
    assert_eq!(rod.rod, ROD_OUT);
    // 公称のテーブルが条件の状態を解いた燃焼度の点
    assert_eq!(rod.burnup, vec![0.0, 20.0]);
    let micro = n.micro.as_ref().unwrap();
    let branch_micro = b.micro.as_ref().unwrap();
    // 公称のテーブルの条件の状態の点（燃焼度 0 と 20）と、分岐のテーブルの点
    for (k, (&p, &u)) in [1usize, 4].iter().zip(&[0usize, 2]).enumerate() {
        let f0 = &n.states[p];
        // L は、追う核種だけが分岐と同じハイブリッドの定数から、線形の補正を足した公称の条件の定数を引いたもの
        let fb = &h.states[k];
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
        let r = &rod.linearization[k];
        for g in 0..2 {
            let sum = |pick: &dyn Fn(usize) -> f64| (0..2).map(|i| d[i] * pick(i)).sum::<f64>();
            let total = sum(&|i| coefficient(i).total[g]);
            let absorption =
                sum(&|i| coefficient(i).total[g] - coefficient(i).scatter[g].iter().sum::<f64>());
            // L = ハイブリッド − 公称 − 補正（丸めは 3 つの項の大きさの範囲）
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
fn untracked_parts_are_the_branch_less_the_hybrid() {
    let (_, _, rod) = residuals(&nominal(), MODERATOR_DENSITY, ROD_OUT, &pair()).unwrap();
    let (b, h) = (branch(), hybrid());
    assert_eq!(rod.untracked.len(), 2);
    for (k, u) in rod.untracked.iter().enumerate() {
        let (fb, fh) = (&b.states[k], &h.states[k]);
        for g in 0..2 {
            assert_eq!(u.total[g], fb.total[g] - fh.total[g]);
            assert_eq!(u.absorption[g], fb.absorption[g] - fh.absorption[g]);
            assert_eq!(u.nu_fission[g], fb.nu_fission[g] - fh.nu_fission[g]);
            assert_eq!(u.fission[g], fb.fission[g] - fh.fission[g]);
            assert_eq!(u.chi[g], fb.chi[g] - fh.chi[g]);
            assert_eq!(u.diffusion[g], fb.diffusion[g] - fh.diffusion[g]);
            for h in 0..2 {
                assert_eq!(u.scatter[g][h], fb.scatter[g][h] - fh.scatter[g][h]);
            }
            for side in 0..6 {
                assert_eq!(u.adf[side][g], fb.adf[side][g] - fh.adf[side][g]);
            }
        }
    }
}

#[test]
fn branch_densities_are_held_at_each_residual_point() {
    // core が L の重みの指標を作る、履歴の分岐の均質化した数密度 N̄ᵢ^j
    let (_, _, rod) = residuals(&nominal(), MODERATOR_DENSITY, ROD_OUT, &pair()).unwrap();
    assert_eq!(rod.densities, branch().micro.unwrap().densities);
}

#[test]
fn the_nominal_condition_with_the_correction_and_both_parts_is_the_branch() {
    // core が履歴の分岐の点（w = 1、ノードの数密度が分岐のもの、L の指標の比が 1）で返す定数は、分岐の格子計算の定数（ADR-0215）
    let n = nominal();
    let b = branch();
    let (_, _, rod) = residuals(&n, MODERATOR_DENSITY, ROD_OUT, &pair()).unwrap();
    let micro = n.micro.as_ref().unwrap();
    for (k, (&p, &u)) in [1usize, 4].iter().zip(&[0usize, 2]).enumerate() {
        let d: Vec<f64> = (0..2)
            .map(|i| b.micro.as_ref().unwrap().densities[k][i] - micro.densities[u][i])
            .collect();
        let corrected = linear(&n.states[p], &micro.states[p], &d);
        let mut r: Residual = rod.untracked[k].clone();
        let l = &rod.linearization[k];
        for g in 0..2 {
            r.total[g] += l.total[g];
            r.absorption[g] += l.absorption[g];
            r.nu_fission[g] += l.nu_fission[g];
            for h in 0..2 {
                r.scatter[g][h] += l.scatter[g][h];
            }
        }
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
        history_correction(&n, MODERATOR_DENSITY, &[(ROD_OUT.into(), pair())]).unwrap();
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
        // ハイブリッドも同じ条件（分岐の条件の検査を確かめる）
        let mut h = hybrid();
        change(&mut h);
        let err = residuals(
            &n,
            MODERATOR_DENSITY,
            ROD_OUT,
            &BranchTables {
                branch: b,
                hybrid: h,
            },
        )
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
    let mut h = hybrid();
    h.axes = b.axes.clone();
    let err = residuals(
        &n2,
        MODERATOR_DENSITY,
        ROD_OUT,
        &BranchTables {
            branch: b,
            hybrid: h,
        },
    )
    .unwrap_err()
    .to_string();
    assert!(err.contains("moderator"), "{err}");
    // 分岐のテーブルは状態 1 つ
    let err = fails(&|b| {
        b.axes[2] = Axis::moderator(vec![565.0, 600.0], vec![0.743, 0.661]).with_nominal(600.0)
    });
    assert!(err.contains("one"), "{err}");
    // 知らない変数
    let err = residuals(&n, "boron", ROD_OUT, &pair())
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
    let mut h = b.clone();
    h.states = vec![few(1.5), few(6.0), few(3.5)];
    // 公称のテーブルは条件の状態を 10 で解いていない
    let err = residuals(
        &n,
        MODERATOR_DENSITY,
        ROD_OUT,
        &BranchTables {
            branch: b,
            hybrid: h,
        },
    )
    .unwrap_err()
    .to_string();
    assert!(err.contains("burnup"), "{err}");
}

#[test]
fn branch_tables_hold_the_same_micro_nuclides_and_heavy_metal() {
    let n = nominal();
    let with = |change: &dyn Fn(&mut Table)| {
        let mut b = branch();
        change(&mut b);
        let mut h = hybrid();
        change(&mut h);
        residuals(
            &n,
            MODERATOR_DENSITY,
            ROD_OUT,
            &BranchTables {
                branch: b,
                hybrid: h,
            },
        )
    };
    let err = with(&|b| b.micro.as_mut().unwrap().nuclides[1] = "Pu-239".into())
        .unwrap_err()
        .to_string();
    assert!(err.contains("nuclides"), "{err}");
    let err = with(&|b| b.heavy_metal = Some(401.0))
        .unwrap_err()
        .to_string();
    assert!(err.contains("heavy metal"), "{err}");
    assert!(with(&|b| b.micro = None).is_err());
}

/// 反応率の微視的断面積の差 `a − b` の核分裂、捕獲、収率の境界より下の核分裂の断面積（割合 × 群の和の σf）
fn rate_parts(a: &MicroNuclide, b: &MicroNuclide) -> (Vec<f64>, Vec<f64>, f64) {
    let below = |x: &MicroNuclide| x.thermal_fission.unwrap() * x.fission.iter().sum::<f64>();
    (
        (0..2).map(|g| a.fission[g] - b.fission[g]).collect(),
        (0..2)
            .map(|g| a.reactions[0].channels[0][g] - b.reactions[0].channels[0][g])
            .collect(),
        below(a) - below(b),
    )
}

#[test]
fn rate_residuals_split_at_the_hybrid_like_the_constants() {
    // 追わない核種の分は分岐 − ハイブリッド、線形化の誤差はハイブリッド − 公称の条件 j の状態の、反応率の微視的断面積の差
    // （U と L と同じハイブリッドで分ける。ADR-0215）
    let n = nominal();
    let mut b = branch();
    let mut h = hybrid();
    // 収率の境界より下の割合も 3 つで違う
    for state in &mut b.micro.as_mut().unwrap().states {
        state[0].thermal_fission = Some(0.85);
    }
    for state in &mut h.micro.as_mut().unwrap().states {
        state[0].thermal_fission = Some(0.875);
    }
    let tables = BranchTables {
        branch: b.clone(),
        hybrid: h.clone(),
    };
    let (_, _, rod) = residuals(&n, MODERATOR_DENSITY, ROD_OUT, &tables).unwrap();
    let micro = n.micro.as_ref().unwrap();
    let (branch_micro, hybrid_micro) = (b.micro.as_ref().unwrap(), h.micro.as_ref().unwrap());
    assert_eq!(rod.micro.untracked.len(), 2);
    assert_eq!(rod.micro.linearization.len(), 2);
    for (k, &p) in [1usize, 4].iter().enumerate() {
        assert_eq!(rod.micro.untracked[k].len(), 2);
        assert_eq!(rod.micro.linearization[k].len(), 2);
        // 微分を持つ Pu-240 も、反応率の σᵢ の差（微分の差でない）
        for i in 0..2 {
            let (nj, bj, hj) = (
                &micro.states[p][i],
                &branch_micro.states[k][i],
                &hybrid_micro.states[k][i],
            );
            for (r, (a, c)) in [
                (&rod.micro.untracked[k][i], (bj, hj)),
                (&rod.micro.linearization[k][i], (hj, nj)),
            ] {
                let (fission, capture, below) = rate_parts(a, c);
                assert_eq!(r.fission, fission);
                assert_eq!(r.reactions[0].channels[0], capture);
                assert_eq!(r.reactions[0].mt, 102);
                assert!(close(r.thermal_fission.unwrap(), below));
            }
            // 2 つの和は分岐 − 公称の条件 j の状態の差（分けない差の検算）
            let (fission, capture, below) = rate_parts(bj, nj);
            let (u, l) = (&rod.micro.untracked[k][i], &rod.micro.linearization[k][i]);
            for g in 0..2 {
                assert!(near(
                    u.fission[g] + l.fission[g],
                    fission[g],
                    bj.fission[g].abs() + hj.fission[g].abs() + nj.fission[g].abs()
                ));
                assert!(near(
                    u.reactions[0].channels[0][g] + l.reactions[0].channels[0][g],
                    capture[g],
                    bj.reactions[0].channels[0][g].abs()
                        + hj.reactions[0].channels[0][g].abs()
                        + nj.reactions[0].channels[0][g].abs()
                ));
            }
            assert!(near(
                u.thermal_fission.unwrap() + l.thermal_fission.unwrap(),
                below,
                below.abs() + bj.fission.iter().sum::<f64>() + nj.fission.iter().sum::<f64>()
            ));
        }
    }
}

#[test]
fn branch_tables_with_other_reactions_or_fission_are_errors() {
    // 分岐かハイブリッドのどちらかの反応の並びか核分裂の有無が、公称のテーブルと違えばエラー
    let n = nominal();
    let changes: Vec<fn(&mut MicroNuclide)> = vec![
        |x| x.reactions[0].mt = 103,
        |x| x.reactions[0].channels.push(vec![0.0; 2]),
        |x| x.thermal_fission = None,
    ];
    let change_table = |t: &mut Table, change: fn(&mut MicroNuclide)| {
        for state in &mut t.micro.as_mut().unwrap().states {
            change(&mut state[0]);
            if state[0].thermal_fission.is_none() {
                state[0].fission = vec![0.0; 2];
                state[0].nu_fission = vec![0.0; 2];
            }
        }
        t.validate().unwrap();
    };
    for change in changes {
        for in_hybrid in [false, true] {
            let mut tables = pair();
            if in_hybrid {
                change_table(&mut tables.hybrid, change);
            } else {
                change_table(&mut tables.branch, change);
            }
            let err = format!(
                "{:#}",
                residuals(&n, MODERATOR_DENSITY, ROD_OUT, &tables).unwrap_err()
            );
            assert!(err.contains("U-235"), "{in_hybrid}: {err}");
        }
    }
}

#[test]
fn rod_values_of_tables_without_a_rod_axis_are_out() {
    let err = residuals(&nominal(), MODERATOR_DENSITY, "aic", &pair())
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
branches = { out = { branch = "b.msgpack", hybrid = "h.msgpack" } }
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
                branches: vec![(
                    ROD_OUT.into(),
                    BranchFiles {
                        branch: "b.msgpack".into(),
                        hybrid: "h.msgpack".into(),
                    }
                )],
            }],
        }
    );
    for bad in [
        // 履歴がない
        "nominal = \"n\"\noutput = \"o\"\n",
        // 分岐のテーブルがない
        "nominal = \"n\"\noutput = \"o\"\n[[history]]\nvariable = \"moderator-density\"\nbranches = {}\n",
        // 知らない項目
        "nominal = \"n\"\noutput = \"o\"\nextra = 1\n[[history]]\nvariable = \"moderator-density\"\nbranches = { out = { branch = \"b\", hybrid = \"h\" } }\n",
        // 出力が公称のテーブル（読んだファイルを上書きしない）
        "nominal = \"n\"\noutput = \"n\"\n[[history]]\nvariable = \"moderator-density\"\nbranches = { out = { branch = \"b\", hybrid = \"h\" } }\n",
        // ハイブリッドのテーブルがない（U と L に分けられない）
        "nominal = \"n\"\noutput = \"o\"\n[[history]]\nvariable = \"moderator-density\"\nbranches = { out = { branch = \"b\" } }\n",
        // 分岐のテーブルだけの形（R を 1 つで持つ前の形）
        "nominal = \"n\"\noutput = \"o\"\n[[history]]\nvariable = \"moderator-density\"\nbranches = { out = \"b\" }\n",
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
    std::fs::write(dir.join("h.msgpack"), hybrid().to_msgpack().unwrap()).unwrap();
    let input = dir.join("histories.toml");
    std::fs::write(
        &input,
        "nominal = \"n.msgpack\"\noutput = \"o.msgpack\"\n[[history]]\nvariable = \"moderator-density\"\n\
         branches = { out = { branch = \"b.msgpack\", hybrid = \"h.msgpack\" } }\n",
    )
    .unwrap();
    run(&input).unwrap();
    let written = Table::from_msgpack(&std::fs::read(dir.join("o.msgpack")).unwrap()).unwrap();
    let mut expected = nominal();
    expected.histories = vec![
        history_correction(&expected, MODERATOR_DENSITY, &[(ROD_OUT.into(), pair())]).unwrap(),
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
    let mut h = hybrid();
    h.axes = b.axes.clone();
    let err = residuals(
        &n,
        MODERATOR_DENSITY,
        ROD_OUT,
        &BranchTables {
            branch: b,
            hybrid: h,
        },
    )
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
    let mut h = hybrid();
    h.axes = b.axes.clone();
    let err = residuals(
        &n,
        MODERATOR_DENSITY,
        ROD_OUT,
        &BranchTables {
            branch: b,
            hybrid: h,
        },
    )
    .unwrap_err()
    .to_string();
    assert!(err.contains("density"), "{err}");
}

#[test]
fn run_does_not_overwrite_the_tables_it_reads() {
    let dir = std::env::temp_dir().join(format!("xlc-residual-same-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let nominal_bytes = nominal().to_msgpack().unwrap();
    let branch_bytes = branch().to_msgpack().unwrap();
    std::fs::write(dir.join("n.msgpack"), &nominal_bytes).unwrap();
    std::fs::write(dir.join("b.msgpack"), &branch_bytes).unwrap();
    let hybrid_bytes = hybrid().to_msgpack().unwrap();
    std::fs::write(dir.join("h.msgpack"), &hybrid_bytes).unwrap();
    let input = dir.join("histories.toml");
    // 同じファイルを別の書き方で指す出力（公称のテーブル、分岐のテーブル、ハイブリッドのテーブル）
    for output in ["./n.msgpack", "b.msgpack", "./h.msgpack"] {
        std::fs::write(
            &input,
            format!(
                "nominal = \"n.msgpack\"\noutput = \"{output}\"\n[[history]]\nvariable = \"moderator-density\"\n\
                 branches = {{ out = {{ branch = \"b.msgpack\", hybrid = \"h.msgpack\" }} }}\n"
            ),
        )
        .unwrap();
        let err = run(&input).unwrap_err().to_string();
        assert!(err.contains("output"), "{err}");
        assert_eq!(std::fs::read(dir.join("n.msgpack")).unwrap(), nominal_bytes);
        assert_eq!(std::fs::read(dir.join("b.msgpack")).unwrap(), branch_bytes);
        assert_eq!(std::fs::read(dir.join("h.msgpack")).unwrap(), hybrid_bytes);
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

/// ハイブリッドのテーブルを変える関数と、エラーの文に入る語
type HybridChange = (fn(&mut Table), &'static str);

#[test]
fn hybrids_solve_the_branch_points_with_the_branch_tracked_densities() {
    // U と L は同じ条件と燃焼度の点の 2 つの定数の差で、ハイブリッドの追う核種は分岐と同じ数密度（ADR-0215）
    let n = nominal();
    let changes: Vec<HybridChange> = vec![
        (
            |h| h.axes[2] = Axis::moderator(vec![565.0], vec![0.743]).with_nominal(565.0),
            "hybrid",
        ),
        (
            |h| {
                h.axes[0] = Axis::Burnup {
                    points: vec![0.0, 10.0],
                }
            },
            "hybrid",
        ),
        (
            |h| h.micro.as_mut().unwrap().densities[1][0] = 5.3e-4,
            "densities",
        ),
        (
            |h| h.micro.as_mut().unwrap().nuclides[1] = "Pu-239".into(),
            "nuclides",
        ),
        (|h| h.heavy_metal = Some(401.0), "heavy metal"),
        (|h| h.micro = None, "micro"),
    ];
    for (change, expected) in changes {
        let mut h = hybrid();
        change(&mut h);
        let tables = BranchTables {
            branch: branch(),
            hybrid: h,
        };
        let err = format!(
            "{:#}",
            residuals(&n, MODERATOR_DENSITY, ROD_OUT, &tables).unwrap_err()
        );
        assert!(err.contains(expected), "{expected:?} not in {err}");
    }
}
