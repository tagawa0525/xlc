use super::*;

#[test]
fn thermal_correction_is_row_sum_minus_elastic() {
    let elastic = vec![4.0, 20.0];
    let matrix = vec![vec![3.9, 0.1], vec![2.0, 27.0]];
    assert_eq!(thermal_correction(&elastic, &matrix), vec![0.0, 9.0]);
}

/// 非分離共鳴（LSSF=1、s 波の1系列）と、MF3 の一定の断面積だけのテープ。
fn unresolved_tape() -> Tape {
    use crate::endf::rec;
    let mf3 = |mt: u32, value: &str| {
        [
            rec(["1.001000+3", "2.360060+2", "0", "0", "0", "0"], 1, 3, mt),
            rec(["0.0", "0.0", "0", "0", "1", "2"], 1, 3, mt),
            rec(["2", "2", "", "", "", ""], 1, 3, mt),
            rec(["1.0-5", value, "2.0+7", value, "", ""], 1, 3, mt),
            rec(["", "", "", "", "", ""], 1, 3, 0),
        ]
        .concat()
    };
    let text = [
        rec(["1.001000+3", "2.360060+2", "0", "0", "1", "0"], 1, 2, 151),
        rec(["1.001000+3", "1.0", "0", "0", "1", "0"], 1, 2, 151),
        rec(["1.0+3", "1.0+5", "2", "2", "0", "1"], 1, 2, 151),
        rec(["0.0", "0.948", "1", "0", "1", "0"], 1, 2, 151),
        rec(["2.360060+2", "0.0", "0", "0", "1", "0"], 1, 2, 151),
        rec(["0.5", "0.0", "2", "0", "18", "2"], 1, 2, 151),
        rec(["0.0", "0.0", "0.0", "1.0", "0.0", "0.0"], 1, 2, 151),
        rec(["1.0+3", "20.0", "0.0", "1.0-3", "2.3-2", "0.0"], 1, 2, 151),
        rec(["1.0+5", "20.0", "0.0", "1.0-3", "2.3-2", "0.0"], 1, 2, 151),
        rec(["", "", "", "", "", ""], 1, 2, 0),
        mf3(1, "11.0"),
        mf3(2, "10.0"),
        mf3(102, "1.0"),
    ]
    .concat();
    Tape::parse(&text).unwrap()
}

/// 分離共鳴（MLBW、1e-5 eV〜20 keV、6.67 eV に1つ）を持つテープ
fn resolved_tape() -> Tape {
    use crate::endf::rec;
    let text = [
        rec(["9.223800+4", "2.360060+2", "0", "0", "1", "0"], 1, 2, 151),
        rec(["9.223800+4", "1.0", "0", "0", "1", "0"], 1, 2, 151),
        rec(["1.0-5", "2.0+4", "1", "2", "0", "0"], 1, 2, 151),
        rec(["0.0", "0.948", "0", "0", "1", "0"], 1, 2, 151),
        rec(["2.360060+2", "0.0", "0", "0", "6", "1"], 1, 2, 151),
        rec(
            ["6.674", "0.5", "2.6e-2", "1.5e-3", "2.3e-2", "0.0"],
            1,
            2,
            151,
        ),
        rec(["", "", "", "", "", ""], 1, 2, 0),
    ]
    .concat();
    Tape::parse(&text).unwrap()
}

#[test]
fn slowing_down_covers_the_resolved_range_in_whole_groups() {
    let groups = GroupStructure::new(vec![1.0e5, 3.0e4, 1.0e4, 10.0, 1.0e-5]).unwrap();
    let slowing = slowing_down_range(&resolved_tape(), 1, &groups, &Weight::lwr())
        .unwrap()
        .unwrap();
    assert_eq!(slowing.awr, 236.006);
    // 範囲の上の散乱源に使うポテンシャル散乱は、分離共鳴の範囲の s 波の散乱半径から
    let potential = resolved(&resolved_tape(), 1).unwrap().unwrap().potential;
    assert_eq!(slowing.potential, potential);
    // 重み関数の熱中性子の区切り 0.1 eV 以上で最も低い群の境界から、分離共鳴の上端 20 keV 以下で
    // 最も高い群の境界まで
    assert_eq!(slowing.range, (10.0, 1.0e4));
}

#[test]
fn resolved_range_has_its_top_and_potential_scattering() {
    use crate::endf::rec;
    let expected = 4.0 * std::f64::consts::PI * 0.948 * 0.948;
    // MLBW: AP
    let mlbw = resolved(&resolved_tape(), 1).unwrap().unwrap();
    assert_eq!(mlbw.top, 2.0e4);
    assert!((mlbw.potential - expected).abs() < 1e-12 * expected);
    // Reich-Moore で AP が 0 なら、l=0 の APL（U-238 と同じ）
    let text = [
        rec(["9.223800+4", "2.360060+2", "0", "0", "1", "0"], 1, 2, 151),
        rec(["9.223800+4", "1.0", "0", "0", "1", "0"], 1, 2, 151),
        rec(["1.0-5", "2.0+4", "1", "3", "0", "0"], 1, 2, 151),
        rec(["0.0", "0.0", "0", "0", "1", "0"], 1, 2, 151),
        rec(["2.360060+2", "0.948", "0", "0", "6", "1"], 1, 2, 151),
        rec(
            ["6.674", "0.5", "1.5e-3", "2.3e-2", "0.0", "0.0"],
            1,
            2,
            151,
        ),
        rec(["", "", "", "", "", ""], 1, 2, 0),
    ]
    .concat();
    let rm = resolved(&Tape::parse(&text).unwrap(), 1).unwrap().unwrap();
    assert!((rm.potential - expected).abs() < 1e-12 * expected);
    // 分離共鳴がなければ持たない
    assert!(resolved(&unresolved_tape(), 1).unwrap().is_none());
}

#[test]
fn rejects_a_zero_s_wave_radius() {
    use crate::endf::{ReichMoore, RmLValue};
    let rm = RangeKind::ReichMoore(ReichMoore {
        spi: 0.0,
        ap: 0.0,
        naps: 0,
        l_values: vec![RmLValue {
            awri: 236.0,
            apl: 0.0,
            l: OrbitalL::S,
            resonances: vec![],
        }],
    });
    assert!(matches!(
        s_wave_radius(&rm),
        Err(ProcessError::SWaveRadius(_))
    ));
}

#[test]
fn rejects_differing_s_wave_radii() {
    use crate::endf::{ParticlePair, RMatrixLimited, RmlChannel, SpinGroup};
    let pair = ParticlePair {
        ma: 1.0,
        mb: 56.0,
        za: 0.0,
        zb: 26.0,
        ia: 0.5,
        ib: 0.5,
        q: 0.0,
        pnt: 1,
        shf: 0,
        mt: 2,
    };
    let group = |aj: f64, ape: f64| SpinGroup {
        aj,
        pj: 0.0,
        channels: vec![RmlChannel {
            pair: 0,
            l: OrbitalL::S,
            sch: aj,
            bnd: 0.0,
            ape,
            apt: ape,
        }],
        resonances: vec![],
    };
    let rml = |second: f64| {
        RangeKind::RMatrixLimited(RMatrixLimited {
            pairs: vec![pair.clone()],
            spin_groups: vec![group(1.0, 0.63), group(2.0, second)],
        })
    };
    assert_eq!(s_wave_radius(&rml(0.63)).unwrap(), 0.63);
    assert!(matches!(
        s_wave_radius(&rml(0.70)),
        Err(ProcessError::SWaveRadii(_))
    ));
}

#[test]
fn pointwise_data_covers_the_range_on_the_total_grid() {
    // 全断面積は部分反応の和なので、その格子は全反応の格子点と不連続点（x の重複）を含む
    let total = Tab1::linear(
        vec![1.0, 5.0, 10.0, 10.0, 20.0, 40.0],
        vec![10.0, 12.0, 14.0, 24.0, 20.0, 18.0],
    );
    let elastic = Tab1::linear(vec![1.0, 40.0], vec![8.0, 8.0]);
    let fission = Tab1::linear(vec![1.0, 40.0], vec![1.0, 2.0]);
    let nubar = Tab1::linear(vec![1.0, 40.0], vec![2.5, 2.6]);
    let tables = BTreeMap::from([(1, total), (2, elastic), (18, fission)]);
    let data = pointwise_within(&tables, Some(&nubar), (4.0, 30.0)).unwrap();
    // 範囲の端と、その間の格子点（不連続点は左右の極限を両方持つ）
    assert_eq!(data.energy, vec![4.0, 5.0, 10.0, 10.0, 20.0, 30.0]);
    assert_eq!(data.total[2], 14.0);
    assert_eq!(data.total[3], 24.0);
    assert!((data.total[0] - 11.5).abs() < 1e-12);
    assert!(data.elastic.iter().all(|&v| v == 8.0));
    let f = |e: f64| 1.0 + (e - 1.0) / 39.0;
    let nu = |e: f64| 2.5 + 0.1 * (e - 1.0) / 39.0;
    for (i, &e) in data.energy.iter().enumerate() {
        assert!((data.fission[i] - f(e)).abs() < 1e-12);
        assert!((data.nu_fission[i] - nu(e) * f(e)).abs() < 1e-12);
    }
}

#[test]
fn pointwise_data_without_fission_is_zero() {
    let tables = BTreeMap::from([
        (1, Tab1::linear(vec![1.0, 40.0], vec![10.0, 10.0])),
        (2, Tab1::linear(vec![1.0, 40.0], vec![8.0, 8.0])),
    ]);
    let data = pointwise_within(&tables, None, (4.0, 30.0)).unwrap();
    assert!(
        data.fission
            .iter()
            .chain(&data.nu_fission)
            .all(|&v| v == 0.0)
    );
}

#[test]
fn no_resolved_range_needs_no_slowing_down() {
    let groups = GroupStructure::new(vec![1.0e5, 1.0e4, 1.0e-5]).unwrap();
    assert!(
        slowing_down_range(&unresolved_tape(), 1, &groups, &Weight::lwr())
            .unwrap()
            .is_none()
    );
}

#[test]
fn unresolved_shielding_lowers_capture_at_small_background() {
    use crate::multigroup::shielded_averages;
    let tape = unresolved_tape();
    let sigma0 = [1.0e10, 1.0];
    let shielding = unresolved_shielding(&tape, 1, 300.0, &sigma0)
        .unwrap()
        .unwrap();
    let tables = crate::pointwise::reconstruct(&tape, 1, 1.0e-3).unwrap();
    let groups = GroupStructure::new(vec![5.0e4, 1.0e4]).unwrap();
    let shielded = shielded_averages(
        &groups,
        &Weight::lwr(),
        &tables,
        None,
        &sigma0,
        Some(&shielding),
        None,
    )
    .unwrap();
    let plain =
        shielded_averages(&groups, &Weight::lwr(), &tables, None, &sigma0, None, None).unwrap();
    // 無限希釈では遮蔽因子は 1、σ0 = 1 b では捕獲が減る
    let capture = |avg: &crate::multigroup::ReactionAverages| avg.sigma[&102][0];
    assert!((capture(&shielded[0]) / capture(&plain[0]) - 1.0).abs() < 1e-5);
    assert!(capture(&shielded[1]) < 0.99 * capture(&plain[1]));
}

#[test]
fn no_unresolved_range_gives_no_shielding() {
    use crate::endf::rec;
    let text = [
        rec(["1.001000+3", "1.0", "0", "0", "0", "0"], 1, 3, 1),
        rec(["0.0", "0.0", "0", "0", "1", "2"], 1, 3, 1),
        rec(["2", "2", "", "", "", ""], 1, 3, 1),
        rec(["1.0-5", "1.0", "2.0+7", "1.0", "", ""], 1, 3, 1),
        rec(["", "", "", "", "", ""], 1, 3, 0),
    ]
    .concat();
    let tape = Tape::parse(&text).unwrap();
    assert_eq!(
        unresolved_shielding(&tape, 1, 300.0, &[10.0]).unwrap(),
        None
    );
}

/// MF7 の MT4（定数と β が1点の小さな表）と、指定があれば MT2 を持つ熱中性子散乱のテープ。
fn thermal_tape(with_elastic: bool) -> Tape {
    use crate::endf::rec;
    let mut lines: Vec<String> = vec![
        rec(["1.001000+3", "9.991673-1", "0", "0", "0", "0"], 1, 7, 4),
        rec(["0.0", "0.0", "0", "0", "6", "0"], 1, 7, 4),
        rec(["4.0", "10.0", "9.9917-1", "5.0", "0.0", "2.0"], 1, 7, 4),
        rec(["0.0", "0.0", "0", "0", "1", "2"], 1, 7, 4),
        rec(["2", "4", "", "", "", ""], 1, 7, 4),
    ];
    for beta in ["0.0", "2.0"] {
        lines.push(rec(["2.936000+2", beta, "0", "0", "1", "2"], 1, 7, 4));
        lines.push(rec(["2", "4", "", "", "", ""], 1, 7, 4));
        lines.push(rec(["0.1", "1.0", "10.0", "0.01", "", ""], 1, 7, 4));
    }
    lines.extend([
        rec(["0.0", "0.0", "0", "0", "1", "1"], 1, 7, 4),
        rec(["1", "2", "", "", "", ""], 1, 7, 4),
        rec(["2.936000+2", "1400.0", "", "", "", ""], 1, 7, 4),
        rec(["", "", "", "", "", ""], 1, 7, 0),
    ]);
    if with_elastic {
        // 非干渉性弾性散乱（LTHR=2）: SB と Debye-Waller 積分
        lines.splice(
            0..0,
            [
                rec(["1.001000+3", "9.991673-1", "2", "0", "0", "0"], 1, 7, 2),
                rec(["2.0", "0.0", "0", "0", "1", "1"], 1, 7, 2),
                rec(["1", "2", "", "", "", ""], 1, 7, 2),
                rec(["2.936000+2", "10.0", "", "", "", ""], 1, 7, 2),
                rec(["", "", "", "", "", ""], 1, 7, 0),
            ],
        );
    }
    Tape::parse(&lines.concat()).unwrap()
}

#[test]
fn reads_incoherent_inelastic_scattering() {
    let data = read_thermal_scattering(&thermal_tape(false)).unwrap();
    assert_eq!(data.beta, vec![0.0, 2.0]);
}

#[test]
fn sab_limits_follow_the_input() {
    // テープの MF7 MT4 の上限 B(4) は 5 eV
    let data = read_thermal_scattering(&thermal_tape(false)).unwrap();
    let thermal = |upper, free_gas_upper| ThermalScattering {
        data: data.clone(),
        upper,
        outside: SabOutside::ShortCollisionTime,
        free_gas_upper,
    };
    assert_eq!(
        thermal(SabUpper::Evaluated, 1.0e3).limits().unwrap(),
        (5.0, 1.0e3)
    );
    assert_eq!(
        thermal(SabUpper::Energy(4.0), 4.0).limits().unwrap(),
        (4.0, 4.0)
    );
    // 評価済みデータの上限より上では S(α, β) を使えない
    for upper in [5.5, 0.0, -1.0, f64::NAN] {
        assert!(
            matches!(
                thermal(SabUpper::Energy(upper), 1.0e3).limits(),
                Err(ProcessError::SabUpper { .. })
            ),
            "{upper}"
        );
    }
    // 自由ガスは S(α, β) の上限から続ける
    for free_gas_upper in [4.9, f64::INFINITY, f64::NAN] {
        assert!(
            matches!(
                thermal(SabUpper::Evaluated, free_gas_upper).limits(),
                Err(ProcessError::FreeGasUpper { .. })
            ),
            "{free_gas_upper}"
        );
    }
}

#[test]
fn thermal_elastic_scattering_is_rejected() {
    assert!(read_thermal_scattering(&thermal_tape(true)).is_err());
}

#[test]
fn scattering_order_sets_the_number_of_moments() {
    let tape = unresolved_tape();
    for order in [0, 2] {
        let settings = Settings {
            groups: GroupStructure::wims69(),
            weight: Weight::lwr(),
            sigma0: vec![1.0e10],
            thermal_cutoff: 4.0,
            tolerance: 1.0e-3,
            pointwise_range: Some((4.0, 1.0e4)),
            scattering_order: order,
        };
        let data = process_nuclide(&tape, "X", &[293.6], None, &settings).unwrap();
        assert_eq!(data.temperatures[0].scatter.len(), order + 1);
    }
}

#[test]
fn processes_each_temperature_and_sigma0() {
    let tape = unresolved_tape();
    let settings = Settings {
        groups: GroupStructure::wims69(),
        weight: Weight::lwr(),
        sigma0: vec![1.0e10, 1.0],
        thermal_cutoff: 4.0,
        tolerance: 1.0e-3,
        pointwise_range: Some((4.0, 1.0e4)),
        scattering_order: 0,
    };
    let data = process_nuclide(&tape, "X", &[293.6, 600.0], None, &settings).unwrap();
    assert_eq!(data.name, "X");
    let temperatures: Vec<f64> = data.temperatures.iter().map(|t| t.temperature).collect();
    assert_eq!(temperatures, vec![293.6, 600.0]);
    for t in &data.temperatures {
        assert_eq!(t.shielded.len(), 2);
        let [dilute, shielded] = &t.shielded[..] else {
            panic!("two σ0 values")
        };
        let last = dilute.total.len() - 1;
        // 熱群の補正は全断面積と弾性散乱の両方に入るので、差は同じ温度に拡げた捕獲の群平均に等しい
        // （一定の 1 b でも、拡げると低いエネルギーで大きくなる）
        let tables = crate::pointwise::reconstruct_at(&tape, 1, 1.0e-3, t.temperature).unwrap();
        let capture =
            crate::multigroup::reaction_averages(&settings.groups, &settings.weight, &tables, None)
                .sigma[&102][last];
        assert!(capture > 1.0);
        let absorption = dilute.total[last] - dilute.elastic[last];
        // σ0 = 1e10 の重み 1/(σt + σ0) と無限希釈の重みは、σt/σ0 ≈ 1e-9 だけ違う
        assert!(
            (absorption - capture).abs() <= 1e-9 * capture,
            "{absorption} vs {capture}"
        );
        // 自由ガスの行の和は温度に拡げた弾性散乱の群平均になり、熱領域の補正は 0 になる
        let elastic =
            crate::multigroup::reaction_averages(&settings.groups, &settings.weight, &tables, None)
                .sigma[&2][last];
        assert!(
            (dilute.elastic[last] - elastic).abs() <= 1e-9 * elastic,
            "{} vs {elastic}",
            dilute.elastic[last]
        );
        // 一定の 10 b を拡げると 10 N(E)（N(E) = (1 + 1/(2a²)) erf(a) + e^(-a²)/(a√π)、a² = A E/kT）。
        // 許容誤差はポイントワイズの断面積の許容誤差 1e-3
        let (awr, kt) = (236.006, BK * t.temperature);
        let effective = |e: f64| {
            let a = (awr * e / kt).sqrt();
            let n = (1.0 + 0.5 / (a * a)) * crate::special::erf(a)
                + (-a * a).exp() / (a * std::f64::consts::PI.sqrt());
            10.0 * n
        };
        let integrate = |f: &dyn Fn(f64) -> f64| {
            crate::multigroup::group_integrals(&settings.groups, &settings.weight, &[], f)[last]
        };
        let expected = integrate(&effective) / integrate(&|_| 1.0);
        assert!(
            (elastic - expected).abs() <= 1e-3 * expected,
            "{elastic} vs {expected}"
        );
        // 非分離共鳴の群（[1.5e4, 2.5e4] eV）では、σ0 が小さいと吸収が自己遮蔽で減る
        let g = settings.groups.find(2.0e4).unwrap();
        let absorption = |s: &Shielded| s.total[g] - s.elastic[g];
        assert!(absorption(shielded) < absorption(dilute));
        assert!(t.chi.iter().all(|&c| c == 0.0));
    }
}
