use super::*;

fn constant(value: f64) -> Tab1 {
    Tab1::linear(vec![1.0e-5, 2.0e7], vec![value, value])
}

#[test]
fn constant_cross_section_averages_to_itself() {
    let groups = GroupStructure::wims69();
    let weight = Weight::lwr();
    let tables = BTreeMap::from([(1, constant(3.0)), (18, constant(0.5))]);
    let nu = constant(2.5);
    let avg = reaction_averages(&groups, &weight, &tables, Some(&nu));
    assert_eq!(avg.flux.len(), 69);
    for g in 0..69 {
        assert!((avg.sigma[&1][g] - 3.0).abs() < 1e-12);
        assert!((avg.nu_fission.as_ref().unwrap()[g] - 1.25).abs() < 1e-12);
    }
}

#[test]
fn no_fission_gives_no_nu_fission() {
    let groups = GroupStructure::wims69();
    let tables = BTreeMap::from([(2, constant(3.0))]);
    let avg = reaction_averages(&groups, &Weight::lwr(), &tables, None);
    assert!(avg.nu_fission.is_none());
}

#[test]
fn threshold_reaction_is_zero_below_threshold() {
    let groups = GroupStructure::new(vec![2.0e6, 1.0e6, 5.0e5]).unwrap();
    let tables = BTreeMap::from([(16, Tab1::linear(vec![1.5e6, 2.0e6], vec![1.0, 1.0]))]);
    let avg = reaction_averages(&groups, &Weight::lwr(), &tables, None);
    assert_eq!(avg.sigma[&16][1], 0.0);
    assert!(avg.sigma[&16][0] > 0.0 && avg.sigma[&16][0] < 1.0);
}
#[test]
fn constant_total_gives_infinite_dilution_for_every_sigma0() {
    let groups = GroupStructure::wims69();
    let weight = Weight::lwr();
    let capture = Tab1::linear(vec![1.0e-5, 1.0, 2.0e7], vec![5.0, 1.0, 0.1]);
    let tables = BTreeMap::from([(1, constant(3.0)), (102, capture)]);
    let infinite = reaction_averages(&groups, &weight, &tables, None);
    let shielded = shielded_averages(
        &groups,
        &weight,
        &tables,
        None,
        &[1.0e10, 10.0, 1.0],
        None,
        None,
    )
    .unwrap();
    assert_eq!(shielded.len(), 3);
    for avg in &shielded {
        for g in 0..69 {
            let (a, b) = (avg.sigma[&102][g], infinite.sigma[&102][g]);
            assert!((a - b).abs() <= 1e-12 * b, "group {g}: {a} vs {b}");
        }
    }
}

#[test]
fn narrow_resonance_flux_matches_analytic_average() {
    // 1/E の重みの群で σx = σt = cE のとき
    //   σ_g = ln((c b + σ0)/(c a + σ0)) / [(1/σ0) ln(b (c a + σ0) / (a (c b + σ0)))]
    let groups = GroupStructure::new(vec![100.0, 10.0]).unwrap();
    let c = 2.0;
    let total = Tab1::linear(vec![10.0, 100.0], vec![10.0 * c, 100.0 * c]);
    let tables = BTreeMap::from([(1, total)]);
    let (a, b) = (10.0_f64, 100.0_f64);
    for s0 in [1.0, 30.0, 1.0e3] {
        let avg = &shielded_averages(&groups, &Weight::lwr(), &tables, None, &[s0], None, None)
            .unwrap()[0];
        let numerator = ((c * b + s0) / (c * a + s0)).ln();
        let denominator = (b * (c * a + s0) / (a * (c * b + s0))).ln() / s0;
        let expected = numerator / denominator;
        let actual = avg.sigma[&1][0];
        assert!(
            (actual - expected).abs() <= 1e-10 * expected,
            "σ0 {s0}: {actual} vs {expected}"
        );
    }
}

/// 20 eV に散乱の大きな共鳴（幅 1 eV）、ほかは散乱 10 b の核種の全断面積・弾性散乱・捕獲
fn scattering_resonance() -> BTreeMap<u32, Tab1> {
    let x: Vec<f64> = (0..=40_000)
        .map(|i| 1.0e-5 * (2.0e7_f64 / 1.0e-5).powf(f64::from(i) / 40_000.0))
        .collect();
    let peak = |e: f64| 1.0 / (1.0 + ((e - 20.0) / 0.5).powi(2));
    let table = |f: &dyn Fn(f64) -> f64| Tab1::linear(x.clone(), x.iter().map(|&e| f(e)).collect());
    BTreeMap::from([
        (1, table(&|e| 10.0 + 3000.0 * peak(e))),
        (2, table(&|e| 10.0 + 2000.0 * peak(e))),
        (102, table(&|e| 1000.0 * peak(e))),
    ])
}

#[test]
fn slowing_down_weight_changes_only_groups_in_its_range() {
    let groups = GroupStructure::new(vec![1.0e3, 100.0, 30.0, 10.0, 1.0]).unwrap();
    let tables = scattering_resonance();
    let slowing = SlowingDown {
        awr: 1.0e4,
        potential: 10.0,
        range: (1.0, 100.0),
    };
    let average = |s: Option<&SlowingDown>| {
        shielded_averages(&groups, &Weight::lwr(), &tables, None, &[50.0], None, s)
            .unwrap()
            .remove(0)
    };
    let (nr, ufg) = (average(None), average(Some(&slowing)));
    // 範囲の上の群は変わらない
    assert!((ufg.sigma[&102][0] - nr.sigma[&102][0]).abs() <= 1e-12 * nr.sigma[&102][0]);
    // 共鳴の群: NR 近似は吸収体の散乱で中性子束を減らしすぎるので、捕獲を小さく見積もる
    assert!(
        ufg.sigma[&102][2] > 1.05 * nr.sigma[&102][2],
        "{} vs {}",
        ufg.sigma[&102][2],
        nr.sigma[&102][2]
    );
}

#[test]
fn flux_is_continuous_at_the_top_of_the_slowing_down_range() {
    // 範囲の上の Bondarenko の中性子束には (σp + σ0) を掛けて上端でつなぐ。σs(E_high) = 10 b と
    // σp = 4 b が違うので、σs(E_high) + σ0 を掛けると上端で 60/54 倍の段差になる
    let (total, elastic) = (constant(10.0), constant(10.0));
    let slowing = SlowingDown {
        awr: 236.0,
        potential: 4.0,
        range: (1.0, 100.0),
    };
    let s0 = 50.0;
    let shape = slowing_down_shape(&total, &elastic, &slowing, &Weight::lwr(), s0).unwrap();
    let flux = slowing_down_flux(&shape, &slowing, s0, |_| 1.0 / (10.0 + s0)).unwrap();
    let high = slowing.range.1;
    let (inside, above) = (flux(high), flux(high.next_up()));
    assert!(
        (inside - above).abs() < 1e-12 * inside,
        "{inside} vs {above}"
    );
}

#[test]
fn below_the_slowing_down_range_the_flux_is_the_weight() {
    // 範囲の下（熱中性子の領域）では σ0 によらず無限希釈の値
    let groups = GroupStructure::new(vec![1.0e3, 100.0, 1.0, 0.01]).unwrap();
    let capture = Tab1::linear(vec![1.0e-5, 2.0e7], vec![300.0, 0.0]);
    let tables = BTreeMap::from([(1, constant(20.0)), (2, constant(10.0)), (102, capture)]);
    let slowing = SlowingDown {
        awr: 236.0,
        potential: 10.0,
        range: (1.0, 100.0),
    };
    let infinite = reaction_averages(&groups, &Weight::lwr(), &tables, None);
    let avg = &shielded_averages(
        &groups,
        &Weight::lwr(),
        &tables,
        None,
        &[1.0],
        None,
        Some(&slowing),
    )
    .unwrap()[0];
    let (a, b) = (avg.sigma[&102][2], infinite.sigma[&102][2]);
    assert!((a - b).abs() < 1e-12 * b, "{a} vs {b}");
}

#[test]
fn slowing_down_weight_keeps_constant_cross_sections() {
    let groups = GroupStructure::new(vec![1.0e3, 100.0, 10.0, 1.0]).unwrap();
    let tables = BTreeMap::from([
        (1, constant(10.0)),
        (2, constant(8.0)),
        (102, constant(2.0)),
    ]);
    let slowing = SlowingDown {
        awr: 236.0,
        potential: 8.0,
        range: (1.0, 100.0),
    };
    let avg = &shielded_averages(
        &groups,
        &Weight::lwr(),
        &tables,
        None,
        &[50.0],
        None,
        Some(&slowing),
    )
    .unwrap()[0];
    for g in 0..3 {
        assert!(
            (avg.sigma[&102][g] - 2.0).abs() < 1e-12,
            "{}",
            avg.sigma[&102][g]
        );
    }
}

#[test]
fn sigma0_must_be_positive_and_finite() {
    let tables = BTreeMap::from([(1, constant(3.0))]);
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let result = shielded_averages(
            &GroupStructure::wims69(),
            &Weight::lwr(),
            &tables,
            None,
            &[bad],
            None,
            None,
        );
        assert!(result.is_err(), "σ0 = {bad}");
    }
}

#[test]
fn shielding_needs_total_cross_section() {
    let tables = BTreeMap::from([(2, constant(3.0))]);
    let result = shielded_averages(
        &GroupStructure::wims69(),
        &Weight::lwr(),
        &tables,
        None,
        &[10.0],
        None,
        None,
    );
    assert!(result.is_err());
}

/// 1群 [1e4, 1e5] eV を覆う非分離共鳴の補正。σ0 は1つで、遮蔽した値と無限希釈の値は一定
fn uniform_shielding(lssf: i64) -> UnresolvedShielding {
    UnresolvedShielding::new(
        1.0e3,
        1.0e6,
        lssf,
        vec![1.0e3, 1.0e6],
        Interpolation::LinLin,
        vec![vec![[8.0, 7.5, 0.0, 0.5]; 2]],
        vec![[10.0, 9.0, 0.0, 1.0]; 2],
    )
    .unwrap()
}

fn one_group_tables() -> BTreeMap<u32, Tab1> {
    BTreeMap::from([
        (1, constant(11.0)),
        (2, constant(9.0)),
        (102, constant(2.0)),
    ])
}

#[test]
fn lssf1_scales_by_shielding_factor() {
    let groups = GroupStructure::new(vec![1.0e5, 1.0e4]).unwrap();
    let shielding = uniform_shielding(1);
    let avg = &shielded_averages(
        &groups,
        &Weight::lwr(),
        &one_group_tables(),
        None,
        &[10.0],
        Some(&shielding),
        None,
    )
    .unwrap()[0];
    // 捕獲 2 × 0.5/1、全断面積 11 × 8/10（中性子束も補正した全断面積で決まるが、一定なので平均は変わらない）
    assert!((avg.sigma[&102][0] - 1.0).abs() < 1e-12);
    assert!((avg.sigma[&1][0] - 8.8).abs() < 1e-12);
}

#[test]
fn lssf0_adds_shielding_difference() {
    let groups = GroupStructure::new(vec![1.0e5, 1.0e4]).unwrap();
    let shielding = uniform_shielding(0);
    let avg = &shielded_averages(
        &groups,
        &Weight::lwr(),
        &one_group_tables(),
        None,
        &[10.0],
        Some(&shielding),
        None,
    )
    .unwrap()[0];
    assert!((avg.sigma[&102][0] - 1.5).abs() < 1e-12);
    assert!((avg.sigma[&2][0] - 7.5).abs() < 1e-12);
}

#[test]
fn flux_uses_shielded_total() {
    // 群の一部だけが非分離共鳴の範囲にあると、中性子束 1/(σt + σ0) の重みが変わる
    let groups = GroupStructure::new(vec![1.0e4, 10.0]).unwrap();
    let shielding = UnresolvedShielding::new(
        1.0e3,
        1.0e6,
        1,
        vec![1.0e3, 1.0e6],
        Interpolation::LinLin,
        vec![vec![[1.0, 1.0, 0.0, 0.0]; 2]],
        vec![[10.0, 10.0, 0.0, 0.0]; 2],
    )
    .unwrap();
    let tables = BTreeMap::from([(1, constant(10.0)), (2, constant(10.0))]);
    let s0 = 1.0;
    let avg = &shielded_averages(
        &groups,
        &Weight::lwr(),
        &tables,
        None,
        &[s0],
        Some(&shielding),
        None,
    )
    .unwrap()[0];
    // 1/E の重みで、[1e3, 1e4] は σt = 1（幅 ln 10）、[10, 1e3] は σt = 10（幅 2 ln 10）
    let (inside, outside) = (1.0 / (1.0 + s0), 2.0 / (10.0 + s0));
    let expected = (1.0 * inside + 10.0 * outside) / (inside + outside);
    assert!((avg.sigma[&1][0] - expected).abs() < 1e-9 * expected);
}

#[test]
fn shielding_must_match_sigma0_count() {
    let shielding = uniform_shielding(1);
    let result = shielded_averages(
        &GroupStructure::new(vec![1.0e5, 1.0e4]).unwrap(),
        &Weight::lwr(),
        &one_group_tables(),
        None,
        &[10.0, 1.0],
        Some(&shielding),
        None,
    );
    assert!(result.is_err());
}
