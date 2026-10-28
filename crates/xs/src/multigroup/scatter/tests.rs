use super::*;
use crate::endf::{AngleTables, Interpolation, LegendreTable};

/// ルジャンドルの次数ごとの、全群の和
fn sums(moments: &[Vec<f64>]) -> Vec<f64> {
    moments.iter().map(|row| row.iter().sum()).collect()
}

#[test]
fn hydrogen_moments_follow_the_kinematics() {
    // A = 1 の弾性散乱は μ_lab = √(E'/E) で E' は [0, E] に一様: ⟨P1⟩ = 2/3, ⟨P2⟩ = 1/4, ⟨P3⟩ = 0
    let groups = GroupStructure::new(vec![1.0e7, 1.0e4, 3.0e3, 1.0e3, 1.0, 1.0e-5]).unwrap();
    let m = two_body_moments(&groups, 5.0e3, 1.0, 0.0, 1.0, 3, None);
    let total = sums(&m);
    for (got, expected) in total.iter().zip([1.0, 2.0 / 3.0, 0.25, 0.0]) {
        assert!((got - expected).abs() < 1e-12, "{total:?}");
    }
    // P0 は一様な振り分けと同じ
    let mut uniform = vec![0.0; groups.len()];
    let (lo, hi) = two_body_range(5.0e3, 1.0, 0.0).unwrap();
    uniform_fractions(&groups, lo, hi, 1.0, &mut uniform);
    for (a, b) in m[0].iter().zip(&uniform) {
        assert!((a - b).abs() < 1e-14, "{:?} vs {uniform:?}", m[0]);
    }
}

#[test]
fn heavy_target_moments_follow_the_kinematics() {
    // 重心系で等方な弾性散乱の実験室系の平均の余弦は 2/(3A)
    let groups = GroupStructure::new(vec![1.0e7, 1.0e3, 9.0e2, 8.5e2, 1.0e-5]).unwrap();
    let m = two_body_moments(&groups, 1.0e3, 10.0, 0.0, 2.0, 1, None);
    let total = sums(&m);
    assert!((total[0] - 2.0).abs() < 1e-12);
    assert!(
        (total[1] / total[0] - 2.0 / 30.0).abs() < 1e-12,
        "{total:?}"
    );
}

/// 重心系の分布 f(μ_c) = 1/2 + 3/2 a μ_c（全エネルギーで同じ）
fn linear_cm(a: f64) -> Mf4 {
    Mf4 {
        frame: Frame::CenterOfMass,
        angular: Angular::Legendre(LegendreTable {
            breaks: vec![2],
            laws: vec![Interpolation::LinLin],
            energies: vec![1.0e-5, 2.0e7],
            coefficients: vec![vec![a], vec![a]],
        }),
    }
}

#[test]
fn isotropic_angular_data_changes_nothing() {
    let groups = GroupStructure::new(vec![1.0e7, 1.0e3, 9.0e2, 8.5e2, 1.0e-5]).unwrap();
    let isotropic = Mf4 {
        frame: Frame::CenterOfMass,
        angular: Angular::Isotropic,
    };
    let plain = two_body_moments(&groups, 1.0e3, 10.0, 0.0, 2.0, 3, None);
    let with = two_body_moments(&groups, 1.0e3, 10.0, 0.0, 2.0, 3, Some(&isotropic));
    for (a, b) in plain.iter().flatten().zip(with.iter().flatten()) {
        assert!((a - b).abs() < 1e-12, "{plain:?} vs {with:?}");
    }
}

#[test]
fn threshold_with_angles_keeps_the_particles() {
    // しきい値ちょうどでは出射エネルギーの幅が 0。角度分布があっても P0 の和は scale
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
    let (awr, q) = (10.0_f64, -1.0e6_f64);
    let threshold = -q * (awr + 1.0) / awr;
    let m = two_body_moments(&groups, threshold, awr, q, 2.0, 1, Some(&linear_cm(0.2)));
    let total = sums(&m);
    assert!((total[0] - 2.0).abs() < 1e-12, "{m:?}");
}

#[test]
fn center_of_mass_distribution_shifts_the_outgoing_energies() {
    // A = 10、E = 1e3 eV の弾性散乱。E' = E (1 + A² + 2A μ_c)/(A+1)² は μ_c に線形なので、
    // 群 [E1, E2] に出る確率は ∫ f dμ_c（f = 1/2 + 3/2 a μ_c）
    let (awr, e, a) = (10.0_f64, 1.0e3_f64, 0.2);
    let groups = GroupStructure::new(vec![1.0e7, 1.0e3, 9.0e2, 8.5e2, 1.0e-5]).unwrap();
    let mf4 = linear_cm(a);
    let m = two_body_moments(&groups, e, awr, 0.0, 1.0, 1, Some(&mf4));
    let mu_c = |ep: f64| ((awr + 1.0).powi(2) * ep / e - 1.0 - awr * awr) / (2.0 * awr);
    let cdf = |mu: f64| {
        let mu = mu.clamp(-1.0, 1.0);
        0.5 * (mu + 1.0) + 0.75 * a * (mu * mu - 1.0)
    };
    let b = groups.bounds();
    for g in 0..groups.len() {
        let expected = cdf(mu_c(b[g])) - cdf(mu_c(b[g + 1]));
        assert!((m[0][g] - expected).abs() < 1e-12, "group {g}: {:?}", m[0]);
    }
    // P1 の和は ∫ f(μ_c) μ_lab(μ_c) dμ_c、μ_lab = (1 + A μ_c)/√(1 + A² + 2A μ_c)
    let nodes = xlc_numeric::quadrature::gauss_legendre(64).unwrap();
    let expected: f64 = nodes
        .iter()
        .map(|&(x, w)| {
            let f = 0.5 + 1.5 * a * x;
            w * f * (1.0 + awr * x) / (1.0 + awr * awr + 2.0 * awr * x).sqrt()
        })
        .sum();
    let total = sums(&m);
    assert!(
        (total[1] - expected).abs() < 1e-12,
        "{total:?} vs {expected}"
    );
}

use crate::endf::{ClampedTab1, IncidentLaw};

fn constant(value: f64) -> Tab1 {
    Tab1::linear(vec![1.0e-5, 2.0e7], vec![value, value])
}

#[test]
fn elastic_ranges() {
    assert_eq!(two_body_range(10.0, 1.0, 0.0), Some((0.0, 10.0)));
    let (lo, hi) = two_body_range(10.0, 238.0, 0.0).unwrap();
    let alpha = (237.0_f64 / 239.0).powi(2);
    assert!((lo - alpha * 10.0).abs() < 1e-12 && (hi - 10.0).abs() < 1e-12);
}

#[test]
fn two_body_range_holds_at_the_rounded_threshold() {
    // しきい値 −Q (A + 1)/A をそのまま計算した入射エネルギーでも、丸めで β² が負（−2.2e-16）になる組がある。
    // しきい値の断面積の中性子を失わないよう、しきい値以上では β² = 0 として出射エネルギーを返す
    let (awr, q): (f64, f64) = (23.885511445655023, -284446.41774354107);
    let threshold = -q * (awr + 1.0) / awr;
    assert!(1.0 + (awr + 1.0) * q / (awr * threshold) < 0.0);
    let expected = threshold / (awr + 1.0).powi(2);
    assert_eq!(
        two_body_range(threshold, awr, q),
        Some((expected, expected))
    );
}

#[test]
fn inelastic_threshold() {
    // しきい値 -Q (A+1)/A
    let (awr, q) = (236.0, -4.49e4);
    let threshold = -q * (awr + 1.0) / awr;
    assert_eq!(two_body_range(threshold * 0.999, awr, q), None);
    let (lo, hi) = two_body_range(threshold * 1.5, awr, q).unwrap();
    assert!(lo < hi && hi < threshold * 1.5);
}

#[test]
fn hydrogen_elastic_matrix() {
    // A=1 では出射エネルギーは [0, E] で一様。重み 1/E のもとで群内に残る割合は解析的に求まる
    let groups = GroupStructure::new(vec![1.0e4, 1.0e3, 1.0e2]).unwrap();
    let tables = BTreeMap::from([(2, constant(20.0))]);
    let m = scattering_matrix(
        &groups,
        &Weight::lwr(),
        1.0,
        &Reactions {
            tables: &tables,
            q_values: &BTreeMap::new(),
            emissions: &BTreeMap::new(),
            angular: &BTreeMap::new(),
        },
        None,
    )
    .unwrap();
    let stay = 20.0 * (1.0 - 0.9 / 10.0_f64.ln());
    assert!((m[0][0] - stay).abs() < 1e-9 * stay, "{m:?}");
    // 群構造の下に出る分も最下位の群に入るので、行の和は σ
    for row in &m {
        assert!((row.iter().sum::<f64>() - 20.0).abs() < 1e-9);
    }
}

#[test]
fn hydrogen_moments_keep_the_mean_cosine() {
    // 重心系で等方な A = 1 の弾性散乱は、入射エネルギーによらず実験室系の平均の余弦が 2/3 なので、
    // どの入射の群でも P1 の行の和は P0 の行の和の 2/3。P0 は scattering_matrix と同じ
    let groups = GroupStructure::new(vec![1.0e4, 1.0e3, 1.0e2]).unwrap();
    let tables = BTreeMap::from([(2, constant(20.0))]);
    let args = (&BTreeMap::new(), &BTreeMap::new());
    let moments = scattering_moments(
        &groups,
        &Weight::lwr(),
        1.0,
        &Reactions {
            tables: &tables,
            q_values: args.0,
            emissions: args.1,
            angular: &BTreeMap::new(),
        },
        None,
        3,
    )
    .unwrap();
    assert_eq!(moments.len(), 4);
    let p0 = scattering_matrix(
        &groups,
        &Weight::lwr(),
        1.0,
        &Reactions {
            tables: &tables,
            q_values: args.0,
            emissions: args.1,
            angular: &BTreeMap::new(),
        },
        None,
    )
    .unwrap();
    assert_eq!(moments[0], p0);
    for (row0, row1) in moments[0].iter().zip(&moments[1]) {
        let ratio = row1.iter().sum::<f64>() / row0.iter().sum::<f64>();
        assert!((ratio - 2.0 / 3.0).abs() < 1e-10, "{ratio}");
    }
}

#[test]
fn discrete_inelastic_uses_q_value() {
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
    let (awr, q) = (236.0, -4.49e4);
    let tables = BTreeMap::from([(2, constant(10.0)), (51, constant(1.0))]);
    let q_values = BTreeMap::from([(51, q)]);
    let m = scattering_matrix(
        &groups,
        &Weight::lwr(),
        awr,
        &Reactions {
            tables: &tables,
            q_values: &q_values,
            emissions: &BTreeMap::new(),
            angular: &BTreeMap::new(),
        },
        None,
    )
    .unwrap();
    // 離散準位の分は高い群から低い群へも移る。行の和は σel + σ51（しきい値より上）
    assert!(m[0][1] > 0.0);
    assert!((m[0].iter().sum::<f64>() - 11.0).abs() < 1e-9);
}

#[test]
fn two_body_mf6_level_uses_q_value() {
    // MF6 の LAW=2（2体反応）で与えられた準位は、MF6 がないときと同じく Q 値の2体反応
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
    let (awr, q) = (118.87, -1.17e6);
    let tables = BTreeMap::from([(2, constant(10.0)), (51, constant(1.0))]);
    let q_values = BTreeMap::from([(51, q)]);
    let two_body = Emission::Mf6(Mf6 {
        awr,
        neutrons: vec![],
        two_body: true,
    });
    let matrix = |emissions: &BTreeMap<u32, Emission>| {
        scattering_matrix(
            &groups,
            &Weight::lwr(),
            awr,
            &Reactions {
                tables: &tables,
                q_values: &q_values,
                emissions,
                angular: &BTreeMap::new(),
            },
            None,
        )
    };
    let with_mf6 = matrix(&BTreeMap::from([(51, two_body.clone())])).unwrap();
    assert_eq!(with_mf6, matrix(&BTreeMap::new()).unwrap());
    // 準位以外の反応が LAW=2 なら、Q 値の2体反応とは扱えない
    let tables = BTreeMap::from([(2, constant(10.0)), (16, constant(1.0))]);
    let err = scattering_matrix(
        &groups,
        &Weight::lwr(),
        awr,
        &Reactions {
            tables: &tables,
            q_values: &q_values,
            emissions: &BTreeMap::from([(16, two_body)]),
            angular: &BTreeMap::new(),
        },
        None,
    )
    .unwrap_err();
    assert!(matches!(err, ScatterError::TwoBody(16)), "{err}");
}

#[test]
fn total_inelastic_without_partials_needs_a_distribution() {
    // MT51〜91 がなければ、MT4 自身を数える（分布がなければエラー）
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6]).unwrap();
    let tables = BTreeMap::from([(4, constant(1.0))]);
    let err = scattering_matrix(
        &groups,
        &Weight::lwr(),
        15.86,
        &Reactions {
            tables: &tables,
            q_values: &BTreeMap::new(),
            emissions: &BTreeMap::new(),
            angular: &BTreeMap::new(),
        },
        None,
    )
    .unwrap_err();
    assert!(matches!(err, ScatterError::NoDistribution(4)), "{err}");
}

#[test]
fn reactions_without_neutrons_need_no_distribution() {
    // MT155 (z,tα) などは MT152〜200 の範囲でも中性子を出さない
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6]).unwrap();
    for mt in [155, 182, 191, 192, 193, 197] {
        let tables = BTreeMap::from([(mt, constant(1.0))]);
        let m = scattering_matrix(
            &groups,
            &Weight::lwr(),
            15.86,
            &Reactions {
                tables: &tables,
                q_values: &BTreeMap::new(),
                emissions: &BTreeMap::new(),
                angular: &BTreeMap::new(),
            },
            None,
        )
        .unwrap();
        assert_eq!(m[0][0], 0.0, "MT{mt}");
    }
}

#[test]
fn free_gas_matrix_satisfies_detailed_balance() {
    // 重み関数の熱領域（0.1 eV 以下）は kT = 0.0253 eV の Maxwell 分布なので、
    // 同じ温度の自由ガスなら群間でも φ_g σ(g→g') = φ_g' σ(g'→g) が成り立つ。
    // 自由原子の断面積 4 b を温度に拡げた 4 N(E) を与える（表の補間誤差は 1e-7 以下）
    // 熱領域だけを使うため WIMS 69群の 0.2 eV 以下に絞る。群構造の外へ散乱した分は端の群に入る（free_gas.rs）ので、最上位の群はバランスの対象にならず、対象の群（0.1 eV 以下）より上に1群以上を残す
    let groups = GroupStructure::new(
        GroupStructure::wims69()
            .bounds()
            .iter()
            .copied()
            .filter(|&e| e <= 0.2)
            .collect(),
    )
    .unwrap();
    let weight = Weight::lwr();
    let tables = BTreeMap::from([(2, broadened_constant(4.0, 15.86, 0.0253))]);
    let thermal = FreeGas {
        kt: 0.0253,
        cutoff: 4.0,
    };
    let m = scattering_matrix(
        &groups,
        &weight,
        15.86,
        &Reactions {
            tables: &tables,
            q_values: &BTreeMap::new(),
            emissions: &BTreeMap::new(),
            angular: &BTreeMap::new(),
        },
        Some(&Thermal::FreeGas(thermal)),
    )
    .unwrap();
    let flux = group_integrals(&groups, &weight, &[], |_| 1.0);
    // 最下位の群は除く。1e-5 eV より下へ散乱した分はこの群に入れる（中性子を失わない）が、
    // 1e-5 eV より下から入射する中性子は積分しないので、この群だけはバランスが厳密でない
    let maxwellian: Vec<usize> = (0..groups.len() - 1)
        .filter(|&g| groups.range(g).1 <= 0.1)
        .collect();
    assert!(maxwellian.len() > 10);
    for &g in &maxwellian {
        // 標的の熱運動で、実効的な散乱断面積は自由原子の値より大きい
        let row: f64 = m[g].iter().sum();
        assert!(row > 4.0 && row < 6.0, "{row}");
        for &h in &maxwellian {
            let (forward, backward) = (flux[g] * m[g][h], flux[h] * m[h][g]);
            assert!(
                (forward - backward).abs() <= 1e-6 * forward.max(backward) + 1e-12,
                "{g}→{h}: {forward} vs {backward}"
            );
        }
    }
    // 熱群から上方散乱がある
    let g = groups.find(0.0253).unwrap();
    assert!(m[g][g - 1] > 0.0);
}

fn sab_matrix(elastic: f64) -> (GroupStructure, Vec<Vec<f64>>, f64) {
    use crate::multigroup::sab_tests::coarse_free_gas_table;
    let groups = GroupStructure::wims69();
    // 重み関数の Maxwell 分布（kT = 0.0253 eV）と同じ温度にする（自由ガスの S は α, β で書けば温度によらない）
    let temperature = 0.0253 / crate::physics::BK;
    let mut data = coarse_free_gas_table();
    data.temperatures = vec![temperature];
    let cutoff = 1.0;
    let sab = ThermalSab::new(
        &data,
        temperature,
        cutoff,
        crate::multigroup::SabOutside::Zero,
    )
    .unwrap();
    let tables = BTreeMap::from([(2, constant(elastic))]);
    let m = scattering_matrix(
        &groups,
        &Weight::lwr(),
        1.0,
        &Reactions {
            tables: &tables,
            q_values: &BTreeMap::new(),
            emissions: &BTreeMap::new(),
            angular: &BTreeMap::new(),
        },
        Some(&Thermal::Sab {
            sab,
            free_gas: FreeGas { kt: 0.0253, cutoff },
        }),
    )
    .unwrap();
    (groups, m, cutoff)
}

/// A = 1 の S(α, β) を 1 eV 未満で使い、その上 `free_gas_upper` (eV) 未満を自由ガスで扱う、一定の弾性散乱断面積
/// 20 b の散乱行列
fn sab_and_free_gas(free_gas_upper: f64) -> Result<Vec<Vec<f64>>, ScatterError> {
    use crate::multigroup::sab_tests::coarse_free_gas_table;
    let kt = 0.0253;
    let temperature = kt / crate::physics::BK;
    let mut data = coarse_free_gas_table();
    data.temperatures = vec![temperature];
    data.effective_temperature = Tab1::linear(vec![temperature], vec![temperature]);
    let sab = ThermalSab::new(
        &data,
        temperature,
        1.0,
        crate::multigroup::SabOutside::ShortCollisionTime,
    )
    .unwrap();
    let tables = BTreeMap::from([(2, constant(20.0))]);
    scattering_matrix(
        &GroupStructure::wims69(),
        &Weight::lwr(),
        1.0,
        &Reactions {
            tables: &tables,
            q_values: &BTreeMap::new(),
            emissions: &BTreeMap::new(),
            angular: &BTreeMap::new(),
        },
        Some(&Thermal::Sab {
            sab,
            free_gas: FreeGas {
                kt,
                cutoff: free_gas_upper,
            },
        }),
    )
}

#[test]
fn free_gas_continues_above_the_sab_cutoff() {
    // S(α, β) の上限から自由ガスの上限までは、標的の熱運動で上方散乱があり、行の和は MF3 の弾性散乱断面積に
    // なる。その上は静止した標的なので上方散乱がない
    let groups = GroupStructure::wims69();
    let m = sab_and_free_gas(4.0).unwrap();
    let weight = Weight::lwr();
    let table = constant(20.0);
    let flux = group_integrals(&groups, &weight, &table.x, |_| 1.0);
    let sigma = group_integrals(&groups, &weight, &table.x, |e| table.eval(e).unwrap());
    let (mut free_gas, mut at_rest) = (0, 0);
    for g in 1..groups.len() {
        let (low, high) = groups.range(g);
        let up: f64 = m[g][..g].iter().sum();
        if low >= 1.0 && high <= 4.0 {
            free_gas += 1;
            assert!(up > 1e-4 * sigma[g] / flux[g], "group {g}: {up}");
            let row: f64 = m[g].iter().sum();
            let expected = sigma[g] / flux[g];
            assert!(
                (row / expected - 1.0).abs() < 1e-9,
                "group {g}: {row} vs {expected}"
            );
        } else if low >= 4.0 {
            at_rest += 1;
            assert_eq!(up, 0.0, "group {g}");
        }
    }
    assert!(free_gas >= 3 && at_rest >= 3);
    // 自由ガスの上限が S(α, β) の上限より低いのはエラー
    assert!(matches!(
        sab_and_free_gas(0.5),
        Err(ScatterError::FreeGasBelowSab { .. })
    ));
}

#[test]
fn sab_matrix_satisfies_detailed_balance() {
    // 重み関数の熱領域（0.1 eV 以下）は kT = 0.0253 eV の Maxwell 分布で、S(α, β) も同じ温度なので、
    // 群間でも φ_g σ(g→g') = φ_g' σ(g'→g) が成り立つ（表の細かさによらない）。
    // 自由ガスとの物理的な一致は、核と群への積分の段階で確かめている（sab のテスト）
    let (groups, m, _) = sab_matrix(1.0);
    let flux = group_integrals(&groups, &Weight::lwr(), &[], |_| 1.0);
    let maxwellian: Vec<usize> = (0..groups.len() - 1)
        .filter(|&g| groups.range(g).1 <= 0.1)
        .collect();
    assert!(maxwellian.len() > 10);
    for &g in &maxwellian {
        for &h in &maxwellian {
            let (forward, backward) = (flux[g] * m[g][h], flux[h] * m[h][g]);
            assert!(
                (forward - backward).abs() <= 1e-6 * forward.max(backward) + 1e-12,
                "{g}→{h}: {forward} vs {backward}"
            );
        }
    }
}

#[test]
fn sab_ignores_mf3_elastic_below_cutoff() {
    // S(α, β) の核は断面積そのものなので、cutoff 未満の行は MF3 の弾性散乱断面積によらない
    let (groups, one, cutoff) = sab_matrix(1.0);
    let (_, ten, _) = sab_matrix(10.0);
    for g in (0..groups.len()).filter(|&g| groups.range(g).1 <= cutoff) {
        assert_eq!(one[g], ten[g], "group {g}");
    }
}

#[test]
fn elastic_uses_free_gas_even_with_mf6_data() {
    // MF6 に MT2 があっても、弾性散乱は2体反応と自由ガスで扱う
    let groups = GroupStructure::new(vec![1.0, 0.01, 1.0e-4]).unwrap();
    let tables = BTreeMap::from([(2, constant(4.0))]);
    let mf6 = BTreeMap::from([(
        2,
        Emission::Mf6(Mf6 {
            awr: 1.0,
            neutrons: vec![],
            two_body: false,
        }),
    )]);
    let thermal = FreeGas {
        kt: 0.0253,
        cutoff: 4.0,
    };
    let m = scattering_matrix(
        &groups,
        &Weight::lwr(),
        1.0,
        &Reactions {
            tables: &tables,
            q_values: &BTreeMap::new(),
            emissions: &mf6,
            angular: &BTreeMap::new(),
        },
        Some(&Thermal::FreeGas(thermal)),
    )
    .unwrap();
    // 2体反応では上方散乱しない
    assert!(m[1][0] > 0.0, "{m:?}");
}

#[test]
fn invalid_free_gas_is_error() {
    let groups = GroupStructure::new(vec![1.0, 0.01]).unwrap();
    let tables = BTreeMap::from([(2, constant(4.0))]);
    for (kt, cutoff) in [(0.0, 4.0), (-0.01, 4.0), (0.0253, f64::NAN)] {
        let thermal = FreeGas { kt, cutoff };
        let result = scattering_matrix(
            &groups,
            &Weight::lwr(),
            15.86,
            &Reactions {
                tables: &tables,
                q_values: &BTreeMap::new(),
                emissions: &BTreeMap::new(),
                angular: &BTreeMap::new(),
            },
            Some(&Thermal::FreeGas(thermal)),
        );
        assert!(result.is_err(), "kt {kt}, cutoff {cutoff}");
    }
}

#[test]
fn neutron_emission_without_distribution_is_error() {
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6]).unwrap();
    let tables = BTreeMap::from([(16, constant(1.0))]);
    let err = scattering_matrix(
        &groups,
        &Weight::lwr(),
        15.86,
        &Reactions {
            tables: &tables,
            q_values: &BTreeMap::new(),
            emissions: &BTreeMap::new(),
            angular: &BTreeMap::new(),
        },
        None,
    )
    .unwrap_err();
    assert!(matches!(err, ScatterError::NoDistribution(16)), "{err}");
}

#[test]
fn mf6_multiplicity_scales_the_row() {
    use crate::endf::{Frame, Interpolation, NeutronEmission, OutgoingDistributions};
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
    let emission = NeutronEmission {
        multiplicity: constant(2.0),
        frame: Frame::Laboratory,
        distributions: OutgoingDistributions {
            breaks: vec![1],
            laws: vec![IncidentLaw {
                law: Interpolation::LinLin,
                unit_base: false,
            }],
            incident: vec![1.0e-5],
            spectra: vec![Tab1::linear(vec![1.0e-5, 1.0e6], vec![1.0e-6, 1.0e-6])],
        },
    };
    let mf6 = BTreeMap::from([(
        16,
        Emission::Mf6(Mf6 {
            awr: 15.86,
            neutrons: vec![emission],
            two_body: false,
        }),
    )]);
    let tables = BTreeMap::from([(16, constant(0.5))]);
    let m = scattering_matrix(
        &groups,
        &Weight::lwr(),
        15.86,
        &Reactions {
            tables: &tables,
            q_values: &BTreeMap::new(),
            emissions: &mf6,
            angular: &BTreeMap::new(),
        },
        None,
    )
    .unwrap();
    // 多重度 2 × σ 0.5 がすべて低い群へ
    assert!(
        (m[0][1] - 1.0).abs() < 1e-6 && m[0][0].abs() < 1e-12,
        "{m:?}"
    );
}

/// θ = 100 keV、U = 1 MeV の蒸発スペクトルの MF5
fn evaporation() -> Mf5 {
    use crate::endf::EnergyLaw;
    Mf5 {
        subsections: vec![(
            constant(1.0),
            EnergyLaw::Evaporation {
                restriction: 1.0e6,
                theta: ClampedTab1::new(constant(1.0e5)).unwrap(),
            },
        )],
    }
}

#[test]
fn mf5_rows_use_the_multiplicity_of_the_reaction() {
    // (n,2n) は2個、(n,nα) は1個の中性子を出す。行の和は 多重度 × σ（しきい値より上）
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
    for (mt, count) in [(16, 2.0), (22, 1.0)] {
        let tables = BTreeMap::from([(mt, constant(0.5))]);
        let emissions = BTreeMap::from([(mt, Emission::Mf5(evaporation()))]);
        let m = scattering_matrix(
            &groups,
            &Weight::lwr(),
            15.86,
            &Reactions {
                tables: &tables,
                q_values: &BTreeMap::new(),
                emissions: &emissions,
                angular: &BTreeMap::new(),
            },
            None,
        )
        .unwrap();
        let row: f64 = m[0].iter().sum();
        assert!((row - count * 0.5).abs() < 1e-9, "MT{mt}: {m:?}");
        assert_eq!(m[1], vec![0.0, 0.0], "MT{mt}");
    }
}

#[test]
fn multiplicities_follow_the_reaction_definitions() {
    for (mt, count) in [
        (11, 2.0),
        (16, 2.0),
        (17, 3.0),
        (22, 1.0),
        (24, 2.0),
        (25, 3.0),
        (28, 1.0),
        (30, 2.0),
        (37, 4.0),
        (41, 2.0),
        (42, 3.0),
        (91, 1.0),
    ] {
        assert_eq!(multiplicity(mt), Some(count), "MT{mt}");
    }
    assert_eq!(multiplicity(5), None);
}

/// 実験室系の分布 f(μ) = 1/2 + 3/2 a μ（全エネルギーで同じ）
fn linear_lab(a: f64) -> Mf4 {
    Mf4 {
        frame: Frame::Laboratory,
        ..linear_cm(a)
    }
}

#[test]
fn continuum_lab_angles_scale_the_spectrum() {
    // MF5 の出射エネルギーと MF4 の実験室系の角度は独立なので、P1 の行は a_1 × P0 の行
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
    let tables = BTreeMap::from([(16, constant(0.5))]);
    let emissions = BTreeMap::from([(16, Emission::Mf5(evaporation()))]);
    let angular = BTreeMap::from([(16, linear_lab(0.3))]);
    let m = scattering_moments(
        &groups,
        &Weight::lwr(),
        15.86,
        &Reactions {
            tables: &tables,
            q_values: &BTreeMap::new(),
            emissions: &emissions,
            angular: &angular,
        },
        None,
        1,
    )
    .unwrap();
    for (p1, p0) in m[1].iter().flatten().zip(m[0].iter().flatten()) {
        assert!((p1 - 0.3 * p0).abs() < 1e-12, "{m:?}");
    }
    assert!(m[0][0].iter().sum::<f64>() > 0.0);
}

#[test]
fn elastic_angles_move_the_outgoing_energies() {
    // 重心系で前方に偏る弾性散乱は、等方より高い出射エネルギーの群に多く出る（行の和は同じ）
    let groups = GroupStructure::new(vec![1.0e7, 1.0e3, 9.0e2, 8.5e2, 1.0e-5]).unwrap();
    let tables = BTreeMap::from([(2, constant(1.0))]);
    let solve = |angular: &BTreeMap<u32, Mf4>| {
        scattering_matrix(
            &groups,
            &Weight::lwr(),
            10.0,
            &Reactions {
                tables: &tables,
                q_values: &BTreeMap::new(),
                emissions: &BTreeMap::new(),
                angular,
            },
            None,
        )
        .unwrap()
    };
    let isotropic = solve(&BTreeMap::new());
    let forward = solve(&BTreeMap::from([(2, linear_cm(0.2))]));
    let g = 1;
    assert!(
        forward[g][g] > isotropic[g][g] * 1.01,
        "{forward:?} vs {isotropic:?}"
    );
    let sum = |m: &Vec<Vec<f64>>| m[g].iter().sum::<f64>();
    assert!((sum(&forward) - sum(&isotropic)).abs() < 1e-12 * sum(&isotropic));
}

#[test]
fn center_of_mass_continuum_angles_are_error() {
    // MF5 の放出の角度分布は実験室系で使う（重心系の分布を黙って実験室系として使わない）
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
    let tables = BTreeMap::from([(16, constant(0.5))]);
    let emissions = BTreeMap::from([(16, Emission::Mf5(evaporation()))]);
    let angular = BTreeMap::from([(16, linear_cm(0.3))]);
    let reactions = Reactions {
        tables: &tables,
        q_values: &BTreeMap::new(),
        emissions: &emissions,
        angular: &angular,
    };
    let err = scattering_moments(&groups, &Weight::lwr(), 15.86, &reactions, None, 1).unwrap_err();
    assert!(
        matches!(err, ScatterError::CenterOfMassContinuum(16)),
        "{err}"
    );
}

#[test]
fn isotropic_center_of_mass_continuum_is_accepted() {
    // O-18 の MF5 の反応は、重心系で等方（LTT=0, LCT=2）の MF4 を持つ。P1 以上は 0 として扱う
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]).unwrap();
    let tables = BTreeMap::from([(16, constant(0.5))]);
    let emissions = BTreeMap::from([(16, Emission::Mf5(evaporation()))]);
    let isotropic = Mf4 {
        frame: Frame::CenterOfMass,
        angular: Angular::Isotropic,
    };
    let angular = BTreeMap::from([(16, isotropic)]);
    let reactions = Reactions {
        tables: &tables,
        q_values: &BTreeMap::new(),
        emissions: &emissions,
        angular: &angular,
    };
    let m = scattering_moments(&groups, &Weight::lwr(), 15.86, &reactions, None, 1).unwrap();
    assert!(m[1].iter().flatten().all(|&v| v == 0.0));
    assert!(m[0][0].iter().sum::<f64>() > 0.0);
}

#[test]
fn tabulated_center_of_mass_kinks_are_integrated_exactly() {
    // 重心系の f(μ_c) が μ_c = 0 で折れる表: (−1, 0.25), (0, 0.25), (1, 1.25)。E' は μ_c に線形なので
    // 群に出る確率は累積分布 F(μ_c) の差
    let (awr, e) = (10.0_f64, 1.0e3_f64);
    let table = Tab1::linear(vec![-1.0, 0.0, 1.0], vec![0.25, 0.25, 1.25]);
    let mf4 = Mf4 {
        frame: Frame::CenterOfMass,
        angular: Angular::Tabulated(AngleTables {
            breaks: vec![2],
            laws: vec![Interpolation::LinLin],
            energies: vec![1.0e-5, 2.0e7],
            tables: vec![table.clone(), table],
        }),
    };
    let groups = GroupStructure::new(vec![1.0e7, 1.0e3, 9.0e2, 8.5e2, 1.0e-5]).unwrap();
    let m = two_body_moments(&groups, e, awr, 0.0, 1.0, 0, Some(&mf4));
    let mu_c = |ep: f64| ((awr + 1.0).powi(2) * ep / e - 1.0 - awr * awr) / (2.0 * awr);
    let cdf = |mu: f64| {
        let mu = mu.clamp(-1.0, 1.0);
        if mu < 0.0 {
            0.25 * (mu + 1.0)
        } else {
            0.25 + 0.25 * mu + 0.5 * mu * mu
        }
    };
    let b = groups.bounds();
    for g in 0..groups.len() {
        let expected = cdf(mu_c(b[g])) - cdf(mu_c(b[g + 1]));
        assert!((m[0][g] - expected).abs() < 1e-12, "group {g}: {:?}", m[0]);
    }
}

#[test]
fn lab_frame_two_body_angles_are_error() {
    let groups = GroupStructure::new(vec![1.0e7, 1.0e3, 1.0e-5]).unwrap();
    let tables = BTreeMap::from([(2, constant(1.0))]);
    let err = scattering_matrix(
        &groups,
        &Weight::lwr(),
        10.0,
        &Reactions {
            tables: &tables,
            q_values: &BTreeMap::new(),
            emissions: &BTreeMap::new(),
            angular: &BTreeMap::from([(2, linear_lab(0.2))]),
        },
        None,
    )
    .unwrap_err();
    assert!(matches!(err, ScatterError::LabTwoBody(2)), "{err}");
}

#[test]
fn mf5_without_known_multiplicity_is_error() {
    let groups = GroupStructure::new(vec![2.0e7, 1.0e6]).unwrap();
    let tables = BTreeMap::from([(5, constant(1.0))]);
    let emissions = BTreeMap::from([(5, Emission::Mf5(evaporation()))]);
    let err = scattering_matrix(
        &groups,
        &Weight::lwr(),
        15.86,
        &Reactions {
            tables: &tables,
            q_values: &BTreeMap::new(),
            emissions: &emissions,
            angular: &BTreeMap::new(),
        },
        None,
    )
    .unwrap_err();
    assert!(matches!(err, ScatterError::Multiplicity(5)), "{err}");
}

/// 一定の自由原子の断面積 `sigma` を温度に拡げた σ N(E)。N(E) = (1 + 1/(2a²)) erf(a) + e^(-a²)/(a√π)、
/// a² = A E/kT。対数で等間隔の 20000 点の表（線形補間の相対誤差は 1e-7 以下）。
fn broadened_constant(sigma: f64, awr: f64, kt: f64) -> Tab1 {
    let n = |e: f64| {
        let a = (awr * e / kt).sqrt();
        (1.0 + 0.5 / (a * a)) * xlc_numeric::error_function::erf(a)
            + (-a * a).exp() / (a * std::f64::consts::PI.sqrt())
    };
    let points = 20000;
    let (lo, hi) = (1.0e-5_f64, 2.0e7_f64);
    let x: Vec<f64> = (0..points)
        .map(|i| lo * (hi / lo).powf(i as f64 / (points - 1) as f64))
        .collect();
    let y = x.iter().map(|&e| sigma * n(e)).collect();
    Tab1::linear(x, y)
}

#[test]
fn free_gas_rows_sum_to_the_broadened_elastic_cross_section() {
    // 温度に拡げた断面積に、規格化した出射分布を掛ける（NJOY THERMR と同じ）。熱領域の共鳴も
    // 拡げた形のまま行の和に入る。散乱はすべて群構造の中に残すので、行の和は群平均の断面積
    let groups = GroupStructure::wims69();
    let weight = Weight::lwr();
    let table = Tab1::linear(
        vec![1.0e-5, 0.2, 0.3, 0.4, 2.0e7],
        vec![4.0, 4.0, 50.0, 4.0, 4.0],
    );
    let thermal = FreeGas {
        kt: 0.0253,
        cutoff: 4.0,
    };
    let m = scattering_matrix(
        &groups,
        &weight,
        235.0,
        &Reactions {
            tables: &BTreeMap::from([(2, table.clone())]),
            q_values: &BTreeMap::new(),
            emissions: &BTreeMap::new(),
            angular: &BTreeMap::new(),
        },
        Some(&Thermal::FreeGas(thermal)),
    )
    .unwrap();
    let flux = group_integrals(&groups, &weight, &table.x, |_| 1.0);
    let sigma = group_integrals(&groups, &weight, &table.x, |e| table.eval(e).unwrap());
    for g in 0..groups.len() {
        let expected = sigma[g] / flux[g];
        let row: f64 = m[g].iter().sum();
        assert!(
            (row - expected).abs() <= 1e-9 * expected,
            "group {g}: {row} vs {expected}"
        );
    }
}
