//! 散乱行列（P0）。
//!
//! 弾性散乱と離散準位の非弾性散乱は、静止した標的との重心系等方の2体反応として扱う。
//! 熱領域の弾性散乱は、指定があれば自由ガス模型で扱う（標的の熱運動による上方散乱を含む）。
//! 連続的な中性子放出反応（MT91, (n,2n) など）は MF6 のエネルギー分布と多重度を使う。
//! MF6 がなく MF4/MF5 で与えられた反応は、MF5 のエネルギー分布と、反応で決まる多重度を使う。

use std::collections::BTreeMap;

use super::outgoing::uniform_fractions;
use super::{
    FreeGas, GroupStructure, GroupedLaw, GroupedOutgoing, ThermalSab, Weight,
    free_gas_group_integrals, group_integrals, group_integrals_multi,
};
use crate::endf::{Frame, Mf5, Mf6, Tab1};

/// 重心系等方の2体反応（Q 値 q、標的の質量比 awr）で、入射エネルギー e の中性子の
/// 出射エネルギーの範囲。出射エネルギーはこの範囲で一様に分布する。しきい値未満なら `None`。
pub fn two_body_range(e: f64, awr: f64, q: f64) -> Option<(f64, f64)> {
    // E' = E/(A+1)² (1 + A²β² + 2Aβμ)、β² = 1 + (A+1)Q/(A E)、μ は重心系の散乱角の余弦
    let beta2 = 1.0 + (awr + 1.0) * q / (awr * e);
    if beta2 < 0.0 {
        return None;
    }
    let beta = beta2.sqrt();
    let c = e / (awr + 1.0).powi(2);
    Some((
        c * (1.0 - awr * beta).powi(2),
        c * (1.0 + awr * beta).powi(2),
    ))
}

/// 連続的な中性子放出反応の出射分布。
#[derive(Debug, Clone)]
pub enum Emission {
    /// MF6（多重度と出射分布）
    Mf6(Mf6),
    /// MF5（出射エネルギー分布。多重度は反応で決まる）
    Mf5(Mf5),
}

/// MF5 で与えられた反応の中性子の多重度（ENDF-6 の反応の定義による）。
/// 定義から決まらない反応（MT5 など）は `None`。
pub fn multiplicity(mt: u32) -> Option<f64> {
    match mt {
        4 | 22 | 23 | 28 | 29 | 32..=36 | 44 | 45 | 51..=91 => Some(1.0),
        11 | 16 | 24 | 30 | 41 => Some(2.0),
        17 | 25 | 42 => Some(3.0),
        37 => Some(4.0),
        _ => None,
    }
}

/// 熱領域の散乱の扱い。どちらも cutoff 未満の入射エネルギーで、弾性散乱の静止標的の2体反応を置き換える。
#[derive(Debug, Clone)]
pub enum Thermal<'a> {
    /// 自由ガス模型（温度に拡げた MT2 の断面積に、標的の熱運動による出射分布を掛ける）
    FreeGas(FreeGas),
    /// S(α, β)（MF3 の弾性散乱断面積は使わず、核そのものが断面積を与える）
    Sab(ThermalSab<'a>),
}

/// 散乱行列を作るときのエラー。
#[derive(Debug, thiserror::Error)]
pub enum ScatterError {
    #[error("free gas needs positive kT and a finite cutoff (kT {kt}, cutoff {cutoff})")]
    FreeGas { kt: f64, cutoff: f64 },
    #[error("MT{0}: unknown neutron multiplicity")]
    Multiplicity(u32),
    #[error("MT{0}: missing Q value")]
    QValue(u32),
    #[error("MT{0}: no outgoing neutron distribution")]
    NoDistribution(u32),
}

/// 散乱行列 σ(g→g') (barn)。[入射の群][出射の群]。
///
/// - `tables`: ポイントワイズの断面積（キーは MT）
/// - `q_values`: 離散準位の非弾性散乱（MT51〜90）の Q 値 (eV)
/// - `emissions`: 連続的な中性子放出反応の出射分布（キーは MT）
/// - `thermal`: 熱領域の弾性散乱の扱い（`None` なら全域で静止標的）
///
/// 中性子を放出するのに分布のない反応があればエラーにする（中性子を黙って失わないため）。
#[tracing::instrument(level = "debug", skip_all)]
pub fn scattering_matrix(
    groups: &GroupStructure,
    weight: &Weight,
    awr: f64,
    tables: &BTreeMap<u32, Tab1>,
    q_values: &BTreeMap<u32, f64>,
    emissions: &BTreeMap<u32, Emission>,
    thermal: Option<&Thermal>,
) -> Result<Vec<Vec<f64>>, ScatterError> {
    if let Some(Thermal::FreeGas(free_gas)) = thermal
        && !(free_gas.kt > 0.0 && free_gas.kt.is_finite() && free_gas.cutoff.is_finite())
    {
        return Err(ScatterError::FreeGas {
            kt: free_gas.kt,
            cutoff: free_gas.cutoff,
        });
    }
    let ng = groups.len();
    let flux = group_integrals(groups, weight, &[], |_| 1.0);
    let has_n2n_partials = tables.keys().any(|mt| (875..=891).contains(mt));
    let has_inelastic_partials = tables.keys().any(|mt| (51..=91).contains(mt));
    let mut matrix = vec![vec![0.0; ng]; ng];
    for (&mt, table) in tables {
        // 部分反応のある MT4 と MT16 は和なので、部分反応の側で数える
        if !emits_neutrons(mt)
            || (mt == 4 && has_inelastic_partials)
            || (mt == 16 && has_n2n_partials)
        {
            continue;
        }
        // 弾性散乱は MF6 の有無にかかわらず2体反応（熱領域では自由ガス）で扱う
        let rows = if mt == 2 {
            elastic_rows(groups, weight, awr, table, thermal)
        } else if let Some(Emission::Mf6(data)) = emissions.get(&mt) {
            emission_rows(groups, weight, table, data)
        } else if let Some(Emission::Mf5(data)) = emissions.get(&mt) {
            let count = multiplicity(mt).ok_or(ScatterError::Multiplicity(mt))?;
            spectrum_rows(groups, weight, table, data, count)
        } else if (51..=90).contains(&mt) {
            let q = *q_values.get(&mt).ok_or(ScatterError::QValue(mt))?;
            two_body_rows(groups, weight, awr, q, table)
        } else {
            return Err(ScatterError::NoDistribution(mt));
        };
        for (sum, row) in matrix.iter_mut().zip(rows) {
            for (s, v) in sum.iter_mut().zip(row) {
                *s += v;
            }
        }
    }
    for (row, phi) in matrix.iter_mut().zip(&flux) {
        for v in row.iter_mut() {
            *v /= phi;
        }
    }
    Ok(matrix)
}

/// 中性子を放出する反応（核分裂を除く）
pub fn emits_neutrons(mt: u32) -> bool {
    matches!(
        mt,
        2 | 4
            | 5
            | 11
            | 16
            | 17
            | 22..=25
            | 28..=30
            | 32..=37
            | 41
            | 42
            | 44
            | 45
            | 51..=91
            // ENDF-6 の MT152〜200 のうち、中性子を出さない 155, 182, 191〜193, 197 を除く
            | 152..=154
            | 156..=181
            | 183..=190
            | 194..=196
            | 198..=200
            | 875..=891
    )
}

/// 弾性散乱の ∫_g σ(E) w(E) P(E→g') dE。`thermal` の cutoff 未満では熱領域の扱い、
/// それ以上では静止標的との2体反応とする。
fn elastic_rows(
    groups: &GroupStructure,
    weight: &Weight,
    awr: f64,
    table: &Tab1,
    thermal: Option<&Thermal>,
) -> Vec<Vec<f64>> {
    let Some(thermal) = thermal else {
        return two_body_rows(groups, weight, awr, 0.0, table);
    };
    let cutoff = match thermal {
        Thermal::FreeGas(free_gas) => free_gas.cutoff,
        Thermal::Sab(sab) => sab.cutoff(),
    };
    let mut breakpoints = table.x.clone();
    breakpoints.extend(two_body_kinks(groups, awr, 0.0));
    breakpoints.push(cutoff);
    group_integrals_multi(groups, weight, &breakpoints, groups.len(), |e, out| {
        let sigma = table.eval(e).unwrap_or(0.0);
        match thermal {
            // S(α, β) の核は断面積そのものなので、MF3 の弾性散乱断面積は使わない
            Thermal::Sab(sab) if e < cutoff => {
                for (o, t) in out.iter_mut().zip(sab.group_integrals(groups, e)) {
                    *o += t;
                }
            }
            // 温度に拡げた断面積はすでに熱運動による実効断面積の増加 N(E) を含むので、
            // 出射分布は和 N(E) で規格化する（NJOY THERMR と同じ）
            Thermal::FreeGas(free_gas) if e < cutoff => {
                let transfer = free_gas_group_integrals(groups, e, awr, free_gas.kt);
                let total: f64 = transfer.iter().sum();
                for (o, t) in out.iter_mut().zip(transfer) {
                    *o += sigma * t / total;
                }
            }
            _ => {
                if let Some((lo, hi)) = two_body_range(e, awr, 0.0) {
                    uniform_fractions(groups, lo, hi, sigma, out);
                }
            }
        }
    })
}

/// 2体反応の ∫_g σ(E) w(E) P(E→g') dE。
fn two_body_rows(
    groups: &GroupStructure,
    weight: &Weight,
    awr: f64,
    q: f64,
    table: &Tab1,
) -> Vec<Vec<f64>> {
    let mut breakpoints = table.x.clone();
    breakpoints.extend(two_body_kinks(groups, awr, q));
    group_integrals_multi(groups, weight, &breakpoints, groups.len(), |e, out| {
        let sigma = table.eval(e).unwrap_or(0.0);
        if sigma == 0.0 {
            return;
        }
        if let Some((lo, hi)) = two_body_range(e, awr, q) {
            uniform_fractions(groups, lo, hi, sigma, out);
        }
    })
}

/// 出射エネルギーの範囲の端が群境界を横切る入射エネルギーと、しきい値。
/// 被積分関数はここで折れ曲がるので、積分区間の境界に入れる。
fn two_body_kinks(groups: &GroupStructure, awr: f64, q: f64) -> Vec<f64> {
    let bounds = groups.bounds();
    let (lowest, highest) = (bounds[groups.len()], bounds[0]);
    let threshold = if q < 0.0 { -q * (awr + 1.0) / awr } else { 0.0 };
    let start = threshold.max(lowest);
    if start >= highest {
        return vec![];
    }
    let mut kinks = vec![start];
    // 範囲の端は E について単調とは限らない（しきい値の近く）ので、標本点で符号の変化を探す
    const SAMPLES: usize = 400;
    let energy = |i: usize| start * (highest / start).powf(i as f64 / SAMPLES as f64);
    let ends = |e: f64| two_body_range(e, awr, q).unwrap_or((e, e));
    for &b in bounds {
        for end in [|r: (f64, f64)| r.0, |r: (f64, f64)| r.1] {
            let f = |e: f64| end(ends(e)) - b;
            for i in 0..SAMPLES {
                let (mut a, mut c) = (energy(i), energy(i + 1));
                if f(a).signum() == f(c).signum() {
                    continue;
                }
                for _ in 0..60 {
                    let m = 0.5 * (a + c);
                    if f(m).signum() == f(a).signum() {
                        a = m;
                    } else {
                        c = m;
                    }
                }
                kinks.push(0.5 * (a + c));
            }
        }
    }
    kinks
}

/// MF6 の中性子放出の ∫_g σ(E) y(E) w(E) P(E→g') dE。
fn emission_rows(
    groups: &GroupStructure,
    weight: &Weight,
    table: &Tab1,
    data: &Mf6,
) -> Vec<Vec<f64>> {
    let ng = groups.len();
    let mut rows = vec![vec![0.0; ng]; ng];
    for emission in &data.neutrons {
        let grouped = match emission.frame {
            Frame::Laboratory => GroupedOutgoing::new(groups, &emission.distributions),
            Frame::CenterOfMass => {
                GroupedOutgoing::center_of_mass(groups, &emission.distributions, data.awr)
            }
        };
        let mut breakpoints = table.x.clone();
        breakpoints.extend(&emission.multiplicity.x);
        breakpoints.extend(grouped.incident());
        let part = group_integrals_multi(groups, weight, &breakpoints, ng, |e, out| {
            let source =
                table.eval(e).unwrap_or(0.0) * emission.multiplicity.eval(e).unwrap_or(0.0);
            if source != 0.0 {
                grouped.accumulate(e, source, out);
            }
        });
        for (sum, row) in rows.iter_mut().zip(part) {
            for (s, v) in sum.iter_mut().zip(row) {
                *s += v;
            }
        }
    }
    rows
}

/// MF5 の出射エネルギー分布と多重度 `count` による ∫_g count σ(E) w(E) Σ_k p_k(E) P_k(E→g') dE。
/// MF5 の分布は実験室系。
fn spectrum_rows(
    groups: &GroupStructure,
    weight: &Weight,
    table: &Tab1,
    data: &Mf5,
    count: f64,
) -> Vec<Vec<f64>> {
    let parts: Vec<(&Tab1, GroupedLaw)> = data
        .subsections
        .iter()
        .map(|(p, law)| (p, GroupedLaw::new(groups, law)))
        .collect();
    let mut breakpoints = table.x.clone();
    for (p, grouped) in &parts {
        breakpoints.extend(&p.x);
        breakpoints.extend(grouped.breakpoints());
    }
    group_integrals_multi(groups, weight, &breakpoints, groups.len(), |e, out| {
        let source = count * table.eval(e).unwrap_or(0.0);
        if source == 0.0 {
            return;
        }
        for (p, grouped) in &parts {
            let pk = p.eval(e).unwrap_or(0.0);
            if pk != 0.0 {
                grouped.accumulate(e, source * pk, out);
            }
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let groups = GroupStructure::new(vec![1.0e4, 1.0e3, 1.0e2]);
        let tables = BTreeMap::from([(2, constant(20.0))]);
        let m = scattering_matrix(
            &groups,
            &Weight::lwr(),
            1.0,
            &tables,
            &BTreeMap::new(),
            &BTreeMap::new(),
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
    fn discrete_inelastic_uses_q_value() {
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]);
        let (awr, q) = (236.0, -4.49e4);
        let tables = BTreeMap::from([(2, constant(10.0)), (51, constant(1.0))]);
        let q_values = BTreeMap::from([(51, q)]);
        let m = scattering_matrix(
            &groups,
            &Weight::lwr(),
            awr,
            &tables,
            &q_values,
            &BTreeMap::new(),
            None,
        )
        .unwrap();
        // 離散準位の分は高い群から低い群へも移る。行の和は σel + σ51（しきい値より上）
        assert!(m[0][1] > 0.0);
        assert!((m[0].iter().sum::<f64>() - 11.0).abs() < 1e-9);
    }

    #[test]
    fn total_inelastic_without_partials_needs_a_distribution() {
        // MT51〜91 がなければ、MT4 自身を数える（分布がなければエラー）
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6]);
        let tables = BTreeMap::from([(4, constant(1.0))]);
        let err = scattering_matrix(
            &groups,
            &Weight::lwr(),
            15.86,
            &tables,
            &BTreeMap::new(),
            &BTreeMap::new(),
            None,
        )
        .unwrap_err();
        assert!(matches!(err, ScatterError::NoDistribution(4)), "{err}");
    }

    #[test]
    fn reactions_without_neutrons_need_no_distribution() {
        // MT155 (z,tα) などは MT152〜200 の範囲でも中性子を出さない
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6]);
        for mt in [155, 182, 191, 192, 193, 197] {
            let tables = BTreeMap::from([(mt, constant(1.0))]);
            let m = scattering_matrix(
                &groups,
                &Weight::lwr(),
                15.86,
                &tables,
                &BTreeMap::new(),
                &BTreeMap::new(),
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
        let groups = GroupStructure::wims69();
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
            &tables,
            &BTreeMap::new(),
            &BTreeMap::new(),
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
        let sab = ThermalSab::new(&data, temperature, cutoff).unwrap();
        let tables = BTreeMap::from([(2, constant(elastic))]);
        let m = scattering_matrix(
            &groups,
            &Weight::lwr(),
            1.0,
            &tables,
            &BTreeMap::new(),
            &BTreeMap::new(),
            Some(&Thermal::Sab(sab)),
        )
        .unwrap();
        (groups, m, cutoff)
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
        let groups = GroupStructure::new(vec![1.0, 0.01, 1.0e-4]);
        let tables = BTreeMap::from([(2, constant(4.0))]);
        let mf6 = BTreeMap::from([(
            2,
            Emission::Mf6(Mf6 {
                awr: 1.0,
                neutrons: vec![],
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
            &tables,
            &BTreeMap::new(),
            &mf6,
            Some(&Thermal::FreeGas(thermal)),
        )
        .unwrap();
        // 2体反応では上方散乱しない
        assert!(m[1][0] > 0.0, "{m:?}");
    }

    #[test]
    fn invalid_free_gas_is_error() {
        let groups = GroupStructure::new(vec![1.0, 0.01]);
        let tables = BTreeMap::from([(2, constant(4.0))]);
        for (kt, cutoff) in [(0.0, 4.0), (-0.01, 4.0), (0.0253, f64::NAN)] {
            let thermal = FreeGas { kt, cutoff };
            let result = scattering_matrix(
                &groups,
                &Weight::lwr(),
                15.86,
                &tables,
                &BTreeMap::new(),
                &BTreeMap::new(),
                Some(&Thermal::FreeGas(thermal)),
            );
            assert!(result.is_err(), "kt {kt}, cutoff {cutoff}");
        }
    }

    #[test]
    fn neutron_emission_without_distribution_is_error() {
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6]);
        let tables = BTreeMap::from([(16, constant(1.0))]);
        let err = scattering_matrix(
            &groups,
            &Weight::lwr(),
            15.86,
            &tables,
            &BTreeMap::new(),
            &BTreeMap::new(),
            None,
        )
        .unwrap_err();
        assert!(matches!(err, ScatterError::NoDistribution(16)), "{err}");
    }

    #[test]
    fn mf6_multiplicity_scales_the_row() {
        use crate::endf::{Frame, Interpolation, NeutronEmission, OutgoingDistributions};
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]);
        let emission = NeutronEmission {
            multiplicity: constant(2.0),
            frame: Frame::Laboratory,
            distributions: OutgoingDistributions {
                breaks: vec![1],
                laws: vec![Interpolation::LinLin],
                incident: vec![1.0e-5],
                spectra: vec![Tab1::linear(vec![1.0e-5, 1.0e6], vec![1.0e-6, 1.0e-6])],
            },
        };
        let mf6 = BTreeMap::from([(
            16,
            Emission::Mf6(Mf6 {
                awr: 15.86,
                neutrons: vec![emission],
            }),
        )]);
        let tables = BTreeMap::from([(16, constant(0.5))]);
        let m = scattering_matrix(
            &groups,
            &Weight::lwr(),
            15.86,
            &tables,
            &BTreeMap::new(),
            &mf6,
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
                    theta: constant(1.0e5),
                },
            )],
        }
    }

    #[test]
    fn mf5_rows_use_the_multiplicity_of_the_reaction() {
        // (n,2n) は2個、(n,nα) は1個の中性子を出す。行の和は 多重度 × σ（しきい値より上）
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]);
        for (mt, count) in [(16, 2.0), (22, 1.0)] {
            let tables = BTreeMap::from([(mt, constant(0.5))]);
            let emissions = BTreeMap::from([(mt, Emission::Mf5(evaporation()))]);
            let m = scattering_matrix(
                &groups,
                &Weight::lwr(),
                15.86,
                &tables,
                &BTreeMap::new(),
                &emissions,
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

    #[test]
    fn mf5_without_known_multiplicity_is_error() {
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6]);
        let tables = BTreeMap::from([(5, constant(1.0))]);
        let emissions = BTreeMap::from([(5, Emission::Mf5(evaporation()))]);
        let err = scattering_matrix(
            &groups,
            &Weight::lwr(),
            15.86,
            &tables,
            &BTreeMap::new(),
            &emissions,
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
            (1.0 + 0.5 / (a * a)) * crate::special::erf(a)
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
            &BTreeMap::from([(2, table.clone())]),
            &BTreeMap::new(),
            &BTreeMap::new(),
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
}
