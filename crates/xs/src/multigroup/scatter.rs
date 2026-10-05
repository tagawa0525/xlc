//! 散乱行列（P0）。
//!
//! 弾性散乱と離散準位の非弾性散乱は、静止した標的との重心系等方の2体反応として扱う。
//! 熱領域の弾性散乱は、指定があれば自由ガス模型で扱う（標的の熱運動による上方散乱を含む）。
//! 連続的な中性子放出反応（MT91, (n,2n) など）は MF6 のエネルギー分布と多重度を使う。

use std::collections::BTreeMap;

use super::outgoing::uniform_fractions;
use super::{
    FreeGas, GroupStructure, GroupedOutgoing, Weight, free_gas_group_integrals, group_integrals,
    group_integrals_multi,
};
use crate::endf::{Frame, Mf6, Tab1};

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

/// 散乱行列 σ(g→g') (barn)。[入射の群][出射の群]。
///
/// - `tables`: ポイントワイズの断面積（キーは MT）
/// - `q_values`: 離散準位の非弾性散乱（MT51〜90）の Q 値 (eV)
/// - `mf6`: MF6 の中性子放出データ（キーは MT）
/// - `thermal`: 熱領域の弾性散乱を自由ガス模型で扱う条件（`None` なら全域で静止標的）
///
/// 中性子を放出するのに分布のない反応があればエラーにする（中性子を黙って失わないため）。
pub fn scattering_matrix(
    groups: &GroupStructure,
    weight: &Weight,
    awr: f64,
    tables: &BTreeMap<u32, Tab1>,
    q_values: &BTreeMap<u32, f64>,
    mf6: &BTreeMap<u32, Mf6>,
    thermal: Option<&FreeGas>,
) -> Result<Vec<Vec<f64>>, String> {
    if let Some(free_gas) = thermal
        && !(free_gas.kt > 0.0 && free_gas.kt.is_finite() && free_gas.cutoff.is_finite())
    {
        return Err(format!(
            "free gas needs positive kT and a finite cutoff (kT {}, cutoff {})",
            free_gas.kt, free_gas.cutoff
        ));
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
        let rows = if let Some(data) = mf6.get(&mt) {
            emission_rows(groups, weight, table, data)
        } else if mt == 2 {
            elastic_rows(groups, weight, awr, table, thermal)
        } else if (51..=90).contains(&mt) {
            let q = *q_values
                .get(&mt)
                .ok_or_else(|| format!("MT{mt}: missing Q value"))?;
            two_body_rows(groups, weight, awr, q, table)
        } else {
            return Err(format!("MT{mt}: no outgoing neutron distribution"));
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
fn emits_neutrons(mt: u32) -> bool {
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

/// 弾性散乱の ∫_g σ(E) w(E) P(E→g') dE。`thermal` の cutoff 未満では自由ガス模型、
/// それ以上では静止標的との2体反応とする。
fn elastic_rows(
    groups: &GroupStructure,
    weight: &Weight,
    awr: f64,
    table: &Tab1,
    thermal: Option<&FreeGas>,
) -> Vec<Vec<f64>> {
    let Some(free_gas) = thermal else {
        return two_body_rows(groups, weight, awr, 0.0, table);
    };
    let mut breakpoints = table.x.clone();
    breakpoints.extend(two_body_kinks(groups, awr, 0.0));
    breakpoints.push(free_gas.cutoff);
    group_integrals_multi(groups, weight, &breakpoints, groups.len(), |e, out| {
        let sigma = table.eval(e).unwrap_or(0.0);
        if sigma == 0.0 {
            return;
        }
        if e < free_gas.cutoff {
            let transfer = free_gas_group_integrals(groups, e, awr, free_gas.kt);
            for (o, t) in out.iter_mut().zip(transfer) {
                *o += sigma * t;
            }
        } else if let Some((lo, hi)) = two_body_range(e, awr, 0.0) {
            uniform_fractions(groups, lo, hi, sigma, out);
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
        assert!(err.contains("MT4"), "{err}");
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
        // 同じ温度の自由ガスなら群間でも φ_g σ(g→g') = φ_g' σ(g'→g) が成り立つ
        let groups = GroupStructure::wims69();
        let weight = Weight::lwr();
        let tables = BTreeMap::from([(2, constant(4.0))]);
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
            Some(&thermal),
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

    #[test]
    fn elastic_uses_free_gas_even_with_mf6_data() {
        // MF6 に MT2 があっても、弾性散乱は2体反応と自由ガスで扱う
        let groups = GroupStructure::new(vec![1.0, 0.01]);
        let tables = BTreeMap::from([(2, constant(4.0))]);
        let mf6 = BTreeMap::from([(
            2,
            Mf6 {
                awr: 1.0,
                neutrons: vec![],
            },
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
            Some(&thermal),
        )
        .unwrap();
        assert!(m[0].iter().sum::<f64>() > 4.0, "{m:?}");
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
                Some(&thermal),
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
        assert!(err.contains("MT16"), "{err}");
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
            Mf6 {
                awr: 15.86,
                neutrons: vec![emission],
            },
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
}
