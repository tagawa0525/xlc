//! 散乱行列（P0）。
//!
//! 弾性散乱と離散準位の非弾性散乱は、静止した標的との重心系等方の2体反応として扱う。
//! 連続的な中性子放出反応（MT91, (n,2n) など）は MF6 のエネルギー分布と多重度を使う。

use std::collections::BTreeMap;

use super::{GroupStructure, Weight};
use crate::endf::{Mf6, Tab1};

/// 重心系等方の2体反応（Q 値 q、標的の質量比 awr）で、入射エネルギー e の中性子の
/// 出射エネルギーの範囲。出射エネルギーはこの範囲で一様に分布する。しきい値未満なら `None`。
pub fn two_body_range(e: f64, awr: f64, q: f64) -> Option<(f64, f64)> {
    todo!()
}

/// 散乱行列 σ(g→g') (barn)。[入射の群][出射の群]。
///
/// - `tables`: ポイントワイズの断面積（キーは MT）
/// - `q_values`: 離散準位の非弾性散乱（MT51〜90）の Q 値 (eV)
/// - `mf6`: MF6 の中性子放出データ（キーは MT）
///
/// 中性子を放出するのに分布のない反応があればエラーにする（中性子を黙って失わないため）。
pub fn scattering_matrix(
    groups: &GroupStructure,
    weight: &Weight,
    awr: f64,
    tables: &BTreeMap<u32, Tab1>,
    q_values: &BTreeMap<u32, f64>,
    mf6: &BTreeMap<u32, Mf6>,
) -> Result<Vec<Vec<f64>>, String> {
    todo!()
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
        )
        .unwrap();
        // 離散準位の分は高い群から低い群へも移る。行の和は σel + σ51（しきい値より上）
        assert!(m[0][1] > 0.0);
        assert!((m[0].iter().sum::<f64>() - 11.0).abs() < 1e-9);
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
        )
        .unwrap_err();
        assert!(err.contains("MT16"), "{err}");
    }

    #[test]
    fn mf6_multiplicity_scales_the_row() {
        use crate::endf::{Interpolation, NeutronEmission, OutgoingDistributions};
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]);
        let emission = NeutronEmission {
            multiplicity: constant(2.0),
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
        )
        .unwrap();
        // 多重度 2 × σ 0.5 がすべて低い群へ
        assert!(
            (m[0][1] - 1.0).abs() < 1e-6 && m[0][0].abs() < 1e-12,
            "{m:?}"
        );
    }
}
