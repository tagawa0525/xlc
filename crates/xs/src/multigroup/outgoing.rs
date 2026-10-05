//! 出射エネルギー分布の群への振り分け。

use super::GroupStructure;
use crate::endf::{Interpolation, OutgoingDistributions, Tab1};

/// 出射エネルギーの確率密度の表を、群ごとに積分する（群構造の外に出る分は数えない）。
pub fn outgoing_group_integrals(groups: &GroupStructure, density: &Tab1) -> Vec<f64> {
    todo!()
}

/// 入射エネルギーごとの出射分布を群に振り分けたもの。入射エネルギーの間は補間する。
#[derive(Debug, Clone)]
pub struct GroupedOutgoing {
    breaks: Vec<usize>,
    laws: Vec<Interpolation>,
    incident: Vec<f64>,
    probabilities: Vec<Vec<f64>>,
}

impl GroupedOutgoing {
    pub fn new(groups: &GroupStructure, distributions: &OutgoingDistributions) -> Self {
        todo!()
    }

    /// 入射エネルギー
    pub fn incident(&self) -> &[f64] {
        &self.incident
    }

    /// 入射エネルギー e での、出射の群ごとの確率に `scale` を掛けて `out` に足す。
    /// 入射エネルギーの範囲外では、端の分布を使う。
    pub fn accumulate(&self, e: f64, scale: f64, out: &mut [f64]) {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integrates_piecewise_linear_density_by_group() {
        let groups = GroupStructure::new(vec![4.0, 2.0, 1.0]);
        let density = Tab1::linear(vec![1.0, 3.0], vec![0.5, 0.5]);
        assert_eq!(outgoing_group_integrals(&groups, &density), vec![0.5, 0.5]);
        // 群構造の外に出る分は数えない
        let density = Tab1::linear(vec![0.0, 3.0, 5.0], vec![0.2, 0.2, 0.2]);
        let r = outgoing_group_integrals(&groups, &density);
        assert!(
            (r[0] - 0.4).abs() < 1e-15 && (r[1] - 0.2).abs() < 1e-15,
            "{r:?}"
        );
    }

    #[test]
    fn linear_density_split_inside_a_segment() {
        let groups = GroupStructure::new(vec![2.0, 1.0, 0.5]);
        let density = Tab1::linear(vec![0.0, 2.0], vec![0.0, 1.0]);
        let r = outgoing_group_integrals(&groups, &density);
        // ∫ E/2 dE
        assert!(
            (r[0] - 0.75).abs() < 1e-15 && (r[1] - 0.1875).abs() < 1e-15,
            "{r:?}"
        );
    }

    #[test]
    fn interpolates_between_incident_energies() {
        let groups = GroupStructure::new(vec![4.0, 2.0, 1.0]);
        let distributions = OutgoingDistributions {
            breaks: vec![2],
            laws: vec![Interpolation::LinLin],
            incident: vec![10.0, 20.0],
            spectra: vec![
                Tab1::linear(vec![2.0, 4.0], vec![0.5, 0.5]),
                Tab1::linear(vec![1.0, 2.0], vec![1.0, 1.0]),
            ],
        };
        let grouped = GroupedOutgoing::new(&groups, &distributions);
        let mut out = vec![0.0; 2];
        grouped.accumulate(12.5, 2.0, &mut out);
        assert_eq!(out, vec![1.5, 0.5]);
        // 範囲外は端の分布
        let mut out = vec![0.0; 2];
        grouped.accumulate(30.0, 1.0, &mut out);
        assert_eq!(out, vec![0.0, 1.0]);
    }
}
