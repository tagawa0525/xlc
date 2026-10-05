//! ENDF の1核種から多群定数を作る処理の全体。

use std::collections::BTreeMap;

use crate::endf::{Mf3, Mf5, Mf6, Tape, parse_nubar};
use crate::library::Nuclide;
use crate::multigroup::{
    FreeGas, GroupStructure, Weight, emits_neutrons, fission_spectrum, reaction_averages,
    scattering_matrix,
};
use crate::physics::BK;
use crate::pointwise::reconstruct;

/// 処理の条件。
#[derive(Debug, Clone)]
pub struct Settings {
    pub groups: GroupStructure,
    pub weight: Weight,
    /// 自由ガス模型の温度 (K)
    pub temperature: f64,
    /// この入射エネルギー (eV) 未満の弾性散乱を自由ガス模型で扱う
    pub thermal_cutoff: f64,
    /// ポイントワイズの断面積の許容誤差
    pub tolerance: f64,
}

/// 熱領域で弾性散乱の実効断面積が 0 K の値と違う分だけ、全断面積を補正する。
///
/// 自由ガスの散乱行列の行の和は、標的の熱運動を含む実効的な弾性散乱断面積になる。
/// 全断面積をそれに合わせて、全断面積 = 吸収 + 散乱 の関係を保つ。
pub fn effective_total(total: &[f64], elastic: &[f64], elastic_matrix: &[Vec<f64>]) -> Vec<f64> {
    todo!()
}

/// テープ中の（ただ1つの）物質を処理する。
pub fn process_nuclide(tape: &Tape, name: &str, settings: &Settings) -> Result<Nuclide, String> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn total_follows_effective_elastic() {
        let total = vec![10.0, 30.0];
        let elastic = vec![4.0, 20.0];
        let matrix = vec![vec![3.9, 0.1], vec![2.0, 27.0]];
        assert_eq!(effective_total(&total, &elastic, &matrix), vec![10.0, 39.0]);
    }
}
