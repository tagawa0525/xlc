//! ポイントワイズの断面積からの多群定数の作成（NJOY の GROUPR に相当）。

mod fission;
mod free_gas;
mod integrate;
mod outgoing;
mod scatter;
mod structure;
mod weight;

use std::collections::BTreeMap;

use crate::endf::Tab1;

pub use fission::fission_spectrum;
pub use free_gas::{FreeGas, free_gas_group_integrals, free_gas_kernel};
pub use integrate::{group_integrals, group_integrals_multi};
pub use outgoing::{GroupedOutgoing, outgoing_group_integrals};
pub use scatter::{emits_neutrons, scattering_matrix, two_body_range};
pub use structure::GroupStructure;
pub use weight::Weight;

/// 無限希釈の群平均断面積。配列は群の順（高エネルギーから）。
#[derive(Debug, Clone, PartialEq)]
pub struct ReactionAverages {
    /// 群ごとの重み関数の積分 ∫_g w dE
    pub flux: Vec<f64>,
    /// MT ごとの σ_g = ∫_g σ w dE / ∫_g w dE (barn)
    pub sigma: BTreeMap<u32, Vec<f64>>,
    /// νσf の群平均（核分裂がなければ `None`）
    pub nu_fission: Option<Vec<f64>>,
}

/// ポイントワイズの断面積 `tables`（lin-lin、キーは MT）を、重み関数 `weight` で群平均する。
/// νσf は ν(E)σf(E) を群平均する。σf には MT18 を使う。
pub fn reaction_averages(
    groups: &GroupStructure,
    weight: &Weight,
    tables: &BTreeMap<u32, Tab1>,
    nubar: Option<&Tab1>,
) -> ReactionAverages {
    let flux = group_integrals(groups, weight, &[], |_| 1.0);
    let average = |breakpoints: &[f64], f: &dyn Fn(f64) -> f64| -> Vec<f64> {
        group_integrals(groups, weight, breakpoints, f)
            .iter()
            .zip(&flux)
            .map(|(integral, phi)| integral / phi)
            .collect()
    };
    let sigma = tables
        .iter()
        .map(|(&mt, t)| (mt, average(&t.x, &|e| t.eval(e).unwrap_or(0.0))))
        .collect();
    let nu_fission = match (nubar, tables.get(&18)) {
        (Some(nu), Some(fission)) => {
            let breakpoints: Vec<f64> = fission.x.iter().chain(&nu.x).copied().collect();
            Some(average(&breakpoints, &|e| {
                nu.eval(e).unwrap_or(0.0) * fission.eval(e).unwrap_or(0.0)
            }))
        }
        _ => None,
    };
    ReactionAverages {
        flux,
        sigma,
        nu_fission,
    }
}

#[cfg(test)]
mod tests {
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
        let groups = GroupStructure::new(vec![2.0e6, 1.0e6, 5.0e5]);
        let tables = BTreeMap::from([(16, Tab1::linear(vec![1.5e6, 2.0e6], vec![1.0, 1.0]))]);
        let avg = reaction_averages(&groups, &Weight::lwr(), &tables, None);
        assert_eq!(avg.sigma[&16][1], 0.0);
        assert!(avg.sigma[&16][0] > 0.0 && avg.sigma[&16][0] < 1.0);
    }
}
