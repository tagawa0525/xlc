//! 核分裂スペクトル χ。

use std::collections::BTreeMap;

use super::{GroupStructure, Weight};
use crate::endf::{Mf5, Tab1};

/// 核分裂源 ν(E)σf(E)w(E) で重み付けした核分裂スペクトル χ_g（和は 1）。
/// σf には MT18 を使う。
pub fn fission_spectrum(
    groups: &GroupStructure,
    weight: &Weight,
    tables: &BTreeMap<u32, Tab1>,
    nubar: &Tab1,
    mf5: &Mf5,
) -> Vec<f64> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::endf::Tape;
    use crate::endf::mf5_tests::two_point_mf5;
    use crate::multigroup::group_integrals;

    #[test]
    fn spectrum_follows_incident_interpolation() {
        // 入射 1e-5 eV では高い群だけ、2e7 eV では低い群だけに出る。その間は E に線形
        let text = two_point_mf5(
            "1",
            ("1.0+6", "2.0+7", "5.263158-8"),
            ("1.0-5", "1.0+6", "1.0-6"),
        );
        let tape = Tape::parse(&text).unwrap();
        let mf5 = Mf5::parse(&tape.section(1, 5, 18).unwrap()).unwrap();
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e-5]);
        let weight = Weight::lwr();
        let constant = Tab1::linear(vec![1.0e-5, 2.0e7], vec![1.0, 1.0]);
        let tables = BTreeMap::from([(18, constant.clone())]);
        let chi = fission_spectrum(&groups, &weight, &tables, &constant, &mf5);

        let t = |e: f64| (e - 1.0e-5) / (2.0e7 - 1.0e-5);
        let high: f64 = group_integrals(&groups, &weight, &[], |e| 1.0 - t(e))
            .iter()
            .sum();
        let all: f64 = group_integrals(&groups, &weight, &[], |_| 1.0).iter().sum();
        assert!((chi[0] - high / all).abs() < 1e-6, "{chi:?}");
        assert!((chi.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }
}
