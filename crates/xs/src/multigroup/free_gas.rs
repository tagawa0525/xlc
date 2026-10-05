//! 自由ガス模型による熱中性子の弾性散乱。
//!
//! 標的核を温度 T の理想気体（Maxwell 分布で運動する自由な原子）とみなし、実験室系の
//! 散乱を等方として、出射エネルギーの分布 S(E→E') を解析的に求める。
//!
//! S は規格化しない。∫S dE' = N(E) は自由ガスの実効断面積と自由原子の断面積の比で、
//! 低エネルギーほど 1 より大きい。規格化すると N(E)/N(E') の分だけ詳細バランスが崩れ、
//! 特に水素で熱中性子のスペクトルが歪む。

use super::GroupStructure;

/// 熱領域の弾性散乱を自由ガス模型で扱う条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FreeGas {
    /// kT (eV)
    pub kt: f64,
    /// この入射エネルギー (eV) 未満で自由ガス模型を使う
    pub cutoff: f64,
}

/// 自由原子の散乱断面積を 1 としたときの、エネルギー移行の確率密度 S(E→E') (1/eV)。
/// `awr` は標的の質量比。
pub fn free_gas_kernel(e: f64, e_out: f64, awr: f64, kt: f64) -> f64 {
    todo!()
}

/// 入射エネルギー e での ∫_g' S(E→E') dE'（群ごと、和は N(E)）。
pub fn free_gas_group_integrals(groups: &GroupStructure, e: f64, awr: f64, kt: f64) -> Vec<f64> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::special::erf;

    const KT: f64 = 0.0253;

    fn assert_rel(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() <= tol * expected.abs(),
            "actual {actual}, expected {expected}"
        );
    }

    #[test]
    fn hydrogen_has_closed_form() {
        // A=1: E' < E では erf(√(E'/kT))/E、E' > E では exp((E-E')/kT) erf(√(E/kT))/E
        let e = 0.05;
        for e_out in [0.001, 0.02, 0.0499] {
            assert_rel(
                free_gas_kernel(e, e_out, 1.0, KT),
                erf((e_out / KT).sqrt()) / e,
                1e-12,
            );
        }
        for e_out in [0.0501, 0.1, 0.3] {
            let expected = ((e - e_out) / KT).exp() * erf((e / KT).sqrt()) / e;
            assert_rel(free_gas_kernel(e, e_out, 1.0, KT), expected, 1e-12);
        }
    }

    #[test]
    fn detailed_balance() {
        // E exp(-E/kT) S(E→E') = E' exp(-E'/kT) S(E'→E)
        for awr in [1.0, 15.86, 236.0] {
            for (e, e_out) in [(0.01, 0.03), (0.05, 0.02), (1.0, 0.98), (2.0, 3.5)] {
                let forward = e * (-e / KT).exp() * free_gas_kernel(e, e_out, awr, KT);
                let backward = e_out * (-e_out / KT).exp() * free_gas_kernel(e_out, e, awr, KT);
                assert_rel(forward, backward, 1e-9);
            }
        }
    }

    #[test]
    fn integrates_to_effective_cross_section() {
        // ∫ S dE' = (1 + 1/(2Aε)) erf(√(Aε)) + exp(-Aε)/√(πAε)、ε = E/kT
        for awr in [1.0, 15.86, 236.0] {
            for e in [0.005, 0.0253, 0.5, 3.0] {
                let expected = effective_ratio(e, awr);
                // 台形公式で積分する（核の幅に比べて十分細かい格子）
                let width = (e * KT / awr).sqrt() + KT;
                let (lo, hi) = ((e - 60.0 * width).max(0.0), e + 60.0 * width);
                let n = 400_000;
                let h = (hi - lo) / n as f64;
                let integral: f64 = (0..=n)
                    .map(|i| {
                        let x = lo + i as f64 * h;
                        let f = free_gas_kernel(e, x, awr, KT);
                        if i == 0 || i == n { 0.5 * f } else { f }
                    })
                    .sum::<f64>()
                    * h;
                assert_rel(integral, expected, 1e-5);
            }
        }
    }

    #[test]
    fn no_overflow_far_below_incident_energy() {
        // E/kT ≈ 160 で E' がずっと低くても有限
        let s = free_gas_kernel(4.0, 0.01, 15.86, KT);
        assert!(s.is_finite() && s >= 0.0, "{s}");
    }

    fn effective_ratio(e: f64, awr: f64) -> f64 {
        let a_eps = awr * e / KT;
        (1.0 + 1.0 / (2.0 * a_eps)) * erf(a_eps.sqrt())
            + (-a_eps).exp() / (std::f64::consts::PI * a_eps).sqrt()
    }

    #[test]
    fn group_integrals_sum_to_effective_ratio() {
        let groups = GroupStructure::wims69();
        for awr in [1.0, 15.86, 236.0] {
            for e in [1.0e-3, 0.0253, 0.9, 3.9] {
                let p = free_gas_group_integrals(&groups, e, awr, KT);
                assert_rel(p.iter().sum::<f64>(), effective_ratio(e, awr), 1e-6);
                assert!(p.iter().all(|&v| v >= 0.0));
            }
        }
        // 重い核ではほとんど同じ群に残る
        let e = 0.0253;
        let g = groups.find(e).unwrap();
        assert!(free_gas_group_integrals(&groups, e, 236.0, KT)[g] > 0.5);
    }
}
