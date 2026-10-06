//! 熱中性子散乱の S(α, β)（MF7 MT4 の非干渉性非弾性散乱）による、エネルギー移行の核。
//!
//! 主散乱体1原子あたりの二重微分断面積は
//!   σ(E→E', μ) = σb/(2kT) √(E'/E) e^(-β/2) S(α, β)、σb = σfree ((A+1)/A)²、
//!   α = (E + E' - 2μ√(EE'))/(A kT)、β = (E' - E)/kT
//! （LASYM=0 の対称な S）。LAT=1 の表は α, β を kT0 = 0.0253 eV で作っているので、表を引くときは
//! kT0 を使う。μ で積分した P0 の核は α の積分になり、
//!   σ0(E→E') = σb A/(4E) e^(-β/2) (kT_表/kT) ∫ S dα_表、α_表 は (√E ∓ √E')²/(A kT_表) の間。

use super::GroupStructure;
use crate::endf::Mf7Incoherent;

/// 1つの温度の S(α, β) による P0 のエネルギー移行の核。
#[derive(Debug, Clone)]
pub struct ThermalSab<'a> {
    data: &'a Mf7Incoherent,
}

impl<'a> ThermalSab<'a> {
    /// 温度 `temperature` (K) の表で、入射エネルギー `cutoff` (eV) 未満に使う核を作る。
    /// 温度は表にあること。表が `cutoff` までの散乱（上方散乱は表の β の上限まで）を覆わなければエラー。
    pub fn new(data: &'a Mf7Incoherent, temperature: f64, cutoff: f64) -> Result<Self, String> {
        let _ = (temperature, cutoff);
        Ok(Self { data })
    }

    /// σ0(E→E') (barn/eV)。
    pub fn kernel(&self, e: f64, e_out: f64) -> f64 {
        let _ = (self.data, e, e_out);
        todo!()
    }

    /// 入射エネルギー e での ∫_g' σ0(E→E') dE'（群ごと）。群構造より下に散乱する分は最下位の群に入れる。
    pub fn group_integrals(&self, groups: &GroupStructure, e: f64) -> Vec<f64> {
        let _ = (groups, e);
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::endf::Interpolation;
    use crate::multigroup::{free_gas_group_integrals, free_gas_kernel};
    use crate::physics::BK;

    const TEMPERATURE: f64 = 293.6;

    /// 自由ガス（A = 1、σfree = 1 b）の対称な S(α, β) = e^(-(α²+β²)/(4α))/√(4πα) を細かい格子で表にしたもの。
    fn free_gas_table() -> Mf7Incoherent {
        let alpha: Vec<f64> = (0..=700)
            .map(|i| 1.0e-4 * 10f64.powf(i as f64 / 100.0))
            .collect();
        let beta: Vec<f64> = (0..=1200).map(|i| 0.25 * f64::from(i)).collect();
        let s = beta
            .iter()
            .flat_map(|&b| {
                alpha.iter().map(move |&a| {
                    (-(a * a + b * b) / (4.0 * a)).exp() / (4.0 * std::f64::consts::PI * a).sqrt()
                })
            })
            .collect();
        Mf7Incoherent {
            lat: 0,
            lasym: 0,
            sigma_free: 1.0,
            awr: 1.0,
            e_max: 10.0,
            alpha,
            alpha_law: Interpolation::LogLin,
            beta,
            beta_law: Interpolation::LogLin,
            temperatures: vec![TEMPERATURE],
            s: vec![s],
            effective_temperature: crate::endf::Tab1::linear(vec![TEMPERATURE], vec![TEMPERATURE]),
        }
    }

    #[test]
    fn free_gas_table_reproduces_free_gas_kernel() {
        // 表の補間（ln S を α, β について線形）の誤差の範囲で、自由ガスの解析式と一致する
        let data = free_gas_table();
        let sab = ThermalSab::new(&data, TEMPERATURE, 1.0).unwrap();
        let kt = BK * TEMPERATURE;
        let groups = GroupStructure::wims69();
        for e in [0.01, 0.0253, 0.1, 0.5] {
            let ours = sab.group_integrals(&groups, e);
            let reference = free_gas_group_integrals(&groups, e, 1.0, kt);
            let total: f64 = reference.iter().sum();
            let ours_total: f64 = ours.iter().sum();
            assert!(
                (ours_total / total - 1.0).abs() < 1e-3,
                "E={e}: {ours_total} vs {total}"
            );
            for (g, (o, r)) in ours.iter().zip(&reference).enumerate() {
                assert!(
                    (o - r).abs() <= 2e-3 * r + 1e-6 * total,
                    "E={e}, g={g}: {o} vs {r}"
                );
            }
            let e_out = 0.7 * e;
            let (k, k_ref) = (sab.kernel(e, e_out), free_gas_kernel(e, e_out, 1.0, kt));
            assert!((k / k_ref - 1.0).abs() < 2e-3, "E={e}: {k} vs {k_ref}");
        }
    }

    #[test]
    fn kernel_satisfies_detailed_balance() {
        // E σ0(E→E') e^(-E/kT) は E と E' について対称
        let data = free_gas_table();
        let sab = ThermalSab::new(&data, TEMPERATURE, 1.0).unwrap();
        let kt = BK * TEMPERATURE;
        for (e, e_out) in [(0.01, 0.05), (0.2, 0.03), (0.0253, 0.026)] {
            let forward = e * sab.kernel(e, e_out) * (-e / kt).exp();
            let backward = e_out * sab.kernel(e_out, e) * (-e_out / kt).exp();
            assert!(
                (forward / backward - 1.0).abs() < 1e-12,
                "{forward} vs {backward}"
            );
        }
    }

    #[test]
    fn rejects_missing_temperature_and_short_table() {
        let data = free_gas_table();
        assert!(ThermalSab::new(&data, 600.0, 1.0).is_err());
        // β の上限 300 は kT ≈ 0.0253 eV で 7.6 eV 分なので、10 eV までの下方散乱を覆わない
        assert!(ThermalSab::new(&data, TEMPERATURE, 10.0).is_err());
    }
}
