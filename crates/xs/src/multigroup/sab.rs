//! 熱中性子散乱の S(α, β)（MF7 MT4 の非干渉性非弾性散乱）による、エネルギー移行の核。
//!
//! 主散乱体1原子あたりの二重微分断面積は
//!   σ(E→E', μ) = σb/(2kT) √(E'/E) e^(-β/2) S(α, β)、σb = σfree ((A+1)/A)²、
//!   α = (E + E' - 2μ√(EE'))/(A kT)、β = (E' - E)/kT
//! （LASYM=0 の対称な S）。LAT=1 の表は α, β を kT0 = 0.0253 eV で作っているので、表を引くときは
//! kT0 を使う。μ で積分した P0 の核は α の積分になり、
//!   σ0(E→E') = σb A/(4E) e^(-β/2) (kT_表/kT) ∫ S dα_表、α_表 は (√E ∓ √E')²/(A kT_表) の間。

use super::GroupStructure;
use super::free_gas::integrate_points;
use super::scatter::{GAUSS_8, legendre};
use crate::endf::{Interpolation, Mf7Incoherent};
use crate::physics::BK;

/// LAT=1 の表で α, β に使う kT0 (eV)（ENDF の定義）
const KT0: f64 = 0.0253;

/// 1つの温度の S(α, β) による P0 のエネルギー移行の核。
#[derive(Debug, Clone)]
pub struct ThermalSab<'a> {
    data: &'a Mf7Incoherent,
    /// 標的の温度の kT と、表を引くときの kT (eV)
    kt: f64,
    kt_table: f64,
    /// [β の番号 × α の数 + α の番号] の ln S（S = 0 は -∞）
    ln_s: Vec<f64>,
    /// 束縛原子の散乱断面積 σb (barn)
    sigma_bound: f64,
    /// この入射エネルギー (eV) 未満で使う
    cutoff: f64,
}

/// S(α, β) の表を散乱に使えないときのエラー。
#[derive(Debug, thiserror::Error)]
pub enum SabError {
    #[error("S(α, β) has an empty α or β grid")]
    EmptyGrid,
    #[error("asymmetric S(α, β) (LASYM≠0) is not supported")]
    Asymmetric,
    #[error("S(α, β) must interpolate ln S linearly in α and β")]
    Interpolation,
    #[error("temperature {temperature} K is not in S(α, β) {available:?}")]
    Temperature {
        temperature: f64,
        available: Vec<f64>,
    },
    #[error(
        "S(α, β) table (α ≤ {alpha_max}, β ≤ {beta_max}) does not cover scattering below {cutoff} eV"
    )]
    Coverage {
        alpha_max: f64,
        beta_max: f64,
        cutoff: f64,
    },
}

impl<'a> ThermalSab<'a> {
    /// 温度 `temperature` (K) の表で、入射エネルギー `cutoff` (eV) 未満に使う核を作る。
    /// 温度は表にあること（温度方向には補間しない。NJOY THERMR も同じ）。
    ///
    /// 表の範囲外の S は使わない。そのかわり、`cutoff` までの下方散乱と、表の β の上限までの
    /// 上方散乱を、表の α, β が覆うことを確かめる（上方散乱の先は e^(-β/2) で無視できる）。
    pub fn new(data: &'a Mf7Incoherent, temperature: f64, cutoff: f64) -> Result<Self, SabError> {
        if data.lasym != 0 {
            return Err(SabError::Asymmetric);
        }
        if data.alpha_law != Interpolation::LogLin || data.beta_law != Interpolation::LogLin {
            return Err(SabError::Interpolation);
        }
        let t = data
            .temperatures
            .iter()
            .position(|&t| (t - temperature).abs() <= 1.0e-9 * temperature)
            .ok_or_else(|| SabError::Temperature {
                temperature,
                available: data.temperatures.clone(),
            })?;
        let kt = BK * temperature;
        let kt_table = if data.lat == 1 { KT0 } else { kt };
        let (Some(&alpha_max), Some(&beta_max)) = (data.alpha.last(), data.beta.last()) else {
            return Err(SabError::EmptyGrid);
        };
        let e_up = cutoff + beta_max * kt_table;
        let alpha_needed = (cutoff.sqrt() + e_up.sqrt()).powi(2) / (data.awr * kt_table);
        if cutoff / kt_table > beta_max || alpha_needed > alpha_max {
            return Err(SabError::Coverage {
                alpha_max,
                beta_max,
                cutoff,
            });
        }
        let ln_s = data.s[t]
            .iter()
            .map(|&v| if v > 0.0 { v.ln() } else { f64::NEG_INFINITY })
            .collect();
        let ratio = (data.awr + 1.0) / data.awr;
        Ok(Self {
            data,
            kt,
            kt_table,
            ln_s,
            sigma_bound: data.sigma_free * ratio * ratio,
            cutoff,
        })
    }

    /// この入射エネルギー (eV) 未満で使う。
    pub fn cutoff(&self) -> f64 {
        self.cutoff
    }

    /// σ0(E→E') (barn/eV)。
    pub fn kernel(&self, e: f64, e_out: f64) -> f64 {
        self.kernel_moments(e, e_out, 0)[0]
    }

    /// ルジャンドルのモーメント σ_l(E→E') = ∫ σ(E→E', μ) P_l(μ) dμ (barn/eV)、l = 0..=order。
    /// μ = (E + E' − α A kT_表)/(2√(EE')) は α の1次式なので、α の積分に P_l(μ(α)) を掛ける。P0 は区間ごとに解析的に、
    /// P1 以上は 8 点のガウス求積で積分する。
    pub fn kernel_moments(&self, e: f64, e_out: f64, order: usize) -> Vec<f64> {
        let mut moments = vec![0.0; order + 1];
        if e <= 0.0 || e_out <= 0.0 {
            return moments;
        }
        let d = self.data;
        let b = ((e_out - e) / self.kt_table).abs();
        let beta = &d.beta;
        if b > beta[beta.len() - 1] {
            return moments;
        }
        let j = beta.partition_point(|&x| x <= b).clamp(1, beta.len() - 1) - 1;
        let u = (b - beta[j]) / (beta[j + 1] - beta[j]);
        let na = d.alpha.len();
        // β について ln S を線形補間した、α の格子点での ln S
        let row = |i: usize| {
            let (lo, hi) = (self.ln_s[j * na + i], self.ln_s[(j + 1) * na + i]);
            if lo == f64::NEG_INFINITY || hi == f64::NEG_INFINITY {
                f64::NEG_INFINITY
            } else {
                lo + u * (hi - lo)
            }
        };
        let scale = d.awr * self.kt_table;
        let (root, root_out) = (e.sqrt(), e_out.sqrt());
        let lo = (root - root_out).powi(2) / scale;
        let hi = ((root + root_out).powi(2) / scale).min(d.alpha[na - 1]);
        // α の区間ごとに exp(ln S) を解析的に積分する。最初の区間の下は、その区間の補間を延ばす
        let first = d.alpha.partition_point(|&x| x <= lo).clamp(1, na - 1) - 1;
        let mu = |alpha: f64| (e + e_out - alpha * scale) / (2.0 * root * root_out);
        for i in first..na - 1 {
            let (a0, a1) = (d.alpha[i], d.alpha[i + 1]);
            let (x0, x1) = (if i == first { lo } else { a0 }, a1.min(hi));
            if x1 <= x0 {
                break;
            }
            let (l0, l1) = (row(i), row(i + 1));
            if l0 == f64::NEG_INFINITY || l1 == f64::NEG_INFINITY {
                continue;
            }
            let slope = (l1 - l0) / (a1 - a0);
            let start = l0 + slope * (x0 - a0);
            let width = x1 - x0;
            moments[0] += if (slope * width).abs() < 1.0e-12 {
                start.exp() * width
            } else {
                start.exp() * (slope * width).exp_m1() / slope
            };
            if order > 0 {
                let half = 0.5 * width;
                for &(x, w) in &GAUSS_8 {
                    let t = half * (x + 1.0);
                    let s = (start + slope * t).exp();
                    for (m, p) in moments.iter_mut().zip(legendre(mu(x0 + t), order)).skip(1) {
                        *m += w * half * s * p;
                    }
                }
            }
            if a1 >= hi {
                break;
            }
        }
        let beta_actual = (e_out - e) / self.kt;
        let factor =
            self.sigma_bound * d.awr / (4.0 * e) * (-0.5 * beta_actual).exp() * self.kt_table
                / self.kt;
        for m in &mut moments {
            *m *= factor;
        }
        moments
    }

    /// 入射エネルギー e での ∫_g' σ0(E→E') dE'（群ごと）。群構造より下に散乱する分は最下位の群に、
    /// 上に散乱する分は最上位の群に入れる（中性子を失わない）。
    pub fn group_integrals(&self, groups: &GroupStructure, e: f64) -> Vec<f64> {
        self.group_moments(groups, e, 0)
            .into_iter()
            .next()
            .unwrap_or_default()
    }

    /// 入射エネルギー e での ∫_g' σ_l(E→E') dE'（[l][群]、l = 0..=order）。
    ///
    /// 群境界と、核が折れ曲がる E' = E ± β_j kT（表の β の格子点）で区間を分け、
    /// 自由ガスと同じく対数的に細かく分けた区間の Gauss-Legendre 則で積分する。
    pub fn group_moments(&self, groups: &GroupStructure, e: f64, order: usize) -> Vec<Vec<f64>> {
        let reach = e + self.data.beta[self.data.beta.len() - 1] * self.kt_table;
        let mut points: Vec<f64> = groups
            .bounds()
            .iter()
            .copied()
            .chain(
                self.data
                    .beta
                    .iter()
                    .flat_map(|&b| [e - b * self.kt_table, e + b * self.kt_table]),
            )
            .chain([0.0, e, reach])
            .filter(|&x| (0.0..=reach).contains(&x))
            .collect();
        points.sort_by(f64::total_cmp);
        points.dedup();
        let kernel = |x: f64| self.kernel_moments(e, x, order);
        integrate_points(groups, &points, order + 1, &kernel)
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::multigroup::{free_gas_group_integrals, free_gas_kernel};

    pub(crate) const TEMPERATURE: f64 = 293.6;

    /// 自由ガス（A = 1、σfree = 1 b）の対称な S(α, β) = e^(-(α²+β²)/(4α))/√(4πα) を細かい格子で表にしたもの。
    /// 小さな α では S が β について幅 √α ほどで急に変わるので、実際の表と同じく β も原点の近くを細かくとる
    pub(crate) fn free_gas_table() -> Mf7Incoherent {
        let alpha: Vec<f64> = (0..=800)
            .map(|i| 1.0e-5 * 10f64.powf(i as f64 / 100.0))
            .collect();
        let beta: Vec<f64> = std::iter::once(0.0)
            .chain((0..=750).map(|i| 1.0e-5 * 10f64.powf(i as f64 / 100.0)))
            .collect();
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

    /// `free_gas_table` を1桁あたり 20 点に粗くした表。補間の誤差は大きいが、詳細バランスのように
    /// 表の細かさによらない性質を、散乱行列の段階で安く確かめるのに使う。
    pub(crate) fn coarse_free_gas_table() -> Mf7Incoherent {
        let fine = free_gas_table();
        let (na, nb) = (fine.alpha.len(), fine.beta.len());
        let keep_alpha: Vec<usize> = (0..na).step_by(5).collect();
        let keep_beta: Vec<usize> = std::iter::once(0).chain((1..nb).step_by(5)).collect();
        let s = keep_beta
            .iter()
            .flat_map(|&j| keep_alpha.iter().map(move |&i| (j, i)))
            .map(|(j, i)| fine.s[0][j * na + i])
            .collect();
        Mf7Incoherent {
            alpha: keep_alpha.iter().map(|&i| fine.alpha[i]).collect(),
            beta: keep_beta.iter().map(|&j| fine.beta[j]).collect(),
            s: vec![s],
            ..fine
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
    fn kernel_moments_match_the_free_gas_double_differential() {
        // 自由ガス（A = 1、σfree = 1 b、σb = 4 b）の σ(E→E', μ) = σb/(2kT) √(E'/E) e^(-β/2) S(α, β) を
        // μ で細かく積分した P_l のモーメントと、表から求めたものが、表の補間の誤差の範囲で一致する
        let data = free_gas_table();
        let sab = ThermalSab::new(&data, TEMPERATURE, 1.0).unwrap();
        let kt = BK * TEMPERATURE;
        let reference = |e: f64, e_out: f64, l: usize| {
            let n = 200_000;
            let beta = (e_out - e) / kt;
            let f = |mu: f64| {
                let alpha = (e + e_out - 2.0 * mu * (e * e_out).sqrt()) / kt;
                let s = (-(alpha * alpha + beta * beta) / (4.0 * alpha)).exp()
                    / (4.0 * std::f64::consts::PI * alpha).sqrt();
                4.0 / (2.0 * kt)
                    * (e_out / e).sqrt()
                    * (-0.5 * beta).exp()
                    * s
                    * crate::multigroup::scatter::legendre(mu, l)[l]
            };
            // 中点則（μ = 1 の近くでも α > 0）
            (0..n)
                .map(|i| f(-1.0 + (i as f64 + 0.5) * 2.0 / n as f64))
                .sum::<f64>()
                * 2.0
                / n as f64
        };
        for (e, e_out) in [(0.05, 0.03), (0.05, 0.07), (0.0253, 0.02), (0.2, 0.15)] {
            let moments = sab.kernel_moments(e, e_out, 3);
            let p0 = reference(e, e_out, 0);
            assert!(
                (moments[0] / p0 - 1.0).abs() < 2e-3,
                "{e} {e_out}: {moments:?} vs {p0}"
            );
            for l in 1..=3 {
                let r = reference(e, e_out, l);
                assert!(
                    (moments[l] - r).abs() < 2e-3 * p0,
                    "{e} {e_out} P{l}: {} vs {r}",
                    moments[l]
                );
            }
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
    fn empty_grid_is_error() {
        let mut data = free_gas_table();
        data.beta.clear();
        assert!(matches!(
            ThermalSab::new(&data, TEMPERATURE, 1.0),
            Err(SabError::EmptyGrid)
        ));
    }

    #[test]
    fn rejects_missing_temperature_and_short_table() {
        let data = free_gas_table();
        assert!(ThermalSab::new(&data, 600.0, 1.0).is_err());
        // β の上限 316 は kT ≈ 0.0253 eV で 8 eV 分なので、10 eV までの下方散乱を覆わない
        assert!(ThermalSab::new(&data, TEMPERATURE, 10.0).is_err());
    }
}
