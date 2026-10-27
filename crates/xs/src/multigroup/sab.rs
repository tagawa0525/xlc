//! 熱中性子散乱の S(α, β)（MF7 MT4 の非干渉性非弾性散乱）による、エネルギー移行の核。
//!
//! 主散乱体1原子あたりの二重微分断面積は
//!   σ(E→E', μ) = σb/(2kT) √(E'/E) e^(-β/2) S(α, β)、σb = σfree ((A+1)/A)²、
//!   α = (E + E' - 2μ√(EE'))/(A kT)、β = (E' - E)/kT
//! （LASYM=0 の対称な S）。LAT=1 の表は α, β を kT0 = 0.0253 eV で作っているので、表を引くときは
//! kT0 を使う。μ で積分した P0 の核は α の積分になり、
//!   σ0(E→E') = σb A/(4E) e^(-β/2) (kT_表/kT) ∫ S dα_表、α_表 は (√E ∓ √E')²/(A kT_表) の間。
//! 表の外は、選べば NJOY THERMR（thermr.f90 の `sig`）と同じ条件で短時間衝突近似で補う（ADR-0133）。

use super::free_gas::{gaussian_law_moments, integrate_points};
use crate::endf::{Interpolation, Mf7Incoherent};
use crate::physics::BK;
use xlc_format::group::GroupStructure;
use xlc_format::library::SabTransfer;
use xlc_numeric::quadrature::{GAUSS_LEGENDRE_8, legendre};

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
    /// 表の外を短時間衝突近似で補うときの、物質の温度と有効温度の比 T/T_eff
    short_collision: Option<f64>,
}

/// S(α, β) の表の外の扱い（ADR-0133）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SabOutside {
    /// NJOY THERMR と同じ条件で短時間衝突近似（MF7 の有効温度）で補う: α か β が表の上限を超えるところと、
    /// ln S を補間する区間の両端の格子点に S = 0 があり、α A ≥ 30 か β ≥ 30 のところ（α, β は表の kT で測る）
    ShortCollisionTime,
    /// 0 とする。表が cutoff までの下方散乱と、表の β の上限までの上方散乱を覆うことを確かめる
    Zero,
}

/// UFG の減速計算が上方散乱を扱う幅（kT の倍数）。lattice の自由ガスの核の上端の側の区間の幅 `TAIL` と同じで、上方散乱の裾
/// e^{−(E−E')/kT} がここで e^{−28} ≈ 7e−13 に下がる（ステージのクレートは互いに依存しないので、値をそろえて書く）
pub const UFG_TAIL: f64 = 28.0;

/// THERMR が S = 0 の格子点に接する区間を短時間衝突近似で補う、α A と β の下限（thermr.f90 の `sig` の test2）
const SHORT_COLLISION_LIMIT: f64 = 30.0;

/// この値以下の S を 0 とみなす（THERMR の sabflg = −225 で、ln S ≤ −225 を 0 とする。thermr.f90 の calcem）
const LN_S_ZERO: f64 = -225.0;

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
    #[error("S(α, β) limit {cutoff} eV is above the evaluated {evaluated} eV")]
    AboveEvaluated { cutoff: f64, evaluated: f64 },
    #[error("no effective temperature of S(α, β) at {0} K")]
    EffectiveTemperature(f64),
    #[error(
        "S(α, β) limit {cutoff} eV is below {lowest} eV, so no S(α, β) scattering reaches the UFG range"
    )]
    TransferRange { cutoff: f64, lowest: f64 },
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
    /// 温度は表にあること（温度方向には補間しない。NJOY THERMR も同じ）。表の外は `outside` で扱う。
    pub fn new(
        data: &'a Mf7Incoherent,
        temperature: f64,
        cutoff: f64,
        outside: SabOutside,
    ) -> Result<Self, SabError> {
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
        if cutoff > data.e_max {
            return Err(SabError::AboveEvaluated {
                cutoff,
                evaluated: data.e_max,
            });
        }
        let short_collision = match outside {
            SabOutside::ShortCollisionTime => {
                let t_eff = data
                    .effective_temperature
                    .eval(temperature)
                    .ok_or(SabError::EffectiveTemperature(temperature))?;
                Some(temperature / t_eff)
            }
            SabOutside::Zero => {
                // 上方散乱は表の β の上限までを確かめる（その先は e^(-β/2) で無視できる）
                let e_up = cutoff + beta_max * kt_table;
                let alpha_needed = (cutoff.sqrt() + e_up.sqrt()).powi(2) / (data.awr * kt_table);
                if cutoff / kt_table > beta_max || alpha_needed > alpha_max {
                    return Err(SabError::Coverage {
                        alpha_max,
                        beta_max,
                        cutoff,
                    });
                }
                None
            }
        };
        let ln_s = data.s[t]
            .iter()
            .map(|&v| {
                let ln = v.ln();
                if ln > LN_S_ZERO {
                    ln
                } else {
                    f64::NEG_INFINITY
                }
            })
            .collect();
        let ratio = (data.awr + 1.0) / data.awr;
        Ok(Self {
            data,
            kt,
            kt_table,
            ln_s,
            sigma_bound: data.sigma_free * ratio * ratio,
            cutoff,
            short_collision,
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
    /// μ = (E + E' − α A kT_表)/(2√(EE')) は α の1次式なので、α の積分に P_l(μ(α)) を掛ける。表の P0 は区間ごとに解析的に、
    /// P1 以上は 8 点のガウス求積で積分する。短時間衝突近似で補う区間は [`gaussian_law_moments`] で積分する。
    pub fn kernel_moments(&self, e: f64, e_out: f64, order: usize) -> Vec<f64> {
        let mut table = vec![0.0; order + 1];
        let mut short = vec![0.0; order + 1];
        if e <= 0.0 || e_out <= 0.0 {
            return table;
        }
        let d = self.data;
        let b = ((e_out - e) / self.kt_table).abs();
        let beta = &d.beta;
        let na = d.alpha.len();
        let alpha_max = d.alpha[na - 1];
        let scale = d.awr * self.kt_table;
        let (root, root_out) = (e.sqrt(), e_out.sqrt());
        let lo = (root - root_out).powi(2) / scale;
        let reach = (root + root_out).powi(2) / scale;
        let mu = |alpha: f64| (e + e_out - alpha * scale) / (2.0 * root * root_out);
        // 表の α の区間 [x0, x1] を、物質の温度の α で短時間衝突近似で積分して `short` に足す
        let to_actual = self.kt_table / self.kt;
        let beta_actual = (e_out - e) / self.kt;
        let mu_actual = |alpha: f64| mu(alpha / to_actual);
        let fill = |x0: f64, x1: f64, short: &mut [f64]| {
            if let Some(ratio) = self.short_collision {
                gaussian_law_moments(
                    x0 * to_actual,
                    x1 * to_actual,
                    beta_actual,
                    ratio,
                    &mu_actual,
                    short,
                );
            }
        };
        if b > beta[beta.len() - 1] {
            fill(lo, reach, &mut short);
        } else {
            if reach > alpha_max {
                fill(lo.max(alpha_max), reach, &mut short);
            }
            let j = beta.partition_point(|&x| x <= b).clamp(1, beta.len() - 1) - 1;
            let u = (b - beta[j]) / (beta[j + 1] - beta[j]);
            // β について ln S を線形補間した、α の格子点での ln S
            let row = |i: usize| {
                let (lo, hi) = (self.ln_s[j * na + i], self.ln_s[(j + 1) * na + i]);
                if lo == f64::NEG_INFINITY || hi == f64::NEG_INFINITY {
                    f64::NEG_INFINITY
                } else {
                    lo + u * (hi - lo)
                }
            };
            let hi = reach.min(alpha_max);
            // α の区間ごとに exp(ln S) を解析的に積分する。最初の区間の下は、その区間の補間を延ばす
            let first = d.alpha.partition_point(|&x| x <= lo).clamp(1, na - 1) - 1;
            for i in first..na - 1 {
                let (a0, a1) = (d.alpha[i], d.alpha[i + 1]);
                let (x0, x1) = (if i == first { lo } else { a0 }, a1.min(hi));
                if x1 <= x0 {
                    break;
                }
                let (l0, l1) = (row(i), row(i + 1));
                if l0 == f64::NEG_INFINITY || l1 == f64::NEG_INFINITY {
                    // S = 0 の格子点に接する区間は、α A ≥ 30 か β ≥ 30 のところだけ補う
                    let start = if b >= SHORT_COLLISION_LIMIT {
                        x0
                    } else {
                        x0.max(SHORT_COLLISION_LIMIT / d.awr)
                    };
                    if start < x1 {
                        fill(start, x1, &mut short);
                    }
                } else {
                    let slope = (l1 - l0) / (a1 - a0);
                    // e^(-β/2) を指数に入れる（kT が表の kT より小さいと、物質の温度の |β| は表の β の上限を超える）
                    let start = l0 + slope * (x0 - a0) - 0.5 * beta_actual;
                    let width = x1 - x0;
                    table[0] += if (slope * width).abs() < 1.0e-12 {
                        start.exp() * width
                    } else {
                        start.exp() * (slope * width).exp_m1() / slope
                    };
                    if order > 0 {
                        let half = 0.5 * width;
                        for &(x, w) in &GAUSS_LEGENDRE_8 {
                            let t = half * (x + 1.0);
                            let s = (start + slope * t).exp();
                            for (m, p) in table.iter_mut().zip(legendre(mu(x0 + t), order)).skip(1)
                            {
                                *m += w * half * s * p;
                            }
                        }
                    }
                }
                if a1 >= hi {
                    break;
                }
            }
        }
        // どちらの分も e^(-β/2) を含む。短時間衝突近似の分は dα_表 = dα / (kT_表/kT) なので kT_表/kT を掛けない
        let common = self.sigma_bound * d.awr / (4.0 * e);
        let factor = common * self.kt_table / self.kt;
        table
            .iter()
            .zip(&short)
            .map(|(t, s)| factor * t + common * s)
            .collect()
    }

    /// UFG の減速計算の S(α, β) の移行の表（ADR-0205）。`low` は UFG の範囲の下限 (eV)。`groups` は散乱断面積（行の積分）を
    /// 求める群構造。
    pub fn ufg_transfer(&self, groups: &GroupStructure, low: f64) -> Result<SabTransfer, SabError> {
        let _ = (groups, low);
        Ok(SabTransfer {
            incident: Vec::new(),
            transfer: Vec::new(),
            density: Vec::new(),
        })
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
    /// 群境界と、核が折れ曲がる E' = E ± β_j kT（表の β の格子点）と、短時間衝突近似で補うかどうかが変わる
    /// E' = E ± 30 kT（表の kT）で区間を分け、
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
                    .chain(std::iter::once(&SHORT_COLLISION_LIMIT))
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
pub(crate) mod tests;
