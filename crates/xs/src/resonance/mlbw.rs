//! 多準位 Breit-Wigner（MLBW）形式の断面積（0 K）。NJOY2016 reconr の csmlbw と同じ手順。

use std::f64::consts::PI;

use super::Sigma;
use super::factors::{channel_radius, penetrability, phase_shift, shift};
use crate::endf::{Mlbw, OrbitalL, Tab1, scattering_radius_at};
use crate::physics::wave_number_factor;

/// 断面積の計算に使う形に前処理した MLBW パラメータ。
#[derive(Debug, Clone)]
pub struct MlbwXs {
    spi: f64,
    /// 位相シフトに使う散乱半径 (10^-12 cm)（`radius_table` がなければ）
    scattering_radius: f64,
    /// エネルギー依存の散乱半径（NRO=1）
    radius_table: Option<Tab1>,
    l_values: Vec<LValue>,
}

/// MLBW の断面積の計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum MlbwError {
    #[error("the scattering radius table does not cover {energy} eV")]
    RadiusTable { energy: f64 },
}

#[derive(Debug, Clone)]
struct LValue {
    l: OrbitalL,
    arat: f64,
    /// 透過因子・シフト因子に使うチャネル半径
    channel_radius: f64,
    resonances: Vec<Resonance>,
}

#[derive(Debug, Clone, Copy)]
struct Resonance {
    er: f64,
    aj: f64,
    gn: f64,
    gg: f64,
    gf: f64,
    /// 共鳴エネルギーでのシフト因子と透過因子
    shift: f64,
    penetrability: f64,
}

impl MlbwXs {
    /// 波数とチャネル半径は、l ごとの AWRI から求める。
    pub fn new(data: &Mlbw) -> Self {
        let cwaven = wave_number_factor();
        let l_values = data
            .l_values
            .iter()
            .map(|lv| {
                let channel_radius = if data.naps == 1 {
                    data.ap
                } else {
                    channel_radius(lv.awri)
                };
                let arat = lv.awri / (lv.awri + 1.0);
                let resonances = lv
                    .resonances
                    .iter()
                    .map(|r| {
                        let rho = cwaven * arat * r.er.abs().sqrt() * channel_radius;
                        Resonance {
                            er: r.er,
                            aj: r.aj,
                            gn: r.gn,
                            gg: r.gg,
                            gf: r.gf,
                            shift: shift(lv.l, rho),
                            penetrability: penetrability(lv.l, rho),
                        }
                    })
                    .collect();
                LValue {
                    l: lv.l,
                    arat,
                    channel_radius,
                    resonances,
                }
            })
            .collect();
        Self {
            spi: data.spi,
            scattering_radius: data.ap,
            radius_table: data.ap_energy.clone(),
            l_values,
        }
    }

    /// エネルギー e (eV) での断面積。散乱半径の表（NRO=1）が e を覆わなければエラー。
    pub fn cross_sections(&self, e: f64) -> Result<Sigma, MlbwError> {
        let radius = scattering_radius_at(self.scattering_radius, self.radius_table.as_ref(), e)
            .ok_or(MlbwError::RadiusTable { energy: e })?;
        let mut sum = Sigma::default();
        for lv in &self.l_values {
            let k = wave_number_factor() * lv.arat * e.abs().sqrt();
            let pifac = PI / (k * k);
            let rho = k * lv.channel_radius;
            let (se, pe) = (shift(lv.l, rho), penetrability(lv.l, rho));
            let phi = phase_shift(lv.l, k * radius);
            // 1 - cos 2φ は低エネルギーで桁落ちするので 2 sin²φ で求める
            let cos2p = 2.0 * phi.sin().powi(2);
            let sin2p = (2.0 * phi).sin();

            let fl = f64::from(lv.l.number());
            let ajmin = ((self.spi - fl).abs() - 0.5).abs();
            let nj = (self.spi + fl + 0.5 - ajmin + 1.0).round() as usize;
            let gj: Vec<f64> = (0..nj)
                .map(|j| (2.0 * (ajmin + j as f64) + 1.0) / (4.0 * self.spi + 2.0))
                .collect();
            // 2つのチャネルスピンから作れる J は gJ を1回しか数えないので、その分の硬球散乱を足す
            let missing = 2.0 * fl + 1.0 - gj.iter().sum::<f64>();

            let mut sig = vec![(0.0, 0.0); nj];
            let (mut capture, mut fission) = (0.0, 0.0);
            for r in &lv.resonances {
                let j = (r.aj - ajmin).round() as usize;
                let shifted = r.er + r.gn * (r.shift - se) / (2.0 * r.penetrability);
                let gn = r.gn * pe / r.penetrability;
                let width = gn + r.gg + r.gf;
                let x = 2.0 * (e - shifted) / width;
                let common = 2.0 * gn / width / (1.0 + x * x);
                sig[j].0 += common;
                sig[j].1 += common * x;
                let rate = common * gj[j] / width;
                fission += rate * r.gf;
                capture += rate * r.gg;
            }
            let elastic: f64 = gj
                .iter()
                .zip(&sig)
                .map(|(g, (s1, s2))| g * ((cos2p - s1).powi(2) + (sin2p + s2).powi(2)))
                .sum::<f64>()
                + 2.0 * missing * cos2p;
            sum.elastic += pifac * elastic;
            sum.fission += 2.0 * pifac * fission;
            sum.capture += 2.0 * pifac * capture;
        }
        sum.total = sum.elastic + sum.fission + sum.capture;
        Ok(sum)
    }
}

#[cfg(test)]
mod tests;
