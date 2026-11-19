//! Breit-Wigner 形式（SLBW、MLBW）に共通する前処理と、エネルギーごとの因子（NJOY2016 reconr の rdf2bw と、
//! csslbw、csmlbw の共通の部分）。

use super::factors::{channel_radius, penetrability, phase_shift, shift};
use crate::endf::{BreitWigner, OrbitalL, Tab1, scattering_radius_at};
use crate::physics::wave_number_factor;

/// Breit-Wigner の断面積の計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum BreitWignerError {
    #[error("the scattering radius table does not cover {energy} eV")]
    RadiusTable { energy: f64 },
}

/// 断面積の計算に使う形に前処理した Breit-Wigner のパラメータ。
#[derive(Debug, Clone)]
pub(super) struct Levels {
    pub(super) spi: f64,
    /// 位相シフトに使う散乱半径 (10^-12 cm)（`radius_table` がなければ）
    scattering_radius: f64,
    /// エネルギー依存の散乱半径（NRO=1）
    radius_table: Option<Tab1>,
    pub(super) l_values: Vec<LValue>,
}

#[derive(Debug, Clone)]
pub(super) struct LValue {
    pub(super) l: OrbitalL,
    arat: f64,
    /// 透過因子・シフト因子に使うチャネル半径
    channel_radius: f64,
    pub(super) resonances: Vec<Resonance>,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct Resonance {
    er: f64,
    pub(super) aj: f64,
    gn: f64,
    pub(super) gg: f64,
    pub(super) gf: f64,
    /// 共鳴エネルギーでのシフト因子と透過因子
    shift: f64,
    penetrability: f64,
}

/// 1 つの l の、エネルギー E での因子。
#[derive(Debug, Clone, Copy)]
pub(super) struct Factors {
    /// π / k²
    pub(super) pifac: f64,
    /// 硬球散乱の位相シフト φ_l
    pub(super) phi: f64,
    /// E でのシフト因子と透過因子
    shift: f64,
    penetrability: f64,
}

impl Levels {
    /// 波数とチャネル半径は、l ごとの AWRI から求める。
    pub(super) fn new(data: &BreitWigner) -> Self {
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

    /// エネルギー e (eV) で位相シフトに使う散乱半径。散乱半径の表（NRO=1）が e を覆わなければエラー。
    pub(super) fn radius(&self, e: f64) -> Result<f64, BreitWignerError> {
        scattering_radius_at(self.scattering_radius, self.radius_table.as_ref(), e)
            .ok_or(BreitWignerError::RadiusTable { energy: e })
    }
}

impl LValue {
    /// エネルギー e (eV) での因子。`radius` は位相シフトに使う散乱半径。
    pub(super) fn factors(&self, e: f64, radius: f64) -> Factors {
        let k = wave_number_factor() * self.arat * e.abs().sqrt();
        let rho = k * self.channel_radius;
        Factors {
            pifac: std::f64::consts::PI / (k * k),
            phi: phase_shift(self.l, k * radius),
            shift: shift(self.l, rho),
            penetrability: penetrability(self.l, rho),
        }
    }
}

impl Resonance {
    /// E でのシフトした共鳴エネルギー E'r と中性子幅 Γn(E) = Γn P(E) / P(|Er|)。
    pub(super) fn at(&self, f: &Factors) -> (f64, f64) {
        let shifted = self.er + self.gn * (self.shift - f.shift) / (2.0 * self.penetrability);
        let gn = self.gn * f.penetrability / self.penetrability;
        (shifted, gn)
    }
}
