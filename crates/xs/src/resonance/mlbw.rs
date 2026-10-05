//! 多準位 Breit-Wigner（MLBW）形式の断面積（0 K）。NJOY2016 reconr の csmlbw と同じ手順。

use std::f64::consts::PI;

use super::Sigma;
use super::factors::{channel_radius, penetrability, phase_shift, shift};
use crate::endf::Mlbw;
use crate::physics::wave_number_factor;

/// 断面積の計算に使う形に前処理した MLBW パラメータ。
#[derive(Debug, Clone)]
pub struct MlbwXs {
    spi: f64,
    l_values: Vec<LValue>,
}

#[derive(Debug, Clone)]
struct LValue {
    l: u32,
    arat: f64,
    /// 透過因子・シフト因子に使うチャネル半径
    channel_radius: f64,
    /// 位相シフトに使う散乱半径
    scattering_radius: f64,
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
                    scattering_radius: data.ap,
                    resonances,
                }
            })
            .collect();
        Self {
            spi: data.spi,
            l_values,
        }
    }

    pub fn cross_sections(&self, e: f64) -> Sigma {
        let mut sum = Sigma::default();
        for lv in &self.l_values {
            let k = wave_number_factor() * lv.arat * e.abs().sqrt();
            let pifac = PI / (k * k);
            let rho = k * lv.channel_radius;
            let (se, pe) = (shift(lv.l, rho), penetrability(lv.l, rho));
            let phi = phase_shift(lv.l, k * lv.scattering_radius);
            // 1 - cos 2φ は低エネルギーで桁落ちするので 2 sin²φ で求める
            let cos2p = 2.0 * phi.sin().powi(2);
            let sin2p = (2.0 * phi).sin();

            let fl = f64::from(lv.l);
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
        sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::endf::{BwLValue, BwResonance};
    use crate::physics::wave_number_factor;
    use crate::resonance::reich_moore_tests::{AP, AWRI, single_level_exact};
    use std::f64::consts::PI;

    fn single(l: u32, spi: f64, aj: f64, er: f64, (gn, gg, gf): (f64, f64, f64)) -> Mlbw {
        Mlbw {
            spi,
            ap: AP,
            naps: 1,
            l_values: vec![BwLValue {
                awri: AWRI,
                l,
                resonances: vec![BwResonance { er, aj, gn, gg, gf }],
            }],
        }
    }

    fn assert_close(actual: f64, expected: f64, what: &str) {
        assert!(
            (actual - expected).abs() <= 1e-9 * expected.abs().max(1e-6),
            "{what}: actual {actual}, expected {expected}"
        );
    }

    #[test]
    fn single_s_wave_level_matches_exact_breit_wigner() {
        // l=0 ではレベルシフトがなく、単一準位の MLBW は厳密な Breit-Wigner 形
        let widths = (1.5e-3, 2.3e-2, 1.0e-2);
        let xs = MlbwXs::new(&single(0, 0.0, 0.5, 6.674, widths));
        for e in [1.0e-3, 1.0, 6.674, 6.7, 100.0] {
            let (a, b) = (
                xs.cross_sections(e),
                single_level_exact(e, widths.0, widths.1, widths.2),
            );
            assert_close(a.elastic, b.elastic, &format!("elastic at {e}"));
            assert_close(a.capture, b.capture, &format!("capture at {e}"));
            assert_close(a.fission, b.fission, &format!("fission at {e}"));
            assert_close(a.total, b.total, &format!("total at {e}"));
        }
    }

    #[test]
    fn missing_spin_states_add_hard_sphere_scattering() {
        // 標的スピン 1/2・l=1 では J=1 が2つのチャネルスピンから作れるが、gJ は1回しか数えないので、
        // 足りない分 (2l+1) - ΣgJ の硬球散乱を足す。共鳴から遠ければ全体で (4π/k²)(2l+1) sin²φ1
        let xs = MlbwXs::new(&single(1, 0.5, 2.0, 1.0e4, (1.0e-9, 1.0e-9, 0.0)));
        let e: f64 = 100.0;
        let k = wave_number_factor() * AWRI / (AWRI + 1.0) * e.sqrt();
        let rho = k * AP;
        let phi = rho - rho.atan();
        let expected = 4.0 * PI / (k * k) * 3.0 * phi.sin().powi(2);
        let sigma = xs.cross_sections(e);
        assert!(
            (sigma.elastic - expected).abs() <= 1e-6 * expected,
            "{} vs {expected}",
            sigma.elastic
        );
    }

    #[test]
    fn level_shift_vanishes_at_resonance_energy() {
        // E = Er ではシフトが 0 なので、捕獲のピークは π g Γn Γγ / (k² Γ²/4)
        let (er, gn, gg) = (2.0e3, 0.3, 0.05);
        let xs = MlbwXs::new(&single(1, 0.0, 1.5, er, (gn, gg, 0.0)));
        let k = wave_number_factor() * AWRI / (AWRI + 1.0) * er.sqrt();
        let g = 2.0;
        let width = gn + gg;
        let expected = PI * g * gn * gg / (k * k * width * width / 4.0);
        assert_close(xs.cross_sections(er).capture, expected, "capture peak");
    }
}
