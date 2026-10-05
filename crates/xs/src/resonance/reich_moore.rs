//! Reich-Moore 形式の断面積（0 K）。
//!
//! NJOY2016 reconr の csrmat と同じ手順で計算する。捕獲は消去されたチャネルとして扱い、
//! 中性子と2つの核分裂チャネルの 3×3 の R 行列を解く。レベルシフトは使わない。

use super::Sigma;
use crate::endf::ReichMoore;

/// 断面積の計算に使う形に前処理した Reich-Moore パラメータ。
#[derive(Debug, Clone)]
pub struct ReichMooreXs {}

impl ReichMooreXs {
    pub fn new(data: &ReichMoore) -> Self {
        todo!()
    }

    /// エネルギー e (eV) での断面積。
    pub fn cross_sections(&self, e: f64) -> Sigma {
        todo!()
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::endf::{RmLValue, RmResonance};
    use crate::physics::wave_number_factor;
    use std::f64::consts::PI;

    pub(crate) const AWRI: f64 = 236.006;
    pub(crate) const AP: f64 = 0.948;

    pub(crate) fn single_level(gn: f64, gg: f64, gf: f64) -> ReichMoore {
        ReichMoore {
            spi: 0.0,
            ap: AP,
            naps: 1,
            l_values: vec![RmLValue {
                awri: AWRI,
                apl: 0.0,
                l: 0,
                resonances: vec![RmResonance {
                    er: 6.674,
                    aj: 0.5,
                    gn,
                    gg,
                    gfa: gf,
                    gfb: 0.0,
                }],
            }],
        }
    }

    /// スピン 0 の標的、l=0 の単一準位での厳密解（Breit-Wigner 型）。
    pub(crate) fn single_level_exact(e: f64, gn: f64, gg: f64, gf: f64) -> Sigma {
        let er: f64 = 6.674;
        let k = wave_number_factor() * AWRI / (AWRI + 1.0) * e.sqrt();
        let kr = wave_number_factor() * AWRI / (AWRI + 1.0) * er.sqrt();
        let pifac = PI / (k * k);
        let phi = k * AP;
        // l=0 の透過因子は ρ なので、中性子幅は ρ に比例する
        let gn_e = gn * (k * AP) / (kr * AP);
        let width = gn_e + gg + gf;
        let d = er - e;
        // U = e^{-2iφ} (1 + iΓn / (d - iΓ/2))
        let denom = d * d + width * width / 4.0;
        let (qr, qi) = (-gn_e * width / 2.0 / denom, gn_e * d / denom);
        let (wr, wi) = (1.0 + qr, qi);
        let (c, s) = ((2.0 * phi).cos(), -(2.0 * phi).sin());
        let ui = c * wi + s * wr;
        // 1 - Re U を差を取らずに計算する（低エネルギーでの桁落ちを避ける）
        let one_minus_ur =
            2.0 * phi.sin().powi(2) - (2.0 * phi).cos() * qr - (2.0 * phi).sin() * qi;
        let total = 2.0 * pifac * one_minus_ur;
        let elastic = pifac * (one_minus_ur.powi(2) + ui * ui);
        let fission = pifac * gn_e * gf / denom;
        Sigma {
            total,
            elastic,
            fission,
            capture: total - elastic - fission,
        }
    }

    fn assert_close(actual: f64, expected: f64, what: &str) {
        assert!(
            (actual - expected).abs() <= 1e-9 * expected.abs().max(1e-6),
            "{what}: actual {actual}, expected {expected}"
        );
    }

    fn assert_sigma(actual: Sigma, expected: Sigma, e: f64) {
        assert_close(actual.total, expected.total, &format!("total at {e}"));
        assert_close(actual.elastic, expected.elastic, &format!("elastic at {e}"));
        assert_close(actual.fission, expected.fission, &format!("fission at {e}"));
        assert_close(actual.capture, expected.capture, &format!("capture at {e}"));
    }

    #[test]
    fn single_level_without_fission_matches_breit_wigner() {
        let xs = ReichMooreXs::new(&single_level(1.5e-3, 2.3e-2, 0.0));
        for e in [1.0e-3, 1.0, 6.0, 6.674, 6.69, 7.0, 100.0, 1.0e4] {
            assert_sigma(
                xs.cross_sections(e),
                single_level_exact(e, 1.5e-3, 2.3e-2, 0.0),
                e,
            );
        }
    }

    #[test]
    fn single_level_with_fission_matches_breit_wigner() {
        let xs = ReichMooreXs::new(&single_level(1.5e-3, 2.3e-2, 1.0e-1));
        for e in [1.0e-3, 1.0, 6.674, 6.8, 100.0] {
            assert_sigma(
                xs.cross_sections(e),
                single_level_exact(e, 1.5e-3, 2.3e-2, 1.0e-1),
                e,
            );
        }
    }

    #[test]
    fn negative_fission_amplitude_gives_same_single_level_result() {
        // 単一準位では核分裂幅の符号（振幅の符号）は結果に影響しない
        let xs = ReichMooreXs::new(&single_level(1.5e-3, 2.3e-2, -1.0e-1));
        let e = 6.7;
        assert_sigma(
            xs.cross_sections(e),
            single_level_exact(e, 1.5e-3, 2.3e-2, 1.0e-1),
            e,
        );
    }

    #[test]
    fn p_wave_potential_scattering_sums_over_all_j() {
        // スピン 0 の標的の l=1 では J=1/2, 3/2 があり、gJ の和は 2l+1 = 3。
        // 共鳴から遠いエネルギーでは、どちらの J も硬球散乱 (4π/k²) gJ sin²φ1 だけになる
        let data = ReichMoore {
            spi: 0.0,
            ap: AP,
            naps: 1,
            l_values: vec![RmLValue {
                awri: AWRI,
                apl: 0.0,
                l: 1,
                resonances: vec![RmResonance {
                    er: 1.0e4,
                    aj: 1.5,
                    gn: 1.0e-9,
                    gg: 1.0e-9,
                    gfa: 0.0,
                    gfb: 0.0,
                }],
            }],
        };
        let e: f64 = 100.0;
        let k = wave_number_factor() * AWRI / (AWRI + 1.0) * e.sqrt();
        let rho = k * AP;
        let phi = rho - rho.atan();
        let expected = 4.0 * PI / (k * k) * 3.0 * phi.sin().powi(2);
        let sigma = ReichMooreXs::new(&data).cross_sections(e);
        assert!((sigma.elastic - expected).abs() <= 1e-6 * expected);
        assert!((sigma.total - expected).abs() <= 1e-6 * expected);
    }
}
