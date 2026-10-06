//! 非分離共鳴の無限希釈の平均断面積（NJOY RECONR csunr2 と同じ単準位 Breit-Wigner の式）。
//!
//! 幅のゆらぎの期待値は、NJOY の求積表ではなく適応積分で求める（`fluctuation`）。

use super::Sigma;
use crate::endf::{Interpolation, Unresolved};

/// 非分離共鳴のパラメータから、無限希釈の平均断面積を計算する。
pub struct UnresolvedXs<'a> {
    urr: &'a Unresolved,
}

impl<'a> UnresolvedXs<'a> {
    pub fn new(urr: &'a Unresolved) -> Self {
        Self { urr }
    }

    /// エネルギー e (eV) での平均断面積。パラメータは lin-lin で補間する。
    pub fn cross_sections(&self, e: f64) -> Sigma {
        let _ = (self.urr, e);
        todo!()
    }
}

/// 平均断面積を計算する格子 [el, eh]。評価者のエネルギー（全系列の和）に、隣の点との比が
/// 大きい区間では NJOY と同じ補助格子の点を加える。格子点の間の断面積は、`interpolation` で補間する。
pub fn unresolved_grid(urr: &Unresolved, el: f64, eh: f64) -> Vec<f64> {
    let _ = (urr, el, eh);
    todo!()
}

/// 格子点の間の断面積の補間則（全系列で同じであること）。
pub fn interpolation(urr: &Unresolved) -> Result<Interpolation, String> {
    let _ = urr;
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::endf::{UrrJValue, UrrLValue, UrrPoint};
    use crate::physics::wave_number_factor;
    use crate::special::erfc;
    use std::f64::consts::PI;

    const AWRI: f64 = 236.0;
    const AP: f64 = 0.948;

    fn point(e: f64, gno: f64) -> UrrPoint {
        UrrPoint {
            e,
            d: 20.0,
            gx: 0.0,
            gno,
            gg: 0.023,
            gf: 0.0,
        }
    }

    fn s_wave(points: Vec<UrrPoint>, int: Interpolation) -> Unresolved {
        Unresolved {
            lssf: 0,
            spi: 0.0,
            ap: AP,
            ap_energy: None,
            naps: 1,
            l_values: vec![UrrLValue {
                awri: AWRI,
                l: 0,
                j_values: vec![UrrJValue {
                    aj: 0.5,
                    int,
                    amux: 0,
                    amun: 1,
                    amuf: 0,
                    points,
                }],
            }],
        }
    }

    /// 平均 1 の χ²_1 での E[1/(x+c)] = c^(-1/2) e^(c/2) √(π/2) erfc(√(c/2))。
    fn inverse_moment(c: f64) -> f64 {
        (PI / (2.0 * c)).sqrt() * (0.5 * c).exp() * erfc((0.5 * c).sqrt())
    }

    #[test]
    fn s_wave_without_fission_matches_closed_form() {
        // GNO は 5e3 eV で 1e-3、2e4 eV で 4e-3 なので、1e4 eV では lin-lin で 2e-3
        let urr = s_wave(
            vec![point(5.0e3, 1.0e-3), point(2.0e4, 4.0e-3)],
            Interpolation::LinLin,
        );
        let e = 1.0e4;
        let sigma = UnresolvedXs::new(&urr).cross_sections(e);

        let k = wave_number_factor() * AWRI / (AWRI + 1.0) * e.sqrt();
        let (d, gg) = (20.0, 0.023);
        let gn = 2.0e-3 * e.sqrt();
        // 中性子幅 gn x、分母 gn x + gg の期待値（c = gg/gn）
        let c = gg / gn;
        let e_inv = inverse_moment(c);
        let e_x = (1.0 - c * e_inv) / gn;
        let e_x2 = (1.0 - c + c * c * e_inv) / gn;
        let scale = 2.0 * PI * PI / (k * k) / d;
        let phase = (k * AP).sin().powi(2);
        let capture = scale * gn * gg * e_x;
        let elastic =
            scale * gn * gn * e_x2 - scale * 2.0 * gn * phase + 4.0 * PI * phase / (k * k);
        for (actual, expected) in [(sigma.capture, capture), (sigma.elastic, elastic)] {
            assert!(
                (actual - expected).abs() <= 1e-8 * expected,
                "{actual} vs {expected}"
            );
        }
        assert_eq!(sigma.fission, 0.0);
        assert!((sigma.total - (sigma.elastic + sigma.capture)).abs() <= 1e-12 * sigma.total);
    }

    #[test]
    fn grid_fills_wide_steps_like_njoy() {
        let urr = s_wave(
            vec![
                point(1.5e3, 1.0e-3),
                point(3.0e3, 1.0e-3),
                point(1.0e4, 1.0e-3),
                point(1.0e5, 1.0e-3),
            ],
            Interpolation::LinLin,
        );
        let expected = vec![
            1.5e3, 1.7e3, 2.0e3, 2.5e3, 3.0e3, 3.5e3, 4.0e3, 5.0e3, 6.0e3, 7.2e3, 8.5e3, 1.0e4,
            1.25e4, 1.5e4, 1.7e4, 2.0e4, 2.5e4, 3.0e4, 3.5e4, 4.0e4, 5.0e4, 6.0e4, 7.2e4, 8.5e4,
            1.0e5,
        ];
        assert_eq!(unresolved_grid(&urr, 1.5e3, 1.0e5), expected);
    }

    #[test]
    fn interpolation_must_agree_across_sequences() {
        let mut urr = s_wave(
            vec![point(1.5e3, 1.0e-3), point(1.0e5, 1.0e-3)],
            Interpolation::LogLog,
        );
        assert_eq!(interpolation(&urr), Ok(Interpolation::LogLog));
        let mut other = urr.l_values[0].j_values[0].clone();
        other.int = Interpolation::LinLin;
        urr.l_values[0].j_values.push(other);
        assert!(interpolation(&urr).is_err());
    }
}
