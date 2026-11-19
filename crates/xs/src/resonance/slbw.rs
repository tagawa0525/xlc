//! 単準位 Breit-Wigner（SLBW）形式の断面積（0 K）。ENDF-102 の付録 D.1.1 の式で、NJOY2016 reconr の csslbw と同じ手順。

use super::Sigma;
use super::breit_wigner::{BreitWignerError, Levels};
use crate::endf::BreitWigner;

/// 断面積の計算に使う形に前処理した SLBW パラメータ。
#[derive(Debug, Clone)]
pub struct SlbwXs {
    levels: Levels,
}

impl SlbwXs {
    /// 波数とチャネル半径は、l ごとの AWRI から求める。
    pub fn new(data: &BreitWigner) -> Self {
        Self {
            levels: Levels::new(data),
        }
    }

    /// エネルギー e (eV) での断面積。弾性散乱は、準位の間の干渉を持たない式のとおり負になりうる（負の扱いは
    /// MF3 のバックグラウンドと足すときに決める）。散乱半径の表（NRO=1）が e を覆わなければエラー。
    pub fn cross_sections(&self, e: f64) -> Result<Sigma, BreitWignerError> {
        let radius = self.levels.radius(e)?;
        let spin_states = 2.0 * (2.0 * self.levels.spi + 1.0);
        let mut sum = Sigma::default();
        for lv in &self.levels.l_values {
            let f = lv.factors(e, radius);
            let sin2 = f.phi.sin().powi(2);
            let sin2p = (2.0 * f.phi).sin();
            // 硬球散乱は (2l+1) の全部の状態で数える（MLBW のような足りない J の補いは要らない）
            let fl = f64::from(lv.l.number());
            let (mut elastic, mut capture, mut fission) = (4.0 * (2.0 * fl + 1.0) * sin2, 0.0, 0.0);
            for r in &lv.resonances {
                let (shifted, gn) = r.at(&f);
                let width = gn + r.gg + r.gf;
                let d = e - shifted;
                let g = (2.0 * r.aj.abs() + 1.0) / spin_states;
                // (D.2) の Γn² − 2 Γn Γ sin²φ + 2 (E − E'r) Γn sin 2φ を Γn でくくる
                let common = g * gn / (d * d + width * width / 4.0);
                elastic += common * (gn - 2.0 * width * sin2 + 2.0 * d * sin2p);
                capture += common * r.gg;
                fission += common * r.gf;
            }
            sum.elastic += f.pifac * elastic;
            sum.capture += f.pifac * capture;
            sum.fission += f.pifac * fission;
        }
        sum.total = sum.elastic + sum.fission + sum.capture;
        Ok(sum)
    }
}

#[cfg(test)]
pub(crate) mod tests;
