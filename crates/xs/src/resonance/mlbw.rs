//! 多準位 Breit-Wigner（MLBW）形式の断面積（0 K）。NJOY2016 reconr の csmlbw と同じ手順。

use super::Sigma;
use super::breit_wigner::{BreitWignerError, Levels};
use crate::endf::BreitWigner;

/// 断面積の計算に使う形に前処理した MLBW パラメータ。
#[derive(Debug, Clone)]
pub struct MlbwXs {
    levels: Levels,
}

impl MlbwXs {
    /// 波数とチャネル半径は、l ごとの AWRI から求める。
    pub fn new(data: &BreitWigner) -> Self {
        Self {
            levels: Levels::new(data),
        }
    }

    /// エネルギー e (eV) での断面積。散乱半径の表（NRO=1）が e を覆わなければエラー。
    pub fn cross_sections(&self, e: f64) -> Result<Sigma, BreitWignerError> {
        let radius = self.levels.radius(e)?;
        let mut sum = Sigma::default();
        for lv in &self.levels.l_values {
            let f = lv.factors(e, radius);
            let pifac = f.pifac;
            // 1 - cos 2φ は低エネルギーで桁落ちするので 2 sin²φ で求める
            let cos2p = 2.0 * f.phi.sin().powi(2);
            let sin2p = (2.0 * f.phi).sin();

            let fl = f64::from(lv.l.number());
            let spi = self.levels.spi;
            let ajmin = ((spi - fl).abs() - 0.5).abs();
            let nj = (spi + fl + 0.5 - ajmin + 1.0).round() as usize;
            let gj: Vec<f64> = (0..nj)
                .map(|j| (2.0 * (ajmin + j as f64) + 1.0) / (4.0 * spi + 2.0))
                .collect();
            // 2つのチャネルスピンから作れる J は gJ を1回しか数えないので、その分の硬球散乱を足す
            let missing = 2.0 * fl + 1.0 - gj.iter().sum::<f64>();

            let mut sig = vec![(0.0, 0.0); nj];
            let (mut capture, mut fission) = (0.0, 0.0);
            for r in &lv.resonances {
                let j = (r.aj - ajmin).round() as usize;
                let (shifted, gn) = r.at(&f);
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
