//! 単準位 Breit-Wigner（SLBW）形式の断面積（0 K）。

use super::Sigma;
use super::breit_wigner::{BreitWignerError, Levels};
use crate::endf::BreitWigner;

/// 断面積の計算に使う形に前処理した SLBW パラメータ。
#[derive(Debug, Clone)]
pub struct SlbwXs {
    levels: Levels,
}

impl SlbwXs {
    pub fn new(data: &BreitWigner) -> Self {
        Self {
            levels: Levels::new(data),
        }
    }

    pub fn cross_sections(&self, e: f64) -> Result<Sigma, BreitWignerError> {
        let _ = self.levels.radius(e)?;
        Ok(Sigma::default())
    }
}

#[cfg(test)]
pub(crate) mod tests;
