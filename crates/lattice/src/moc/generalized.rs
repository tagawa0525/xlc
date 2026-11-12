//! 一般化随伴の固定源の問題（ADR-0247 の段 3。#0279）。

use super::{Detected, Material, MocError, Options, Start};
use crate::tracking::Tracking;

/// 固定源の反復の収束の方法
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GeneralizedMethod {
    /// 射影つきの源の反復（加速なし）
    SourceIteration,
}

/// 一般化随伴の固定源の反復の方法と打ち切り
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeneralizedOptions {
    pub method: GeneralizedMethod,
    /// 中性子束の変化の上限（群ごとの、領域の |Γ| の最大値に対する比）
    pub tolerance: f64,
    /// スイープの回数の上限
    pub max_iterations: usize,
}

/// 1 つの源の一般化随伴
#[derive(Debug, Clone, PartialEq)]
pub struct Generalized {
    /// [領域][群] のスカラー中性子束
    pub flux: Vec<Vec<f64>>,
    /// スイープの回数
    pub iterations: usize,
}

/// 一般化随伴の固定源の問題を解く（ADR-0247）
#[allow(clippy::too_many_arguments)]
pub fn solve_generalized(
    _tracking: &Tracking,
    _materials: &[&Material],
    _options: &Options,
    _symmetry: Option<&Detected>,
    _forward: &Start,
    _adjoint: &Start,
    _sources: &[Vec<Vec<f64>>],
    _settings: &GeneralizedOptions,
) -> Result<Vec<Generalized>, MocError> {
    Err(MocError::GeneralizedNotConverged { iterations: 0 })
}

#[cfg(test)]
mod tests;
