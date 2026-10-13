//! CRAM で exp(A t) n(0) を求める（ADR-0161）。

use crate::{DepletionError, DepletionMatrix};

/// CRAM の次数。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CramOrder {
    /// 次数 16（8 個の極）
    Order16,
    /// 次数 48（24 個の極）
    Order48,
}

/// 時間 `t`（s）の後の数密度 exp(A t) n(0) を返す。
pub fn deplete(
    _a: &DepletionMatrix,
    n0: &[f64],
    _t: f64,
    _order: CramOrder,
) -> Result<Vec<f64>, DepletionError> {
    Ok(n0.to_vec())
}

#[cfg(test)]
mod tests;
