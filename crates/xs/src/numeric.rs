//! 数値計算（特殊関数、適応積分、級数）のエラー。

/// 特殊関数や数値積分が値を返せない理由。
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum NumericError {
    #[error("{0} is outside the domain")]
    Domain(String),
    #[error("{0} did not converge")]
    NotConverged(String),
}
