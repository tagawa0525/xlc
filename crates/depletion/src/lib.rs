//! 燃焼のソルバ（ADR-0161、ADR-0162）。
//!
//! 燃焼行列 A と初期の数密度 n(0) から、時間 t の後の数密度 n(t) = exp(A t) n(0) を CRAM で求める。lattice の燃焼と
//! core のミクロ燃焼が使う。時間の進め方（予測子・修正子）と反応率の作り方は、呼ぶ側が持つ。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

mod cram;
mod matrix;

pub use cram::{CramOrder, deplete};
pub use matrix::DepletionMatrix;

/// 燃焼の計算の入力が不正か、解けない。
#[derive(Debug, thiserror::Error)]
pub enum DepletionError {
    #[error("depletion matrix entry ({row}, {col}) is outside a {n} x {n} matrix")]
    OutOfRange { row: usize, col: usize, n: usize },
    #[error("depletion matrix entry ({row}, {col}) is not finite: {value}")]
    NonFiniteEntry { row: usize, col: usize, value: f64 },
    #[error("nuclide {col} produces nuclide {row} at a negative rate {value}")]
    NegativeProduction { row: usize, col: usize, value: f64 },
    #[error("nuclide {nuclide} has a positive diagonal {value} (a negative loss rate)")]
    PositiveLoss { nuclide: usize, value: f64 },
    #[error("{found} initial densities for a {expected}-nuclide matrix")]
    LengthMismatch { expected: usize, found: usize },
    #[error("initial density of nuclide {nuclide} is not finite: {value}")]
    NonFiniteDensity { nuclide: usize, value: f64 },
    #[error("depletion time must be finite and non-negative: {0}")]
    InvalidTime(f64),
    #[error("depletion matrix times the time step overflows at entry ({row}, {col})")]
    Overflow { row: usize, col: usize },
    #[error("sparse LU of the CRAM system failed: {0}")]
    Factorization(String),
    #[error("CRAM gave a non-finite density for nuclide {0}")]
    NonFiniteResult(usize),
}
