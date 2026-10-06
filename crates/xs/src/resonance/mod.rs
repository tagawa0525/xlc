//! 共鳴パラメータからの断面積の計算。

mod complex;
mod factors;
mod fluctuation;
mod line_shape;
mod mlbw;
mod quadrature;
mod r_matrix_limited;
mod reich_moore;
mod unresolved;

pub use mlbw::MlbwXs;
pub use r_matrix_limited::{RMatrixLimitedXs, RmlError};
pub use reich_moore::ReichMooreXs;
pub use unresolved::{
    ShieldedSigma, UnresolvedXs, UrrError, interpolation as unresolved_interpolation,
    unresolved_grid,
};

/// 共鳴領域の断面積 (barn)。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Sigma {
    pub total: f64,
    pub elastic: f64,
    pub fission: f64,
    pub capture: f64,
}

#[cfg(test)]
pub(crate) use r_matrix_limited::tests as r_matrix_limited_tests;
#[cfg(test)]
pub(crate) use reich_moore::tests as reich_moore_tests;
