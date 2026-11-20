//! 共鳴パラメータからの断面積の計算。

mod breit_wigner;
mod coulomb;
mod factors;
mod fluctuation;
mod ladder;
mod linear;
mod mlbw;
mod quadrature;
mod r_matrix_limited;
mod reich_moore;
mod slbw;
mod unresolved;

pub use breit_wigner::BreitWignerError;
pub use ladder::{Background, Bondarenko, Sampling, Scales};
pub use mlbw::MlbwXs;
pub use r_matrix_limited::{RMatrixLimitedXs, RmlError};
pub use reich_moore::{ReichMooreError, ReichMooreXs};
pub use slbw::SlbwXs;
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
#[cfg(test)]
pub(crate) use slbw::tests as slbw_tests;
