//! 共鳴パラメータからの断面積の計算。

mod factors;
mod fluctuation;
mod line_shape;
mod mlbw;
mod quadrature;
mod reich_moore;
mod unresolved;

pub use mlbw::MlbwXs;
pub use reich_moore::ReichMooreXs;
pub use unresolved::{
    ShieldedSigma, UnresolvedXs, UrrError, interpolation as unresolved_interpolation,
    unresolved_grid,
};

/// エネルギーごとに共鳴領域の断面積を計算できるもの。
pub trait CrossSections {
    fn cross_sections(&self, e: f64) -> Sigma;
}

impl CrossSections for ReichMooreXs {
    fn cross_sections(&self, e: f64) -> Sigma {
        ReichMooreXs::cross_sections(self, e)
    }
}

impl CrossSections for MlbwXs {
    fn cross_sections(&self, e: f64) -> Sigma {
        MlbwXs::cross_sections(self, e)
    }
}

/// 共鳴領域の断面積 (barn)。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Sigma {
    pub total: f64,
    pub elastic: f64,
    pub fission: f64,
    pub capture: f64,
}

#[cfg(test)]
pub(crate) use reich_moore::tests as reich_moore_tests;
