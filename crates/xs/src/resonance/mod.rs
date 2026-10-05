//! 共鳴パラメータからの断面積の計算。

mod mlbw;
mod reich_moore;

pub use mlbw::MlbwXs;
pub use reich_moore::ReichMooreXs;

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
