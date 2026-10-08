//! 減速計算の格子（xs のライブラリの重みと、lattice のピンごとの超微細群の計算で共通）。

/// 1つの窓（散乱で届くレサジーの幅 Δ）に入れる格子点の最少の数。格子はポイントワイズの断面積の格子に、これより
/// 粗いところだけ点を足す（xs の `finer_grids_do_not_change_the_flux` で収束を確かめた値）
pub const WINDOW_POINTS: usize = 16;
