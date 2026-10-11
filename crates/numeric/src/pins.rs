//! ピン出力の再構成のノードの中の均質な中性子束（ADR-0096）。core の再構成と、lattice の角の不連続因子
//! （ADR-0131）の両方で使う。

mod basis;
pub mod corners;
pub mod node;

pub use basis::FACE_POINTS;

/// セルの平均の Gauss-Legendre の点の数（各方向）。ピンのセル（1.26 cm）の中で κh ≲ 1 なので、6 点の誤差は 1e-12 より
/// 小さい。lattice の形状関数の分母と core の再構成で同じ点を使わないと、定数を作った体系でセルの出力を再現しない
pub const CELL_POINTS: usize = 6;
