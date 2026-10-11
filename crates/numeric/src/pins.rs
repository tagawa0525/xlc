//! ピン出力の再構成のノードの中の均質な中性子束（ADR-0096）。core の再構成と、lattice の角の不連続因子
//! （ADR-0125）の両方で使う。

mod basis;
pub mod corners;
pub mod node;

pub use basis::FACE_POINTS;
