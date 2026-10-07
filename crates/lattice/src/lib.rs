//! 格子計算: 多群ライブラリから、ピンセルの輸送計算（MOC）で少数群定数を作る。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

pub mod cmfd;
pub mod cp;
pub mod equivalent;
pub mod fewgroup;
pub mod geometry;
pub mod input;
mod lu;
pub mod macroscopic;
pub mod material;
pub mod moc;
pub mod power;
pub mod shielding;
pub mod sph;
pub mod symmetry;
pub mod tracking;
pub mod ufg;
