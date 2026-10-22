//! 格子計算: 多群ライブラリから、ピンセルの輸送計算（MOC）で少数群定数を作る。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

pub mod axial;
pub mod b1;
pub mod balance;
pub mod branches;
pub mod cmfd;
pub mod cp;
pub mod depletion;
pub mod equivalent;
pub mod fewgroup;
pub mod geometry;
pub mod history;
pub mod input;
pub mod macroscopic;
pub mod material;
pub mod medium;
pub mod moc;
pub mod power;
pub mod shielding;
pub mod sph;
pub mod symmetry;
pub mod track_table;
pub mod tracking;
pub mod ufg;
pub mod water;
