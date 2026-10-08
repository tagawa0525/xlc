//! ENDF 形式の評価済み核データから多群断面積ライブラリを作る。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

pub mod doppler;
pub mod endf;
pub mod multigroup;
pub mod numeric;
pub mod physics;
pub mod pointwise;
pub mod process;
pub mod resonance;
pub mod special;
