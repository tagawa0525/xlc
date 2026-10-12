//! 炉心計算: 少数群定数から、3次元の拡散方程式をノード展開法（NEM）で解く。
//!
//! パッケージ名は `xlc-core` とする（`core` は Rust 標準のクレートと衝突するため）。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

pub mod case;
pub mod input;
pub mod nem;
pub mod pins;
pub mod search;
pub mod state;
