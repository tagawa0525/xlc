//! ENDF/B の decay、nfy、中性子のサブライブラリから全核種の燃焼鎖を作る（ADR-0164、ADR-0166）。
//!
//! `endf` が評価を読み、`build` が ADR-0164 の規則で鎖を組み立てる。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

pub mod build;
pub mod elements;
pub mod endf;
