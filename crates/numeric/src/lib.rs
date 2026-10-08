//! 段階の間で共有する計算（ADR-0002、ADR-0003）。
//!
//! 段階のクレート（xs, lattice, core）は互いに依存しないので、2つ以上の段階が使う計算はここに置く。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

pub mod anm;
pub mod nem;
pub mod quadrature;
