//! 段階の間で受け渡すデータファイルの形式（全体計画の「ファイルの約束事」）。
//!
//! 段階のクレート（xs, lattice, core）は互いに依存せず、データファイルの型はここに置く。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

pub mod group;
pub mod library;
pub mod msgpack;
