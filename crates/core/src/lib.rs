//! 炉心計算: 少数群定数から、3次元の拡散方程式をノード展開法（NEM）で解く。
//!
//! パッケージ名は `xlc-core` とする（`core` は Rust 標準のクレートと衝突するため）。

pub mod input;
pub mod nem;
