//! 履歴の補正の残り R（ADR-0215）。公称の燃焼のテーブルと、履歴の分岐の燃焼の履歴を自分の条件で解いた分岐のテーブル（`[history]`）
//! から、燃焼度の点 u ごとに
//!
//! R(u) = f(履歴の分岐の組成、条件 j) − f(公称の組成、条件 j) − Σᵢ cᵢ(条件 j、u) (N̄ᵢ^j(u) − N̄ᵢ^0(u))
//!
//! を作り、公称のテーブルの `histories` に書き足す。第 2 項と係数 cᵢ は公称のテーブルの条件 j の状態の値で、cᵢ は core の線形の補正と
//! 同じ（数密度の微分を持つ核種は ∂Σ/∂N̄ᵢ、ほかは σᵢ。吸収は全断面積の係数から散乱の行の和を引いたもの。ADR-0216、ADR-0221）。
//! N̄ᵢ^j と N̄ᵢ^0 はそれぞれのテーブルの均質化した数密度。補正を持たない項目（χ、拡散係数、ADF）の R は 2 つの定数の差。
//!
//! 使い方: lattice histories <入力ファイル>（ファイルのパスは入力ファイルのあるディレクトリから）
//!
//! ```toml
//! nominal = "table-nominal.msgpack"            # 公称のテーブル（燃焼度の軸、公称の値、ミクロ燃焼のデータを持つ）
//! output = "table-nominal-histories.msgpack"   # R を書き足した公称のテーブル
//!
//! [[history]]
//! variable = "moderator-density"               # 履歴の変数
//! branches = { out = "table-density-branch.msgpack" }   # 制御棒の値ごとの、履歴の分岐の自分の条件のテーブル
//! ```
//!
//! 履歴の分岐のテーブルは、燃焼度の軸と、点が 1 つの分岐の軸（制御棒の軸は持てば値が 1 つ）を持ち、その状態（条件 j）が公称の
//! テーブルの格子点で、変数の軸（減速材の密度は `moderator`）だけが公称と違うこと。燃焼度の点は、公称のテーブルが条件 j の状態を
//! 解いた燃焼度の点と同じこと。変数の値は、公称の燃焼が公称の格子点の値 H⁰、履歴の分岐が条件 j の格子点の値 Hᵇ。

use std::path::{Path, PathBuf};

use anyhow::bail;
use xlc_format::fewgroup::FewGroup;
use xlc_format::table::{HistoryCorrection, HistoryRod, MicroNuclide, Table};

/// `lattice histories` の入力。
#[derive(Debug, Clone, PartialEq)]
pub struct ResidualInput {
    /// 公称のテーブル
    pub nominal: PathBuf,
    /// R を書き足した公称のテーブルの出力
    pub output: PathBuf,
    /// 履歴の変数ごとの、制御棒の値と履歴の分岐の自分の条件のテーブル
    pub histories: Vec<HistoryEntry>,
}

/// 1 つの履歴の変数の入力。
#[derive(Debug, Clone, PartialEq)]
pub struct HistoryEntry {
    pub variable: String,
    /// (制御棒の値, 履歴の分岐の自分の条件のテーブル)
    pub branches: Vec<(String, PathBuf)>,
}

impl ResidualInput {
    /// TOML の入力を読む
    pub fn parse(_text: &str) -> anyhow::Result<Self> {
        bail!("not implemented")
    }
}

/// 公称の状態の少数群定数 `state` に、線形の補正 Σᵢ cᵢ Δᵢ を足したもの（core の `micro::correct` と同じ式）。`micro` は追う核種の
/// 係数の元（数密度の微分を持てばその微分、なければ微視的断面積）、`delta` は核種ごとの数密度の差 (1/(barn cm))
pub fn linear(state: &FewGroup, _micro: &[MicroNuclide], _delta: &[f64]) -> FewGroup {
    state.clone()
}

/// 公称のテーブル `nominal` の、変数 `variable` の履歴の補正を、制御棒の値ごとの履歴の分岐の自分の条件のテーブル `branches` から作る
pub fn history_correction(
    _nominal: &Table,
    _variable: &str,
    _branches: &[(String, Table)],
) -> anyhow::Result<HistoryCorrection> {
    bail!("not implemented")
}

/// 公称のテーブル `nominal` の制御棒の値 `rod` の残りを、履歴の分岐の自分の条件のテーブル `branch` から作る。返すのは (H⁰, Hᵇ, 残り)
pub fn residuals(
    _nominal: &Table,
    _variable: &str,
    _rod: &str,
    _branch: &Table,
) -> anyhow::Result<(f64, f64, HistoryRod)> {
    bail!("not implemented")
}

/// 入力ファイル `path` を読み、R を書き足した公称のテーブルを書く
pub fn run(_path: &Path) -> anyhow::Result<()> {
    bail!("not implemented")
}

#[cfg(test)]
mod tests;
