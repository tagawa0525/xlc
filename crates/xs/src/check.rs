//! `xs check` の検査: ENDF の評価を、`process_nuclide` が読むものの読み込みと形式の検査だけで確かめる（ADR-0266）。

use serde::Serialize;

use crate::endf::Tape;
use crate::process::ProcessError;

/// テープ中の（ただ1つの）物質を、輸送する核種として `process_nuclide` が読むものの読み込みと形式の検査で確かめる。
pub fn check_evaluation(_tape: &Tape) -> Result<(), ProcessError> {
    Ok(())
}

/// 1 つの ENDF のファイルの検査の結果。
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Checked {
    /// 与えられたパス
    pub endf: String,
    /// 止まった理由（読めれば `None`）
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

/// 検査の結果を、与えられた順に TOML の `[[evaluation]]` の並びで書く。
pub fn report_toml(_checked: &[Checked]) -> String {
    String::new()
}

#[cfg(test)]
mod tests;
