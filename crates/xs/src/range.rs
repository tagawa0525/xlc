//! `xs range` が書き、`xs nuclide` と `xs assemble` が読む、ポイントワイズの範囲のファイル（ADR-0138）。

use anyhow::anyhow;
use serde::Deserialize;

/// `to_toml` と `from_toml` が書く、ポイントワイズの範囲のファイルの中身。
#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct RangeFile {
    /// 範囲がなければ空の配列
    pointwise_range: Vec<f64>,
}

/// ポイントワイズの範囲 (eV) を TOML に書く。範囲がなければ空の配列を書く（項目は欠かさない。ADR-0129）。
/// 浮動小数点数は、読み直して同じ値になる最短の表記（`{:?}`）で書く。
pub fn to_toml(range: Option<(f64, f64)>) -> String {
    match range {
        Some((low, high)) => format!("pointwise-range = [{low:?}, {high:?}]\n"),
        None => "pointwise-range = []\n".into(),
    }
}

/// `to_toml` の形を読む。
pub fn from_toml(text: &str) -> anyhow::Result<Option<(f64, f64)>> {
    let file: RangeFile = toml::from_str(text)?;
    match file.pointwise_range[..] {
        [] => Ok(None),
        [low, high] if low > 0.0 && high > low && high.is_finite() => Ok(Some((low, high))),
        ref other => Err(anyhow!(
            "pointwise-range must be empty or [low, high] with 0 < low < high, got {other:?}"
        )),
    }
}

#[cfg(test)]
mod tests;
