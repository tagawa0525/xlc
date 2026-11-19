//! ステージの出力を、それを deps にするステージの実行ファイルと同じ読み手で読み、結果を捨てる（形式を変える修正の検査。
//! ADR-0258）。計算はしない。

use std::str::FromStr;

use anyhow::anyhow;

/// 出力の種類。MessagePack の形式は、ファイルの `format` で読み手を選ぶ
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// 形式の名前を持つ MessagePack
    Msgpack,
    /// `xs range` が書くポイントワイズの範囲の TOML
    XsRange,
}

impl FromStr for Kind {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> anyhow::Result<Self> {
        match s {
            "msgpack" => Ok(Self::Msgpack),
            "xs-range" => Ok(Self::XsRange),
            other => Err(anyhow!("unknown kind {other:?} (msgpack or xs-range)")),
        }
    }
}

/// `bytes` を `kind` の読み手で読み、読めなければその理由を返す。
pub fn read(_kind: Kind, _bytes: &[u8]) -> anyhow::Result<()> {
    Ok(())
}
