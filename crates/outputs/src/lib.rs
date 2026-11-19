//! ステージの出力を、それを deps にするステージの実行ファイルと同じ読み手で読み、結果を捨てる（形式を変える修正の検査。
//! ADR-0258）。計算はしない。

use std::str::FromStr;

use anyhow::{Context, anyhow, bail};
use serde::Deserialize;
use xlc_core::depletion::CoreHistory;
use xlc_depletion::chain::Chain;
use xlc_format::library::Library;
use xlc_format::medium::Medium;
use xlc_format::table::Table;
use xlc_lattice::depletion::History;

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
pub fn read(kind: Kind, bytes: &[u8]) -> anyhow::Result<()> {
    match kind {
        Kind::Msgpack => read_msgpack(bytes),
        Kind::XsRange => {
            let text = std::str::from_utf8(bytes).context("the range file is not UTF-8")?;
            xlc_xs::range::from_toml(text).map(|_| ())
        }
    }
}

/// 形式の名前だけを読む。ほかの項目は読み飛ばす（読み手が全体を読むので、ここで木を作らない）
#[derive(Deserialize)]
struct Head {
    format: String,
}

/// 名前の一致は読み手が確かめ直すので、ここでは読み手を選ぶだけにする。名前を変えた形式の古い出力は、どの読み手にも
/// 当たらずに止まる
fn read_msgpack(bytes: &[u8]) -> anyhow::Result<()> {
    let head: Head = rmp_serde::from_slice(bytes).context("cannot read the format name")?;
    let format = head.format.as_str();
    if format == xlc_format::library::FORMAT {
        Library::from_msgpack(bytes)?;
    } else if format == xlc_format::table::FORMAT {
        Table::from_msgpack(bytes)?;
    } else if format == xlc_format::medium::FORMAT {
        Medium::from_msgpack(bytes)?;
    } else if format == xlc_depletion::chain::FORMAT {
        Chain::from_msgpack(bytes)?;
    } else if format == xlc_lattice::depletion::FORMAT {
        History::from_msgpack(bytes)?;
    } else if format == xlc_core::depletion::FORMAT {
        CoreHistory::from_msgpack(bytes)?;
    } else {
        bail!("no Rust reader takes the format {format:?}");
    }
    Ok(())
}
