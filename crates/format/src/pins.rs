//! ピン出力（セルごとの核分裂率）のファイル。`lattice` の `pin-powers` と、`core` の再構成（ADR-0096）が書き、
//! 検証の道具が読む。
//!
//! フィールド名付きの map 形式の MessagePack（ADR-0011、ADR-0013）: `format` ("xlc-pin-powers") と、
//! `rows`（[行（上から）][列（左から）] の値）。

use serde::Serialize;

use crate::msgpack::{self, MsgpackError};

const FORMAT: &str = "xlc-pin-powers";

/// セルごとの出力。行は上から、各行は左から。
#[derive(Debug, Clone, PartialEq)]
pub struct PinPowers {
    pub rows: Vec<Vec<f64>>,
}

impl PinPowers {
    /// MessagePack（フィールド名付きの map 形式）で書く。
    pub fn to_msgpack(&self) -> Result<Vec<u8>, MsgpackError> {
        #[derive(Serialize)]
        struct File<'a> {
            format: &'static str,
            rows: &'a [Vec<f64>],
        }
        msgpack::to_named(&File {
            format: FORMAT,
            rows: &self.rows,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pin_powers_are_written_by_rows() {
        #[derive(serde::Deserialize, Debug, PartialEq)]
        #[serde(deny_unknown_fields)]
        struct File {
            format: String,
            rows: Vec<Vec<f64>>,
        }
        let rows = vec![vec![1.0, 0.0], vec![0.5, 1.5]];
        let bytes = PinPowers { rows: rows.clone() }.to_msgpack().unwrap();
        let value: File = msgpack::from_value(msgpack::read_value(&bytes).unwrap()).unwrap();
        assert_eq!(
            value,
            File {
                format: "xlc-pin-powers".into(),
                rows,
            }
        );
    }
}
