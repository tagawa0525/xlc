//! ピン出力（セルごとの核分裂率）のファイル。`lattice` の `pin-powers` と、`core` の再構成（ADR-0096）が書き、
//! 検証のツールが読む。
//!
//! フィールド名付きの map 形式の MessagePack（ADR-0011、ADR-0013）: `format` ("xlc-pin-powers") と、
//! `rows`（[行（上から）][列（左から）] の値）。
//!
//! 軸方向の層ごとの出力（`core` の `axial-pin-powers` が書き、examples/vera/openmc_model_3d.py の `layers.msgpack` と
//! 同じ形）: `format` ("xlc-axial-pin-powers")、`bounds`（層の境界の z (cm)、下から）、`axial`（層ごとの径方向に積分した
//! 出力の線密度）、`layers`（[層（下から）][行（上から）][列（左から）] の線密度）。

use serde::Serialize;

use crate::msgpack::{self, MsgpackError};

const FORMAT: &str = "xlc-pin-powers";
const AXIAL_FORMAT: &str = "xlc-axial-pin-powers";

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

/// 軸方向の層ごとの出力。
#[derive(Debug, Clone, PartialEq)]
pub struct AxialPinPowers {
    pub bounds: Vec<f64>,
    pub axial: Vec<f64>,
    pub layers: Vec<Vec<Vec<f64>>>,
}

impl AxialPinPowers {
    /// MessagePack（フィールド名付きの map 形式）で書く。
    pub fn to_msgpack(&self) -> Result<Vec<u8>, MsgpackError> {
        #[derive(Serialize)]
        struct File<'a> {
            format: &'static str,
            bounds: &'a [f64],
            axial: &'a [f64],
            layers: &'a [Vec<Vec<f64>>],
        }
        msgpack::to_named(&File {
            format: AXIAL_FORMAT,
            bounds: &self.bounds,
            axial: &self.axial,
            layers: &self.layers,
        })
    }
}

#[cfg(test)]
mod tests;
