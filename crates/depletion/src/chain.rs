//! 燃焼チェーン（ADR-0164、ADR-0166）。xs の段階の `chain`（crates/chain）が ENDF から作って書き、lattice と core が読んで
//! 燃焼行列を組み立てる。
//!
//! フィールド名付きの map 形式の MessagePack（ADR-0011、ADR-0013）:
//!
//! ```text
//! {format: "xlc-chain",
//!  nuclides: [{name: "U-235", z, a, state (LISO), excitation (eV),
//!              decay: {half_life (s), energy: {light, electromagnetic, heavy} (eV),
//!                      modes: [{mode: "beta-,n", ratio, products: [{nuclide (添字), count}]}]}、安定なら nil,
//!              reactions: [{mt, channels: [{final_state: {mf, lfs}（状態に分けないなら nil）, products}]}],
//!              fission: {yield_source: "U-238", energies: [eV、昇順], yields: [[{nuclide, value}]（エネルギーの順）]}
//!                       （核分裂しないなら nil）}]}
//! ```
//!
//! 核種は (z, a, state) の昇順に並び、生成物は核種の添字で指す。

use serde::{Deserialize, Serialize};
use xlc_format::msgpack::MsgpackError;

use crate::{DepletionError, DepletionMatrix};

/// チェーンの読み書きと、燃焼行列の組み立てのエラー。
#[derive(Debug, thiserror::Error)]
pub enum ChainError {
    #[error(transparent)]
    Msgpack(#[from] MsgpackError),
    #[error("expected format xlc-chain, found {0}")]
    Format(String),
    #[error("invalid chain: {0}")]
    Invalid(String),
    #[error("invalid reaction rates: {0}")]
    Rates(String),
    #[error(transparent)]
    Depletion(#[from] DepletionError),
}

/// 燃焼チェーン。核種は (z, a, state) の昇順。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Chain {
    pub nuclides: Vec<Nuclide>,
}

/// チェーンの 1 核種（1 つの状態）。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Nuclide {
    /// ADR-0015 の形（`U-235`、`Am-242m`）
    pub name: String,
    pub z: u32,
    pub a: u32,
    /// 準安定状態の番号（ENDF の LISO。基底状態は 0）
    pub state: u32,
    /// 励起エネルギー (eV)
    pub excitation: f64,
    /// 安定なら `None`（nil を書く。項目は省けない）
    #[serde(deserialize_with = "nil_or")]
    pub decay: Option<Decay>,
    pub reactions: Vec<Reaction>,
    /// 核分裂しないなら `None`（nil を書く。項目は省けない）
    #[serde(deserialize_with = "nil_or")]
    pub fission: Option<Fission>,
}

/// 崩壊。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Decay {
    /// 半減期 (s)
    pub half_life: f64,
    pub energy: DecayEnergy,
    /// 分岐比の和は 1
    pub modes: Vec<DecayMode>,
}

/// 1 回の崩壊の平均のエネルギー (eV)。
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DecayEnergy {
    /// 軽い粒子（β と電子）
    pub light: f64,
    /// 電磁（γ と X 線）
    pub electromagnetic: f64,
    /// 重い粒子（α、中性子、陽子、核分裂片）
    pub heavy: f64,
}

/// 崩壊の型の 1 つ。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DecayMode {
    /// 型を順に `,` でつないだもの（`beta-`、`ec/beta+`、`IT`、`alpha`、`n`、`sf`、`p`）
    pub mode: String,
    pub ratio: f64,
    /// 生まれる核種。自発核分裂と、娘がチェーンにない型は空（親は減る）
    pub products: Vec<Product>,
}

/// 1 回の遷移で生まれる核種と個数。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Product {
    /// 核種の添字
    pub nuclide: usize,
    pub count: u32,
}

/// 中性子の反応の 1 つ（核分裂を除く）。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Reaction {
    /// 反応を代表する MT（102、16、17、37、103、107）
    pub mt: u32,
    /// 生成物の状態ごとの経路。状態に分けない反応は 1 つ
    pub channels: Vec<Channel>,
}

/// 反応の、生成物の状態ごとの経路。反応率は経路ごとに与える。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Channel {
    /// 状態ごとの割合を与える ENDF のファイル（MF9 か MF10）と状態（LFS）。状態に分けない反応は `None`
    #[serde(deserialize_with = "nil_or")]
    pub final_state: Option<FinalState>,
    /// 生まれる核種（生成物と、(n,p) の H-1、(n,α) の He-4）。生成物がチェーンにないときは空
    pub products: Vec<Product>,
}

/// 状態ごとの割合の出所。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FinalState {
    pub mf: u32,
    pub lfs: u32,
}

/// 核分裂の収率。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Fission {
    /// 収率を取った核種（自分の収率がなければ、ADR-0164 の規則で選んだ別の核種）
    pub yield_source: String,
    /// 収率の入射エネルギー (eV)、昇順
    pub energies: Vec<f64>,
    /// エネルギーごとの独立収率
    pub yields: Vec<Vec<Yield>>,
}

/// 1 つの生成物の収率。
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Yield {
    pub nuclide: usize,
    pub value: f64,
}

/// 1 核種の反応率 (1/s)。形はチェーンの核種に合わせる（`Chain::zero_rates`）。
#[derive(Debug, Clone, PartialEq)]
pub struct NuclideRates {
    /// [反応][経路]
    pub reactions: Vec<Vec<f64>>,
    pub fission: f64,
    /// 収率のエネルギーの点ごとの重み（和は 1。`Fission::yield_bounds` の区間の核分裂率の割合）。核分裂しない核種は空
    pub yield_weights: Vec<f64>,
}

/// `Option` を、値か nil からだけ読む（serde の derive は、省いた `Option` を `None` にするため）
fn nil_or<'de, D, T>(d: D) -> Result<Option<T>, D::Error>
where
    D: serde::Deserializer<'de>,
    T: Deserialize<'de>,
{
    Option::deserialize(d)
}

impl Chain {
    /// MessagePack（フィールド名付きの map 形式）で書く。読むときと同じ検査をする。
    pub fn to_msgpack(&self) -> Result<Vec<u8>, ChainError> {
        Err(ChainError::Invalid("not written".into()))
    }

    /// MessagePack から読み、`check` で確かめる。
    pub fn from_msgpack(_bytes: &[u8]) -> Result<Self, ChainError> {
        Err(ChainError::Invalid("not read".into()))
    }

    /// 名前、並び、添字の範囲、分岐比の和、収率の形を確かめる。
    pub fn check(&self) -> Result<(), ChainError> {
        Ok(())
    }

    /// 反応率がすべて 0 の、チェーンの形に合った反応率。
    pub fn zero_rates(&self) -> Vec<NuclideRates> {
        Vec::new()
    }

    /// 燃焼行列の項（行が生まれる核種、列が親。`DepletionMatrix::from_triplets` の形）。
    pub fn matrix_entries(
        &self,
        _rates: &[NuclideRates],
    ) -> Result<Vec<(usize, usize, f64)>, ChainError> {
        Ok(Vec::new())
    }

    /// 燃焼行列。
    pub fn matrix(&self, rates: &[NuclideRates]) -> Result<DepletionMatrix, ChainError> {
        Ok(DepletionMatrix::from_triplets(
            self.nuclides.len(),
            &self.matrix_entries(rates)?,
        )?)
    }
}

impl Fission {
    /// 収率のエネルギーの点ごとの区間の、内側の境界 (eV)。隣の点の幾何平均で区切る（ADR-0164）。点が k 個なら k − 1 個
    pub fn yield_bounds(&self) -> Vec<f64> {
        Vec::new()
    }
}

#[cfg(test)]
mod tests;
