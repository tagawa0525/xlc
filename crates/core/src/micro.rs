//! ノードのミクロ燃焼（ADR-0216）。
//!
//! 追う核種（テーブルの `micro` の `nuclides`）と、補助の核種（追う核種へチェーンの崩壊で流れ込む核種のうち、追う核種の反応、
//! 核分裂の収率、補助の核種の崩壊で生まれうるもの）を、格子計算の燃焼と同じチェーンで解く。追う核種の反応率は、テーブルの
//! 微視的断面積 × ノードの少数群の中性子束で、中性子束の大きさはノードの出力密度で決める（追う核種の核分裂率 × チェーンの Q 値の和が
//! 出力密度になる。格子計算の燃焼と同じ規格化。ADR-0173）。補助の核種は崩壊だけ。核分裂の収率の重みは、テーブルの収率の境界より下の
//! 核分裂の割合から ADR-0183 の点の選び方で作る。
//!
//! ノードの少数群定数には線形の補正 Σᵢ σᵢ (Nᵢ − N̄ᵢ⁰) を足す（全断面積、吸収、νΣf、Σf、散乱行列。ADR-0216）。

use xlc_depletion::DepletionMatrix;
use xlc_depletion::chain::Chain;
use xlc_format::fewgroup::FewGroup;
use xlc_format::table::MicroNuclide;

use crate::state::MicroAt;

/// 追う核種と補助の核種の燃焼行列の材料（チェーンから 1 度だけ作る）
#[derive(Debug, Clone, PartialEq)]
pub struct MicroChain {
    /// 追う核種の名前（テーブルの `nuclides` の順）
    pub nuclides: Vec<String>,
    /// 解く核種（追う核種と補助の核種）のチェーンの名前（チェーンの順）
    pub solved: Vec<String>,
}

/// 追う核種 1 つの 1 原子あたりの反応率 (1/s)
#[derive(Debug, Clone, PartialEq)]
pub struct TrackedRates {
    /// [チェーンの反応][経路]
    pub reactions: Vec<Vec<f64>>,
    pub fission: f64,
    /// チェーンの収率のエネルギーの点ごとの重み（核分裂しない核種は空）
    pub yield_weights: Vec<f64>,
}

impl MicroChain {
    /// チェーン `chain` から、追う核種 `nuclides` とそのテーブルの微視的断面積 `reference`（反応の並びと核分裂の有無を見る）の
    /// 燃焼行列の材料を作る。追う核種がチェーンにないか、反応の並びや核分裂の有無がチェーンと違えばエラー
    pub fn new(
        _chain: &Chain,
        _nuclides: &[String],
        _reference: &[MicroNuclide],
    ) -> Result<Self, String> {
        Err("not implemented".into())
    }

    /// 補助の核種の数
    pub fn auxiliary(&self) -> usize {
        0
    }

    /// 追う核種の数密度 `tracked`（テーブルの順）から、解く核種の数密度（補助の核種は 0）
    pub fn initial(&self, _tracked: &[f64]) -> Vec<f64> {
        Vec::new()
    }

    /// 解く核種の数密度 `n` のうち、追う核種のもの（テーブルの順）
    pub fn tracked(&self, _n: &[f64]) -> Vec<f64> {
        Vec::new()
    }

    /// ノードの微視的断面積 `micro`、少数群の中性子束の形 `flux`、解く核種の数密度 `n` (1/(barn cm)) で、追う核種の核分裂の出力が
    /// 出力密度 `power_density` (W/cm³) になるように規格化した、追う核種ごとの反応率
    pub fn rates(
        &self,
        _micro: &MicroAt,
        _flux: &[f64],
        _n: &[f64],
        _power_density: f64,
    ) -> Result<Vec<TrackedRates>, String> {
        Err("not implemented".into())
    }

    /// 追う核種の反応率 `rates` の、解く核種の燃焼行列
    pub fn matrix(&self, _rates: &[TrackedRates]) -> Result<DepletionMatrix, String> {
        Err("not implemented".into())
    }
}

/// 収率の境界より下の核分裂の割合 `thermal` の、収率のエネルギーの点 `energies` ごとの重み（ADR-0183。格子計算の燃焼の
/// `yield_weights` と同じ点の選び方）
pub fn yield_weights(energies: &[f64], _thermal: f64) -> Vec<f64> {
    vec![0.0; energies.len()]
}

/// 少数群定数 `f` に、線形の補正 Σᵢ σᵢ (Nᵢ − N̄ᵢ⁰) を足す（`n` は追う核種の数密度、テーブルの順）
pub fn correct(_f: &mut FewGroup, _micro: &MicroAt, _n: &[f64]) {}

#[cfg(test)]
mod tests;
