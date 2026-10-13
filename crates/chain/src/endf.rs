//! 燃焼鎖に使う ENDF の評価の読み込み（decay の MF8/MT457、中性子の MF3 の MT と MF8 の生成物の状態、nfy の MT454）。
//!
//! 行とレコードは xlc-xs の `endf::Tape` で読む（ADR-0166）。

use anyhow::Result;
use xlc_xs::endf::Tape;

/// 評価の対象（MF1/MT451）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Material {
    pub mat: i32,
    pub z: u32,
    pub a: u32,
    /// 準安定状態の番号（LISO）
    pub state: u32,
    /// 励起エネルギー (eV)（ELIS）
    pub excitation: f64,
}

/// 崩壊の 1 段。RTYP の 1 桁にあたる。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Step {
    BetaMinus,
    EcBetaPlus,
    IsomericTransition,
    Alpha,
    Neutron,
    SpontaneousFission,
    Proton,
}

/// 崩壊の型の 1 つ（MF8/MT457 の RTYP、RFS、BR）。
#[derive(Debug, Clone, PartialEq)]
pub struct DecayMode {
    pub steps: Vec<Step>,
    /// 娘の準安定状態の番号（RFS）
    pub final_state: u32,
    pub ratio: f64,
}

/// MF8/MT457 の崩壊のデータ。
#[derive(Debug, Clone, PartialEq)]
pub struct DecayData {
    /// NST = 1
    pub stable: bool,
    /// 半減期 (s)。安定なら 0
    pub half_life: f64,
    /// 平均の崩壊エネルギー (eV): 軽い粒子、電磁、重い粒子
    pub energy: [f64; 3],
    pub modes: Vec<DecayMode>,
}

/// decay のサブライブラリの 1 評価。
#[derive(Debug, Clone, PartialEq)]
pub struct DecayEvaluation {
    pub material: Material,
    pub decay: DecayData,
}

/// 中性子の反応の生成物の 1 状態（MF8 の MT = 反応の、LMF = 9 か 10 の項目）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FinalStateEntry {
    /// 生成物の ZA（1000 Z + A）
    pub za: u32,
    /// 励起エネルギー (eV)（ELFS）
    pub excitation: f64,
    /// 割合を与えるファイル（LMF）
    pub mf: u32,
    pub lfs: u32,
}

/// 中性子のサブライブラリの 1 評価のうち、鎖に使うもの。
#[derive(Debug, Clone, PartialEq)]
pub struct NeutronEvaluation {
    pub material: Material,
    /// MF3 の MT（昇順）
    pub mts: Vec<u32>,
    /// MF8 が生成物の状態を MF9 か MF10 で挙げる反応の MT と、その状態（MT の昇順）
    pub final_states: Vec<(u32, Vec<FinalStateEntry>)>,
}

/// 核分裂の 1 つの生成物の独立収率。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RawYield {
    pub za: u32,
    /// 準安定状態の番号（FPS）
    pub state: u32,
    pub value: f64,
}

/// nfy の 1 評価（MT454 の独立収率）。
#[derive(Debug, Clone, PartialEq)]
pub struct YieldEvaluation {
    pub material: Material,
    /// 入射エネルギー (eV) ごとの収率（テープの順）
    pub tables: Vec<(f64, Vec<RawYield>)>,
}

/// MF1/MT451 から評価の対象を読む。テープの最初の物質だけを読む。
pub fn read_material(_tape: &Tape) -> Result<Material> {
    anyhow::bail!("not read")
}

/// decay の評価のテープを読む。
pub fn read_decay(_text: &str) -> Result<DecayEvaluation> {
    anyhow::bail!("not read")
}

/// 中性子の評価のテープを読む。`reaction_mts` の MF8 だけ生成物の状態を読む。
pub fn read_neutron(_text: &str, _reaction_mts: &[u32]) -> Result<NeutronEvaluation> {
    anyhow::bail!("not read")
}

/// nfy の評価のテープを読む。
pub fn read_yields(_text: &str) -> Result<YieldEvaluation> {
    anyhow::bail!("not read")
}

#[cfg(test)]
mod tests;
