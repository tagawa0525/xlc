//! ADR-0164 の規則で、読んだ評価から全核種の燃焼チェーンを組み立てる。

use anyhow::Result;
use xlc_depletion::chain::Chain;

use crate::endf::{DecayEvaluation, NeutronEvaluation, YieldEvaluation};

/// チェーンの反応の MT（代表の MT。MF8 の生成物の状態もこの MT で読む）
pub const REACTION_MTS: [u32; 6] = [16, 17, 37, 102, 103, 107];

/// 組み立てたチェーンと、データにないものを規則で補った箇所の警告。
#[derive(Debug, Clone, PartialEq)]
pub struct Built {
    pub chain: Chain,
    pub warnings: Vec<String>,
}

/// 評価からチェーンを組み立てる。評価の順は問わない。
pub fn build(
    _decays: &[DecayEvaluation],
    _neutrons: &[NeutronEvaluation],
    _yields: &[YieldEvaluation],
) -> Result<Built> {
    anyhow::bail!("not built")
}

#[cfg(test)]
mod tests;
