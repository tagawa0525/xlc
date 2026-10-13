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
use xlc_format::library::is_nuclide_name;
use xlc_format::msgpack::{self, MsgpackError};

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

/// 形式の名前
const FORMAT: &str = "xlc-chain";

/// 分岐比の和と収率の重みの和の、1 からのずれの許容。チェーンを作るときに比で割り直す（ADR-0164）ので、和の丸め
/// （数十項で数 ε）より十分大きく、割り直しを忘れたチェーン（ずれは評価で最大 9e-7）を止める大きさ
const SUM_TOLERANCE: f64 = 1e-12;

/// ファイルの中身（フィールド名付きの map 形式で書く）
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ChainFile {
    format: String,
    nuclides: Vec<Nuclide>,
}

impl Chain {
    /// MessagePack（フィールド名付きの map 形式）で書く。読むときと同じ検査をして、壊れたチェーンを書かない。
    pub fn to_msgpack(&self) -> Result<Vec<u8>, ChainError> {
        self.check()?;
        Ok(msgpack::to_named(&ChainFile {
            format: FORMAT.into(),
            nuclides: self.nuclides.clone(),
        })?)
    }

    /// MessagePack から読み、`check` で確かめる。
    pub fn from_msgpack(bytes: &[u8]) -> Result<Self, ChainError> {
        let value = msgpack::read_value(bytes)?;
        require_maps(&value)?;
        let file: ChainFile = msgpack::from_value(value)?;
        if file.format != FORMAT {
            return Err(ChainError::Format(file.format));
        }
        let chain = Self {
            nuclides: file.nuclides,
        };
        chain.check()?;
        Ok(chain)
    }

    /// 名前、並び、添字の範囲、分岐比の和、収率の形を確かめる。
    pub fn check(&self) -> Result<(), ChainError> {
        let n = self.nuclides.len();
        let invalid = |message: String| Err(ChainError::Invalid(message));
        let products_in_range =
            |products: &[Product]| products.iter().all(|p| p.nuclide < n && p.count > 0);
        for (k, nuclide) in self.nuclides.iter().enumerate() {
            let name = &nuclide.name;
            if !is_nuclide_name(name) {
                return invalid(format!(
                    "nuclide name `{name}` is not of the form U-235, Am-242m"
                ));
            }
            if k > 0 {
                let previous = &self.nuclides[k - 1];
                if (previous.z, previous.a, previous.state) >= (nuclide.z, nuclide.a, nuclide.state)
                {
                    return invalid(format!(
                        "{name} is not after {} in (z, a, state)",
                        previous.name
                    ));
                }
            }
            if !(nuclide.excitation.is_finite() && nuclide.excitation >= 0.0) {
                return invalid(format!("{name}: invalid excitation {}", nuclide.excitation));
            }
            if let Some(decay) = &nuclide.decay {
                if !(decay.half_life.is_finite() && decay.half_life > 0.0) {
                    return invalid(format!("{name}: invalid half-life {}", decay.half_life));
                }
                let e = decay.energy;
                if ![e.light, e.electromagnetic, e.heavy]
                    .iter()
                    .all(|v| v.is_finite() && *v >= 0.0)
                {
                    return invalid(format!("{name}: invalid decay energy {e:?}"));
                }
                if decay.modes.is_empty() {
                    return invalid(format!("{name}: an unstable nuclide without decay modes"));
                }
                let mut sum = 0.0;
                for mode in &decay.modes {
                    if !(mode.ratio.is_finite() && mode.ratio > 0.0) {
                        return invalid(format!(
                            "{name}: invalid ratio {} of {}",
                            mode.ratio, mode.mode
                        ));
                    }
                    if !products_in_range(&mode.products) {
                        return invalid(format!(
                            "{name}: a product of {} is outside the chain",
                            mode.mode
                        ));
                    }
                    sum += mode.ratio;
                }
                if (sum - 1.0).abs() > SUM_TOLERANCE {
                    return invalid(format!("{name}: decay ratios sum to {sum}"));
                }
            }
            for reaction in &nuclide.reactions {
                let mt = reaction.mt;
                if reaction.channels.is_empty() {
                    return invalid(format!("{name}: MT{mt} has no channels"));
                }
                // 状態に分けた反応は、どの経路も割合の出所を持つ（反応率を経路ごとに与えるため）
                if reaction.channels.len() > 1
                    && reaction.channels.iter().any(|c| c.final_state.is_none())
                {
                    return invalid(format!(
                        "{name}: MT{mt} has several channels without final states"
                    ));
                }
                if !reaction
                    .channels
                    .iter()
                    .all(|c| products_in_range(&c.products))
                {
                    return invalid(format!("{name}: a product of MT{mt} is outside the chain"));
                }
            }
            if let Some(fission) = &nuclide.fission {
                let energies = &fission.energies;
                if energies.is_empty()
                    || !energies.iter().all(|e| e.is_finite() && *e > 0.0)
                    || !energies.windows(2).all(|w| w[0] < w[1])
                {
                    return invalid(format!(
                        "{name}: yield energies must be positive and ascending, got {energies:?}"
                    ));
                }
                if fission.yields.len() != energies.len() {
                    return invalid(format!(
                        "{name}: {} yield tables for {} energies",
                        fission.yields.len(),
                        energies.len()
                    ));
                }
                for table in &fission.yields {
                    if !table
                        .iter()
                        .all(|y| y.nuclide < n && y.value.is_finite() && y.value > 0.0)
                    {
                        return invalid(format!(
                            "{name}: a yield is outside the chain or not positive"
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    /// 反応率がすべて 0 の、チェーンの形に合った反応率。収率の重みは最初のエネルギーの点に 1 を置く。
    pub fn zero_rates(&self) -> Vec<NuclideRates> {
        self.nuclides
            .iter()
            .map(|nuclide| NuclideRates {
                reactions: nuclide
                    .reactions
                    .iter()
                    .map(|r| vec![0.0; r.channels.len()])
                    .collect(),
                fission: 0.0,
                yield_weights: nuclide.fission.as_ref().map_or_else(Vec::new, |f| {
                    let mut weights = vec![0.0; f.energies.len()];
                    weights[0] = 1.0;
                    weights
                }),
            })
            .collect()
    }

    /// 燃焼行列の項（行が生まれる核種、列が親。`DepletionMatrix::from_triplets` の形）。崩壊は ln 2 / 半減期に分岐比を、
    /// 反応は経路の反応率を、核分裂は核分裂率に、収率のエネルギーの点ごとの重みで足した収率を掛ける。
    pub fn matrix_entries(
        &self,
        rates: &[NuclideRates],
    ) -> Result<Vec<(usize, usize, f64)>, ChainError> {
        if rates.len() != self.nuclides.len() {
            return Err(ChainError::Rates(format!(
                "{} rates for {} nuclides",
                rates.len(),
                self.nuclides.len()
            )));
        }
        let mut entries = Vec::new();
        let gain = |entries: &mut Vec<(usize, usize, f64)>,
                    col: usize,
                    products: &[Product],
                    rate: f64| {
            for p in products {
                entries.push((p.nuclide, col, f64::from(p.count) * rate));
            }
        };
        for (j, (nuclide, rate)) in self.nuclides.iter().zip(rates).enumerate() {
            let name = &nuclide.name;
            let bad = |message: String| Err(ChainError::Rates(format!("{name}: {message}")));
            let valid = |r: f64| r.is_finite() && r >= 0.0;
            if let Some(decay) = &nuclide.decay {
                let lambda = std::f64::consts::LN_2 / decay.half_life;
                entries.push((j, j, -lambda));
                for mode in &decay.modes {
                    gain(&mut entries, j, &mode.products, mode.ratio * lambda);
                }
            }
            if rate.reactions.len() != nuclide.reactions.len() {
                return bad(format!(
                    "{} reaction rates for {} reactions",
                    rate.reactions.len(),
                    nuclide.reactions.len()
                ));
            }
            for (reaction, channel_rates) in nuclide.reactions.iter().zip(&rate.reactions) {
                if channel_rates.len() != reaction.channels.len() {
                    return bad(format!(
                        "{} channel rates for the {} channels of MT{}",
                        channel_rates.len(),
                        reaction.channels.len(),
                        reaction.mt
                    ));
                }
                for (channel, &r) in reaction.channels.iter().zip(channel_rates) {
                    if !valid(r) {
                        return bad(format!("invalid rate {r} of MT{}", reaction.mt));
                    }
                    entries.push((j, j, -r));
                    gain(&mut entries, j, &channel.products, r);
                }
            }
            if !valid(rate.fission) {
                return bad(format!("invalid fission rate {}", rate.fission));
            }
            match &nuclide.fission {
                None => {
                    if rate.fission != 0.0 || !rate.yield_weights.is_empty() {
                        return bad(
                            "fission rate or yield weights for a nuclide without fission".into(),
                        );
                    }
                }
                Some(fission) => {
                    let weights = &rate.yield_weights;
                    if weights.len() != fission.energies.len() || !weights.iter().all(|&w| valid(w))
                    {
                        return bad(format!(
                            "yield weights {weights:?} for {} energies",
                            fission.energies.len()
                        ));
                    }
                    let sum: f64 = weights.iter().sum();
                    if (sum - 1.0).abs() > SUM_TOLERANCE {
                        return bad(format!("yield weights sum to {sum}"));
                    }
                    entries.push((j, j, -rate.fission));
                    for (&w, table) in weights.iter().zip(&fission.yields) {
                        if w == 0.0 {
                            continue;
                        }
                        for y in table {
                            entries.push((y.nuclide, j, rate.fission * w * y.value));
                        }
                    }
                }
            }
        }
        Ok(entries)
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
        self.energies
            .windows(2)
            .map(|w| (w[0] * w[1]).sqrt())
            .collect()
    }
}

/// 構造体が入る位置の値がすべて map であることを確かめる（serde は構造体を配列からも読むため）
fn require_maps(value: &rmpv::Value) -> Result<(), MsgpackError> {
    use msgpack::{array, map};
    let products = |owner: &rmpv::Value| -> Result<(), MsgpackError> {
        for product in array(map(owner)?, "products") {
            map(product)?;
        }
        Ok(())
    };
    for nuclide in array(map(value)?, "nuclides") {
        let entries = map(nuclide)?;
        for (key, item) in entries {
            match key.as_str() {
                Some("decay") if !item.is_nil() => {
                    let decay = map(item)?;
                    for (k, v) in decay {
                        if k.as_str() == Some("energy") {
                            map(v)?;
                        }
                    }
                    for mode in array(decay, "modes") {
                        products(mode)?;
                    }
                }
                Some("reactions") => {
                    for reaction in item.as_array().map_or(&[][..], Vec::as_slice) {
                        for channel in array(map(reaction)?, "channels") {
                            products(channel)?;
                            for (k, v) in map(channel)? {
                                if k.as_str() == Some("final_state") && !v.is_nil() {
                                    map(v)?;
                                }
                            }
                        }
                    }
                }
                Some("fission") if !item.is_nil() => {
                    for table in array(map(item)?, "yields") {
                        for y in table.as_array().map_or(&[][..], Vec::as_slice) {
                            map(y)?;
                        }
                    }
                }
                _ => {}
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
