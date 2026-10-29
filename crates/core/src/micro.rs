//! ノードのミクロ燃焼（ADR-0216）。
//!
//! 追う核種（テーブルの `micro` の `nuclides`）と、補助の核種（追う核種へチェーンの崩壊で流れ込む核種のうち、追う核種の反応、
//! 核分裂の収率、補助の核種の崩壊で生まれうるもの）を、格子計算の燃焼と同じチェーンで解く。生まれえない先行核は数密度が 0 のまま
//! なので解かない。追う核種の反応率は、テーブルの微視的断面積 × ノードの少数群の中性子束で、中性子束の大きさはノードの出力密度で
//! 決める（追う核種の核分裂率 × チェーンの Q 値の和が出力密度になる。格子計算の燃焼と同じ規格化。ADR-0173）。補助の核種は崩壊
//! だけ。核分裂の収率の重みは、テーブルの収率の境界より下の核分裂の割合から ADR-0183 の点の選び方で作る。
//!
//! ノードの少数群定数には線形の補正 Σᵢ σᵢ (Nᵢ − N̄ᵢ⁰) を足す（全断面積、吸収、νΣf、Σf、散乱行列。ADR-0216）。

use std::collections::{BTreeMap, BTreeSet};

use xlc_depletion::DepletionMatrix;
use xlc_depletion::chain::Chain;
use xlc_format::fewgroup::FewGroup;
use xlc_format::table::MicroNuclide;

use crate::state::MicroAt;

/// 1 eV (J)
const ELECTRON_VOLT: f64 = 1.602_176_634e-19;
/// 1 barn (cm²)。数密度 1/(barn cm) を 1/cm³ にするのにも使う
const BARN: f64 = 1.0e-24;

/// 核分裂の収率の熱の点 (eV)。格子計算の燃焼（crates/lattice/src/depletion.rs の `YIELD_THERMAL`）と同じ値（ADR-0183）
const YIELD_THERMAL: f64 = 0.0253;
/// 核分裂の収率の高速の点 (eV)。格子計算の燃焼の `YIELD_FAST` と同じ値（ADR-0183）
const YIELD_FAST: f64 = 5.0e5;

/// 追う核種の 1 つの、燃焼行列の項の材料（生成物は解く核種の番号。解かない核種への生成物は落とす）
#[derive(Debug, Clone, PartialEq)]
struct TrackedTerms {
    /// 解く核種の中の番号
    column: usize,
    /// [チェーンの反応][経路] の生成物（番号、個数）
    reactions: Vec<Vec<Vec<(usize, f64)>>>,
    /// 核分裂: Q 値 (eV)、収率のエネルギー (eV)、エネルギーごとの生成物（番号、収率）
    fission: Option<FissionTerms>,
}

#[derive(Debug, Clone, PartialEq)]
struct FissionTerms {
    q: f64,
    energies: Vec<f64>,
    yields: Vec<Vec<(usize, f64)>>,
}

/// 追う核種と補助の核種の燃焼行列の材料（チェーンから 1 度だけ作る）
#[derive(Debug, Clone, PartialEq)]
pub struct MicroChain {
    /// 追う核種の名前（テーブルの `nuclides` の順）
    pub nuclides: Vec<String>,
    /// 解く核種（追う核種と補助の核種）のチェーンの名前（チェーンの順）
    pub solved: Vec<String>,
    /// 崩壊の項（行、列、値 (1/s)）。反応率によらない
    decay: Vec<(usize, usize, f64)>,
    /// 追う核種ごと（テーブルの順）
    tracked: Vec<TrackedTerms>,
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
    /// 燃焼行列の材料を作る。追う核種がチェーンにないか重なるか、反応の並び（MT と経路の数）がチェーンと違うか、テーブルが核分裂
    /// するのにチェーンが収率を持たなければエラー
    pub fn new(
        chain: &Chain,
        nuclides: &[String],
        reference: &[MicroNuclide],
    ) -> Result<Self, String> {
        if nuclides.len() != reference.len() {
            return Err(format!(
                "{} micro cross sections for {} nuclides",
                reference.len(),
                nuclides.len()
            ));
        }
        let index: BTreeMap<&str, usize> = chain
            .nuclides
            .iter()
            .enumerate()
            .map(|(i, n)| (n.name.as_str(), i))
            .collect();
        let mut tracked = Vec::with_capacity(nuclides.len());
        for (name, micro) in nuclides.iter().zip(reference) {
            let &i = index
                .get(name.as_str())
                .ok_or_else(|| format!("micro nuclide {name} is not in the chain"))?;
            if tracked.contains(&i) {
                return Err(format!("micro nuclide {name} appears twice"));
            }
            let chain_shape: Vec<(u32, usize)> = chain.nuclides[i]
                .reactions
                .iter()
                .map(|r| (r.mt, r.channels.len()))
                .collect();
            let table_shape: Vec<(u32, usize)> = micro
                .reactions
                .iter()
                .map(|r| (r.mt, r.channels.len()))
                .collect();
            if chain_shape != table_shape {
                return Err(format!(
                    "micro nuclide {name}: the table reactions {table_shape:?} differ from the chain reactions \
                     {chain_shape:?} (MT and channels. ADR-0172)"
                ));
            }
            if micro.thermal_fission.is_some() && chain.nuclides[i].fission.is_none() {
                return Err(format!(
                    "micro nuclide {name} fissions in the table but has no fission yields in the chain"
                ));
            }
            tracked.push(i);
        }
        let solved = solved_nuclides(chain, &tracked);
        let position: BTreeMap<usize, usize> =
            solved.iter().enumerate().map(|(k, &i)| (i, k)).collect();
        let at = |products: &mut dyn Iterator<Item = (usize, f64)>| -> Vec<(usize, f64)> {
            products
                .filter_map(|(i, v)| position.get(&i).map(|&k| (k, v)))
                .collect()
        };
        let mut decay = Vec::new();
        for (k, &i) in solved.iter().enumerate() {
            if let Some(d) = &chain.nuclides[i].decay {
                let lambda = std::f64::consts::LN_2 / d.half_life;
                decay.push((k, k, -lambda));
                for mode in &d.modes {
                    for (row, count) in at(&mut mode
                        .products
                        .iter()
                        .map(|p| (p.nuclide, f64::from(p.count))))
                    {
                        decay.push((row, k, mode.ratio * lambda * count));
                    }
                }
            }
        }
        let tracked = tracked
            .iter()
            .map(|&i| {
                let nuclide = &chain.nuclides[i];
                TrackedTerms {
                    column: position[&i],
                    reactions: nuclide
                        .reactions
                        .iter()
                        .map(|r| {
                            r.channels
                                .iter()
                                .map(|c| {
                                    at(&mut c
                                        .products
                                        .iter()
                                        .map(|p| (p.nuclide, f64::from(p.count))))
                                })
                                .collect()
                        })
                        .collect(),
                    fission: nuclide.fission.as_ref().map(|f| FissionTerms {
                        q: f.q,
                        energies: f.energies.clone(),
                        yields: f
                            .yields
                            .iter()
                            .map(|table| at(&mut table.iter().map(|y| (y.nuclide, y.value))))
                            .collect(),
                    }),
                }
            })
            .collect();
        Ok(Self {
            nuclides: nuclides.to_vec(),
            solved: solved
                .iter()
                .map(|&i| chain.nuclides[i].name.clone())
                .collect(),
            decay,
            tracked,
        })
    }

    /// 補助の核種の数
    pub fn auxiliary(&self) -> usize {
        self.solved.len() - self.tracked.len()
    }

    /// 追う核種の数密度 `tracked`（テーブルの順）から、解く核種の数密度（補助の核種は 0）
    pub fn initial(&self, tracked: &[f64]) -> Vec<f64> {
        let mut n = vec![0.0; self.solved.len()];
        for (t, &v) in self.tracked.iter().zip(tracked) {
            n[t.column] = v;
        }
        n
    }

    /// 解く核種の数密度 `n` のうち、追う核種のもの（テーブルの順）
    pub fn tracked(&self, n: &[f64]) -> Vec<f64> {
        self.tracked.iter().map(|t| n[t.column]).collect()
    }

    /// ノードの微視的断面積 `micro`、少数群の中性子束の形 `flux`、解く核種の数密度 `n` (1/(barn cm)) で、追う核種の核分裂の出力が
    /// 出力密度 `power_density` (W/cm³) になるように規格化した、追う核種ごとの反応率。負の数密度は出力に数えない（格子計算の燃焼と
    /// 同じ。ADR-0173）
    pub fn rates(
        &self,
        micro: &MicroAt,
        flux: &[f64],
        n: &[f64],
        power_density: f64,
    ) -> Result<Vec<TrackedRates>, String> {
        if micro.nuclides.len() != self.tracked.len() || n.len() != self.solved.len() {
            return Err(format!(
                "{} micro cross sections and {} densities for {} tracked and {} solved nuclides",
                micro.nuclides.len(),
                n.len(),
                self.tracked.len(),
                self.solved.len()
            ));
        }
        let collapse = |sigma: &[f64]| -> f64 { sigma.iter().zip(flux).map(|(s, f)| s * f).sum() };
        // 中性子束の形のままの 1 原子あたりの率（barn × 中性子束）
        let mut power = 0.0;
        let mut rates = Vec::with_capacity(self.tracked.len());
        for (t, sigma) in self.tracked.iter().zip(&micro.nuclides) {
            let fission = collapse(&sigma.fission);
            if let Some(f) = &t.fission {
                power += n[t.column].max(0.0) / BARN * fission * BARN * f.q * ELECTRON_VOLT;
            }
            rates.push(TrackedRates {
                reactions: sigma
                    .reactions
                    .iter()
                    .map(|r| r.channels.iter().map(|c| collapse(c) * BARN).collect())
                    .collect(),
                fission: if t.fission.is_some() {
                    fission * BARN
                } else {
                    0.0
                },
                yield_weights: match (&t.fission, sigma.thermal_fission) {
                    (Some(f), Some(thermal)) => yield_weights(&f.energies, thermal),
                    // 核分裂しない核種の収率の重みは使わない。和が 1 になるように最初の点に置く（`Chain::zero_rates` と同じ）
                    (Some(f), None) => {
                        let mut w = vec![0.0; f.energies.len()];
                        w[0] = 1.0;
                        w
                    }
                    (None, _) => Vec::new(),
                },
            });
        }
        if !(power > 0.0 && power.is_finite()) {
            return Err(format!(
                "the tracked nuclides of the node have no fission power ({power}); the micro depletion \
                 normalizes the flux by them"
            ));
        }
        let scale = power_density / power;
        for r in &mut rates {
            r.fission *= scale;
            for c in r.reactions.iter_mut().flatten() {
                *c *= scale;
            }
        }
        Ok(rates)
    }

    /// 追う核種の反応率 `rates` の、解く核種の燃焼行列（`Chain::matrix_entries` と同じ項。解かない核種への生成物は落とし、
    /// 消滅は対角に残す）
    pub fn matrix(&self, rates: &[TrackedRates]) -> Result<DepletionMatrix, String> {
        if rates.len() != self.tracked.len() {
            return Err(format!(
                "{} rates for {} tracked nuclides",
                rates.len(),
                self.tracked.len()
            ));
        }
        let mut entries = self.decay.clone();
        for (t, r) in self.tracked.iter().zip(rates) {
            let j = t.column;
            for (channels, rate) in t.reactions.iter().zip(&r.reactions) {
                for (products, &rate) in channels.iter().zip(rate) {
                    if rate > 0.0 {
                        entries.push((j, j, -rate));
                        entries.extend(products.iter().map(|&(row, count)| (row, j, count * rate)));
                    }
                }
            }
            if let Some(f) = &t.fission
                && r.fission > 0.0
            {
                entries.push((j, j, -r.fission));
                for (&w, products) in r.yield_weights.iter().zip(&f.yields) {
                    if w != 0.0 {
                        entries.extend(
                            products
                                .iter()
                                .map(|&(row, value)| (row, j, r.fission * w * value)),
                        );
                    }
                }
            }
        }
        DepletionMatrix::from_triplets(self.solved.len(), &entries).map_err(|e| e.to_string())
    }
}

/// 解く核種（チェーンの番号の昇順）: 追う核種 `tracked` と、追う核種へ崩壊で流れ込む核種（崩壊の逆向きのクロージャ）のうち、
/// 追う核種の反応と収率、そこから崩壊でたどれるもの（補助の核種は崩壊だけなので、ほかの経路では生まれない）
fn solved_nuclides(chain: &Chain, tracked: &[usize]) -> Vec<usize> {
    let n = chain.nuclides.len();
    // 崩壊の逆向き: 娘 → 親
    let mut parents: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (i, nuclide) in chain.nuclides.iter().enumerate() {
        for mode in nuclide.decay.iter().flat_map(|d| &d.modes) {
            for p in &mode.products {
                parents[p.nuclide].push(i);
            }
        }
    }
    let mut feeding: BTreeSet<usize> = tracked.iter().copied().collect();
    let mut stack: Vec<usize> = tracked.to_vec();
    while let Some(i) = stack.pop() {
        for &p in &parents[i] {
            if feeding.insert(p) {
                stack.push(p);
            }
        }
    }
    // 追う核種から、反応、収率、崩壊でたどる（崩壊の先が追う核種へ流れ込む核種の中にあるものだけ）
    let mut reached: BTreeSet<usize> = tracked.iter().copied().collect();
    let mut stack: Vec<usize> = tracked.to_vec();
    let is_tracked: BTreeSet<usize> = tracked.iter().copied().collect();
    while let Some(i) = stack.pop() {
        let nuclide = &chain.nuclides[i];
        let mut next: Vec<usize> = nuclide
            .decay
            .iter()
            .flat_map(|d| &d.modes)
            .flat_map(|m| &m.products)
            .map(|p| p.nuclide)
            .collect();
        if is_tracked.contains(&i) {
            next.extend(
                nuclide
                    .reactions
                    .iter()
                    .flat_map(|r| &r.channels)
                    .flat_map(|c| &c.products)
                    .map(|p| p.nuclide),
            );
            next.extend(
                nuclide
                    .fission
                    .iter()
                    .flat_map(|f| &f.yields)
                    .flatten()
                    .map(|y| y.nuclide),
            );
        }
        for k in next {
            if feeding.contains(&k) && reached.insert(k) {
                stack.push(k);
            }
        }
    }
    reached.into_iter().collect()
}

/// 収率の境界より下の核分裂の割合 `thermal` の、収率のエネルギーの点 `energies` ごとの重み（ADR-0183）。点の選び方は格子計算の燃焼
/// （crates/lattice/src/depletion.rs の `yield_weights`。OpenMC 0.15.2 の `FissionYieldCutoffHelper`）と同じ: 熱と高速の点がちょうど
/// あればそれを使い、どちらかがなければ、境界が最初の点以下なら最初の点、最後の点以上なら最後の点だけ、そうでなければ境界より下で
/// 熱に、境界以上で高速に最も近い点
pub fn yield_weights(energies: &[f64], thermal: f64) -> Vec<f64> {
    // 熱と高速の点の幾何平均（112.47 eV）。格子計算の燃焼の `yield_cutoff` と同じ
    let cutoff = (YIELD_THERMAL * YIELD_FAST).sqrt();
    let mut weights = vec![0.0; energies.len()];
    let exact = |e: f64| energies.iter().position(|&p| p == e);
    let (t, f) = match (exact(YIELD_THERMAL), exact(YIELD_FAST)) {
        (Some(t), Some(f)) => (t, f),
        (t, f) => {
            if cutoff <= energies[0] {
                weights[0] = 1.0;
                return weights;
            }
            let last = energies.len() - 1;
            if cutoff >= energies[last] {
                weights[last] = 1.0;
                return weights;
            }
            let split = energies.partition_point(|&e| e < cutoff);
            // 同じ近さなら前の点（Python の min と同じ）
            let nearest = |from: usize, to: usize, target: f64| {
                (from..to).fold(from, |best, k| {
                    if (energies[k] - target).abs() < (energies[best] - target).abs() {
                        k
                    } else {
                        best
                    }
                })
            };
            (
                t.unwrap_or_else(|| nearest(0, split, YIELD_THERMAL)),
                f.unwrap_or_else(|| nearest(split, energies.len(), YIELD_FAST)),
            )
        }
    };
    weights[t] += thermal;
    weights[f] += 1.0 - thermal;
    weights
}

/// 少数群定数 `f` に、線形の補正 Σᵢ σᵢ (Nᵢ − N̄ᵢ⁰) を足す（`n` は追う核種の数密度、テーブルの順）。吸収の σ は全断面積から散乱の
/// 行の和を引いたもの。χ、拡散係数、ADF、形状関数と k∞ は補正しない（ADR-0216）
pub fn correct(f: &mut FewGroup, micro: &MicroAt, n: &[f64]) {
    for ((sigma, &nominal), &density) in micro.nuclides.iter().zip(&micro.densities).zip(n) {
        let d = density - nominal;
        if d == 0.0 {
            continue;
        }
        for (g, row) in sigma.scatter.iter().enumerate() {
            f.total[g] += d * sigma.total[g];
            f.absorption[g] += d * (sigma.total[g] - row.iter().sum::<f64>());
            f.nu_fission[g] += d * sigma.nu_fission[g];
            f.fission[g] += d * sigma.fission[g];
            for (a, s) in f.scatter[g].iter_mut().zip(row) {
                *a += d * s;
            }
        }
    }
}

#[cfg(test)]
mod tests;
