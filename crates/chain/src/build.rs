//! ADR-0164 の規則で、読んだ評価から全核種の燃焼チェーンを組み立てる。

use std::collections::{BTreeMap, HashMap};

use anyhow::{Result, anyhow, bail};
use xlc_depletion::chain::{
    Chain, Channel, Decay, DecayEnergy, DecayMode as ChainMode, FinalState, Fission, Nuclide,
    Product, Reaction, Yield,
};

use crate::elements::nuclide_name;
use crate::endf::{
    DecayEvaluation, FinalStateEntry, Material, NeutronEvaluation, RawYield, Step, YieldEvaluation,
};

/// チェーンの反応の MT（代表の MT。MF8 の生成物の状態もこの MT で読む）
pub const REACTION_MTS: [u32; 6] = [16, 17, 37, 102, 103, 107];

/// 組み立てたチェーンと、データにないものを規則で補った箇所の警告。
#[derive(Debug, Clone, PartialEq)]
pub struct Built {
    pub chain: Chain,
    pub warnings: Vec<String>,
}

/// チェーンの反応: 代表の MT、MF3 でその反応とみなす MT、(ΔZ, ΔA)、生まれる軽い核種 (Z, A)
struct ReactionSpec {
    mt: u32,
    /// MT の範囲（両端を含む）
    mts: &'static [(u32, u32)],
    dz: i64,
    da: i64,
    light: Option<(u32, u32)>,
}

/// OpenMC の `Chain.from_endf` の既定と同じ反応（ADR-0164）。部分の MT（875〜891 など）だけがある評価も含める
const REACTIONS: [ReactionSpec; 6] = [
    ReactionSpec {
        mt: 16,
        mts: &[(16, 16), (875, 891)],
        dz: 0,
        da: -1,
        light: None,
    },
    ReactionSpec {
        mt: 17,
        mts: &[(17, 17)],
        dz: 0,
        da: -2,
        light: None,
    },
    ReactionSpec {
        mt: 37,
        mts: &[(37, 37)],
        dz: 0,
        da: -3,
        light: None,
    },
    ReactionSpec {
        mt: 102,
        mts: &[(102, 102)],
        dz: 0,
        da: 1,
        light: None,
    },
    ReactionSpec {
        mt: 103,
        mts: &[(103, 103), (600, 649)],
        dz: -1,
        da: 0,
        light: Some((1, 1)),
    },
    ReactionSpec {
        mt: 107,
        mts: &[(107, 107), (800, 849)],
        dz: -2,
        da: -3,
        light: Some((2, 4)),
    },
];

/// 核分裂の MT（全、第 1〜3 チャンス、第 4 チャンス）
const FISSION_MTS: [u32; 5] = [18, 19, 20, 21, 38];

/// 崩壊の分岐比の和の、1 からのずれの許容。ENDF/B-VIII.0 の decay の最大は 9e-7（Ba-114）で、丸めの範囲。これを
/// 超えるずれは、型の抜けなどの評価の誤りとして止める
const DECAY_SUM_TOLERANCE: f64 = 1e-5;

/// MF8 の生成物の状態の励起エネルギーを、decay の状態の励起エネルギーに合わせるときの相対の許容。同じ評価の値を
/// コピーしたものは一致し、別の準位とは 1e-3 より大きく離れる
const EXCITATION_TOLERANCE: f64 = 1e-3;

type Key = (u32, u32, u32);

fn key(m: &Material) -> Key {
    (m.z, m.a, m.state)
}

fn name(key: Key) -> Result<String> {
    nuclide_name(key.0, key.1, key.2).ok_or_else(|| anyhow!("no element symbol for Z = {}", key.0))
}

fn step_name(step: Step) -> &'static str {
    match step {
        Step::BetaMinus => "beta-",
        Step::EcBetaPlus => "ec/beta+",
        Step::IsomericTransition => "IT",
        Step::Alpha => "alpha",
        Step::Neutron => "n",
        Step::SpontaneousFission => "sf",
        Step::Proton => "p",
    }
}

/// 核種の添字の引き当てと、軽い核種
struct Index {
    keys: HashMap<Key, usize>,
    h1: usize,
    he4: usize,
}

impl Index {
    fn get(&self, key: Key) -> Option<usize> {
        self.keys.get(&key).copied()
    }
}

/// (Z, A) に変化を足す。負になれば `None`
fn shifted(z: u32, a: u32, dz: i64, da: i64) -> Option<(u32, u32)> {
    let z = u32::try_from(i64::from(z) + dz).ok()?;
    let a = u32::try_from(i64::from(a) + da).ok()?;
    Some((z, a))
}

/// 崩壊の型の生成物。自発核分裂を含む型は空。娘がチェーンになければ娘を除き（中性子そのものは追わない）、警告を足す
fn decay_products(
    parent: &str,
    from: Key,
    mode: &crate::endf::DecayMode,
    index: &Index,
    warnings: &mut Vec<String>,
) -> Result<Vec<Product>> {
    if mode.steps.contains(&Step::SpontaneousFission) {
        return Ok(Vec::new());
    }
    let (mut z, mut a) = (i64::from(from.0), i64::from(from.1));
    let (mut alphas, mut protons) = (0, 0);
    for step in &mode.steps {
        match step {
            Step::BetaMinus => z += 1,
            Step::EcBetaPlus => z -= 1,
            Step::IsomericTransition | Step::SpontaneousFission => {}
            Step::Alpha => {
                z -= 2;
                a -= 4;
                alphas += 1;
            }
            Step::Neutron => a -= 1,
            Step::Proton => {
                z -= 1;
                a -= 1;
                protons += 1;
            }
        }
    }
    let mut products = Vec::new();
    match (u32::try_from(z), u32::try_from(a)) {
        (Ok(0), _) => {}
        (Ok(z), Ok(a)) => {
            let daughter = (z, a, mode.final_state);
            match index.get(daughter) {
                Some(i) => products.push(Product { nuclide: i, count: 1 }),
                None => warnings.push(format!(
                    "{parent}: the daughter {} of {} is not in the decay data; the branch produces no daughter",
                    name(daughter)?,
                    mode_name(&mode.steps)
                )),
            }
        }
        _ => bail!(
            "{parent}: decay {} leaves Z = {z}, A = {a}",
            mode_name(&mode.steps)
        ),
    }
    if protons > 0 {
        products.push(Product {
            nuclide: index.h1,
            count: protons,
        });
    }
    if alphas > 0 {
        products.push(Product {
            nuclide: index.he4,
            count: alphas,
        });
    }
    Ok(products)
}

fn mode_name(steps: &[Step]) -> String {
    steps
        .iter()
        .map(|&s| step_name(s))
        .collect::<Vec<_>>()
        .join(",")
}

fn decay(
    parent: &str,
    evaluation: &DecayEvaluation,
    index: &Index,
    warnings: &mut Vec<String>,
) -> Result<Option<Decay>> {
    let data = &evaluation.decay;
    if data.stable {
        return Ok(None);
    }
    if data.half_life == 0.0 {
        warnings.push(format!(
            "{parent}: unstable without a half-life (modes: {}); treated as stable",
            data.modes
                .iter()
                .map(|m| mode_name(&m.steps))
                .collect::<Vec<_>>()
                .join("; ")
        ));
        return Ok(None);
    }
    let sum: f64 = data.modes.iter().map(|m| m.ratio).sum();
    if !sum.is_finite() || (sum - 1.0).abs() > DECAY_SUM_TOLERANCE {
        bail!("{parent}: decay ratios sum to {sum}");
    }
    let mut modes = Vec::new();
    for mode in &data.modes {
        if mode.ratio < 0.0 {
            bail!("{parent}: negative decay ratio {}", mode.ratio);
        }
        if mode.ratio == 0.0 {
            continue;
        }
        modes.push(ChainMode {
            mode: mode_name(&mode.steps),
            ratio: mode.ratio / sum,
            products: decay_products(parent, key(&evaluation.material), mode, index, warnings)?,
        });
    }
    Ok(Some(Decay {
        half_life: data.half_life,
        energy: DecayEnergy {
            light: data.energy[0],
            electromagnetic: data.energy[1],
            heavy: data.energy[2],
        },
        modes,
    }))
}

/// 反応の生成物の状態を、decay の同じ (Z, A) の状態に励起エネルギーで合わせる。decay にない準位は、即発の γ で
/// 基底状態に落ちるので基底状態にし、`true` を返す（基底状態もなければ `None`）
fn final_state(
    parent: &str,
    mt: u32,
    entry: &FinalStateEntry,
    expected: (u32, u32),
    evaluations: &[&DecayEvaluation],
    index: &Index,
) -> Result<(Option<usize>, bool)> {
    if entry.za != 1000 * expected.0 + expected.1 {
        bail!(
            "{parent}: MF8/MT{mt} names product ZA {}, expected {}",
            entry.za,
            1000 * expected.0 + expected.1
        );
    }
    let matched = evaluations.iter().find(|e| {
        let m = &e.material;
        (m.z, m.a) == expected
            && (m.excitation - entry.excitation).abs() <= EXCITATION_TOLERANCE * entry.excitation
    });
    Ok(match matched {
        Some(e) => (index.get(key(&e.material)), false),
        None => (index.get((expected.0, expected.1, 0)), true),
    })
}

fn reactions(
    parent: &str,
    evaluation: &NeutronEvaluation,
    evaluations: &[&DecayEvaluation],
    index: &Index,
    warnings: &mut Vec<String>,
) -> Result<Vec<Reaction>> {
    let m = &evaluation.material;
    let mut reactions = Vec::new();
    for spec in &REACTIONS {
        if !evaluation
            .mts
            .iter()
            .any(|mt| spec.mts.iter().any(|&(lo, hi)| (lo..=hi).contains(mt)))
        {
            continue;
        }
        let (z, a) = shifted(m.z, m.a, spec.dz, spec.da)
            .ok_or_else(|| anyhow!("{parent}: MT{} leaves no nucleus", spec.mt))?;
        let light: Vec<Product> = spec
            .light
            .map(|(lz, la)| {
                let nuclide = if (lz, la) == (1, 1) {
                    index.h1
                } else {
                    index.he4
                };
                Product { nuclide, count: 1 }
            })
            .into_iter()
            .collect();
        let states = evaluation
            .final_states
            .iter()
            .find(|(mt, _)| *mt == spec.mt)
            .map(|(_, entries)| entries);
        let channels = match states {
            Some(entries) => {
                let mut grounded = Vec::new();
                let mut channels = Vec::new();
                for entry in entries {
                    let (product, to_ground) =
                        final_state(parent, spec.mt, entry, (z, a), evaluations, index)?;
                    if to_ground {
                        grounded.push(format!("LFS {} at {} eV", entry.lfs, entry.excitation));
                    }
                    let mut products: Vec<Product> = product
                        .map(|nuclide| Product { nuclide, count: 1 })
                        .into_iter()
                        .collect();
                    if product.is_none() {
                        warnings.push(format!(
                            "{parent}: the product {} of MT{} is not in the decay data; the reaction produces no residual",
                            name((z, a, 0))?,
                            spec.mt
                        ));
                    }
                    products.extend(light.iter().copied());
                    channels.push(Channel {
                        final_state: Some(FinalState {
                            mf: entry.mf,
                            lfs: entry.lfs,
                        }),
                        products,
                    });
                }
                if !grounded.is_empty() {
                    warnings.push(format!(
                        "{parent}: MF8/MT{} final states not in the decay data go to the ground state of {}: {}",
                        spec.mt,
                        name((z, a, 0))?,
                        grounded.join(", ")
                    ));
                }
                channels
            }
            None => {
                let mut products = Vec::new();
                match index.get((z, a, 0)) {
                    Some(i) => products.push(Product { nuclide: i, count: 1 }),
                    None => warnings.push(format!(
                        "{parent}: the product {} of MT{} is not in the decay data; the reaction produces no residual",
                        name((z, a, 0))?,
                        spec.mt
                    )),
                }
                products.extend(light.iter().copied());
                vec![Channel {
                    final_state: None,
                    products,
                }]
            }
        };
        reactions.push(Reaction {
            mt: spec.mt,
            channels,
        });
    }
    Ok(reactions)
}

/// 収率の生成物をチェーンの核種にする。準安定状態がなければ基底状態、基底状態もなければ β− の向き（ADR-0164）
fn yield_product(za: u32, state: u32, index: &Index) -> Result<(usize, bool)> {
    let (mut z, a) = (za / 1000, za % 1000);
    if let Some(i) = index.get((z, a, state)) {
        return Ok((i, false));
    }
    loop {
        if let Some(i) = index.get((z, a, 0)) {
            return Ok((i, true));
        }
        z += 1;
        if z > 118 {
            bail!("no nuclide of mass {a} in the decay data for the fission product ZA {za}");
        }
    }
}

/// nfy の評価を、チェーンの収率にする
fn fission_yields(
    source: &str,
    evaluation: &YieldEvaluation,
    index: &Index,
    nuclides: &[Nuclide],
    warnings: &mut Vec<String>,
) -> Result<Fission> {
    let mut tables: Vec<&(f64, Vec<RawYield>)> = evaluation.tables.iter().collect();
    tables.sort_by(|x, y| x.0.total_cmp(&y.0));
    if tables.windows(2).any(|w| w[0].0 == w[1].0) {
        bail!("{source}: repeated yield energies");
    }
    let mut moved: Vec<String> = Vec::new();
    let mut energies = Vec::new();
    let mut yields = Vec::new();
    for (energy, raw) in tables {
        let mut by_nuclide: BTreeMap<usize, f64> = BTreeMap::new();
        for y in raw {
            if !(y.value.is_finite() && y.value >= 0.0) {
                bail!("{source}: invalid yield {} of ZA {}", y.value, y.za);
            }
            if y.value == 0.0 {
                continue;
            }
            let (i, replaced) = yield_product(y.za, y.state, index)?;
            if replaced {
                let label = format!(
                    "{} -> {}",
                    name((y.za / 1000, y.za % 1000, y.state))?,
                    nuclides[i].name
                );
                if !moved.contains(&label) {
                    moved.push(label);
                }
            }
            *by_nuclide.entry(i).or_insert(0.0) += y.value;
        }
        energies.push(*energy);
        yields.push(
            by_nuclide
                .into_iter()
                .map(|(nuclide, value)| Yield { nuclide, value })
                .collect(),
        );
    }
    if !moved.is_empty() {
        warnings.push(format!(
            "{source}: yields to {} states outside the decay data moved: {}",
            moved.len(),
            moved.join(", ")
        ));
    }
    Ok(Fission {
        yield_source: source.to_string(),
        energies,
        yields,
    })
}

/// 収率のない核分裂する核種に使う収率の評価: 原子番号の最も近い元素（同じなら軽い方）の、質量数の最も近い核種（同じなら
/// 軽い方、同じ質量数なら基底状態から）（ADR-0164）
fn nearest_yields(target: Key, yields: &[&YieldEvaluation]) -> Option<usize> {
    let distance = |m: &Material| {
        (
            m.z.abs_diff(target.0),
            m.z,
            m.a.abs_diff(target.1),
            m.a,
            m.state,
        )
    };
    (0..yields.len()).min_by_key(|&i| distance(&yields[i].material))
}

/// 評価からチェーンを組み立てる。評価の順は問わない。
pub fn build(
    decays: &[DecayEvaluation],
    neutrons: &[NeutronEvaluation],
    yields: &[YieldEvaluation],
) -> Result<Built> {
    let mut warnings = Vec::new();
    // 中性子そのもの（Z = 0）はチェーンに入れない
    let mut evaluations: Vec<&DecayEvaluation> =
        decays.iter().filter(|d| d.material.z > 0).collect();
    evaluations.sort_by_key(|d| key(&d.material));
    if let Some(w) = evaluations
        .windows(2)
        .find(|w| key(&w[0].material) == key(&w[1].material))
    {
        bail!(
            "{} appears twice in the decay data",
            name(key(&w[0].material))?
        );
    }
    let keys: HashMap<Key, usize> = evaluations
        .iter()
        .enumerate()
        .map(|(i, e)| (key(&e.material), i))
        .collect();
    let light = |k: Key| {
        keys.get(&k)
            .copied()
            .ok_or_else(|| anyhow!("{} is not in the decay data", name(k).unwrap_or_default()))
    };
    let index = Index {
        h1: light((1, 1, 0))?,
        he4: light((2, 4, 0))?,
        keys: keys.clone(),
    };

    let mut nuclides = Vec::with_capacity(evaluations.len());
    for evaluation in &evaluations {
        let k = key(&evaluation.material);
        let parent = name(k)?;
        let decay = decay(&parent, evaluation, &index, &mut warnings)?;
        nuclides.push(Nuclide {
            name: parent,
            z: k.0,
            a: k.1,
            state: k.2,
            excitation: evaluation.material.excitation,
            decay,
            reactions: Vec::new(),
            fission: None,
        });
    }

    let mut yield_evaluations: Vec<&YieldEvaluation> = yields.iter().collect();
    yield_evaluations.sort_by_key(|y| key(&y.material));
    for y in &yield_evaluations {
        if index.get(key(&y.material)).is_none() {
            bail!(
                "fission yields for {}, which is not in the decay data",
                name(key(&y.material))?
            );
        }
    }
    let mut fissionable = Vec::new();
    // 中性子そのものの評価（n + n）もチェーンに入れない
    for evaluation in neutrons.iter().filter(|n| n.material.z > 0) {
        let k = key(&evaluation.material);
        let i = index.get(k).ok_or_else(|| {
            anyhow!(
                "neutron data for {}, which is not in the decay data",
                name(k).unwrap_or_default()
            )
        })?;
        if !nuclides[i].reactions.is_empty() {
            bail!("{} appears twice in the neutron data", nuclides[i].name);
        }
        let parent = nuclides[i].name.clone();
        nuclides[i].reactions =
            reactions(&parent, evaluation, &evaluations, &index, &mut warnings)?;
        if evaluation.mts.iter().any(|mt| FISSION_MTS.contains(mt)) {
            fissionable.push(i);
        }
    }
    fissionable.sort_unstable();
    for i in fissionable {
        let k = (nuclides[i].z, nuclides[i].a, nuclides[i].state);
        let source = nearest_yields(k, &yield_evaluations)
            .ok_or_else(|| anyhow!("no fission yields for {}", nuclides[i].name))?;
        let evaluation = yield_evaluations[source];
        let source_name = name(key(&evaluation.material))?;
        if key(&evaluation.material) != k {
            warnings.push(format!(
                "{}: no fission yields; using those of {source_name}",
                nuclides[i].name
            ));
        }
        let fission = fission_yields(&source_name, evaluation, &index, &nuclides, &mut warnings)?;
        nuclides[i].fission = Some(fission);
    }

    let chain = Chain { nuclides };
    chain.check()?;
    Ok(Built { chain, warnings })
}

#[cfg(test)]
mod tests;
