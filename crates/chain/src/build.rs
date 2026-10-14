//! ADR-0164 の規則で、読んだ評価から全核種の燃焼チェーンを組み立てる。

use std::collections::{BTreeMap, HashMap};

use anyhow::{Result, anyhow, bail};
use xlc_depletion::chain::{
    Chain, Channel, Decay, DecayEnergy, DecayMode as ChainMode, FinalState, Fission, Nuclide,
    Product, Reaction, Yield, nuclide_name,
};
use xlc_format::library::{CHAIN_REACTIONS, ChainReaction};

use crate::endf::{
    DecayEvaluation, FinalStateEntry, Material, NeutronEvaluation, RawYield, Step, YieldEvaluation,
};

/// 組み立てたチェーンと、データにないものを規則で補った箇所の警告。
#[derive(Debug, Clone, PartialEq)]
pub struct Built {
    pub chain: Chain,
    pub warnings: Vec<String>,
}

/// チェーンの反応: 代表の MT と MF3 でその反応とみなす MT（xs と共通。ADR-0172）、(ΔZ, ΔA)、生まれる軽い核種 (Z, A)
struct ReactionSpec {
    reaction: ChainReaction,
    dz: i64,
    da: i64,
    light: Option<(u32, u32)>,
}

/// OpenMC の `Chain.from_endf` の既定と同じ反応（ADR-0164）。部分の MT（875〜891 など）だけがある評価も含める。MT の範囲は xs と
/// 共通（`CHAIN_REACTIONS`。ADR-0172）
const REACTIONS: [ReactionSpec; 6] = [
    ReactionSpec {
        reaction: CHAIN_REACTIONS[0],
        dz: 0,
        da: -1,
        light: None,
    },
    ReactionSpec {
        reaction: CHAIN_REACTIONS[1],
        dz: 0,
        da: -2,
        light: None,
    },
    ReactionSpec {
        reaction: CHAIN_REACTIONS[2],
        dz: 0,
        da: -3,
        light: None,
    },
    ReactionSpec {
        reaction: CHAIN_REACTIONS[3],
        dz: 0,
        da: 1,
        light: None,
    },
    ReactionSpec {
        reaction: CHAIN_REACTIONS[4],
        dz: -1,
        da: 0,
        light: Some((1, 1)),
    },
    ReactionSpec {
        reaction: CHAIN_REACTIONS[5],
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

/// MF8 の生成物の状態と、対応させた decay の準安定状態の励起エネルギーの相対の差が、これを超えれば警告に出す。同じ
/// 評価の値をコピーしたものは 1e-3 の内で一致する
const EXCITATION_TOLERANCE: f64 = 1e-3;

type Key = (u32, u32, u32);

fn key(m: &Material) -> Key {
    (m.z, m.a, m.state)
}

/// 核種の名前。名前のない (Z, A)（A < Z など）は警告とエラーの文に使えるよう、数で書く（チェーンの核種なら `Chain::check` が止める）
fn name(key: Key) -> Result<String> {
    Ok(nuclide_name(key.0, key.1, key.2)
        .unwrap_or_else(|| format!("(Z = {}, A = {}, state {})", key.0, key.1, key.2)))
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
            // decay にない準安定状態は基底状態へ（収率と MF8 の状態と同じ規則。ADR-0164）
            let ground = (z, a, 0);
            match (index.get(daughter), index.get(ground)) {
                (Some(i), _) => products.push(Product { nuclide: i, count: 1 }),
                (None, Some(i)) if mode.final_state > 0 => {
                    warnings.push(format!(
                        "{parent}: the daughter {} of {} is not in the decay data; sent to {}",
                        name(daughter)?,
                        mode_name(&mode.steps),
                        name(ground)?
                    ));
                    products.push(Product { nuclide: i, count: 1 });
                }
                _ => warnings.push(format!(
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
    Ok(merged(products))
}

/// 同じ核種の生成物を 1 つにまとめ、個数を足す（娘が放出した軽い核種と同じとき。Li-8 の β−α の He-4）
fn merged(products: Vec<Product>) -> Vec<Product> {
    let mut out: Vec<Product> = Vec::with_capacity(products.len());
    for p in products {
        match out.iter_mut().find(|q| q.nuclide == p.nuclide) {
            Some(q) => q.count += p.count,
            None => out.push(p),
        }
    }
    out
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
        if !data.modes.is_empty() {
            bail!("{parent}: a stable evaluation lists decay modes");
        }
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

/// MF8 の生成物の状態（LMF = 9、10）を、decay の状態に対応させる（ADR-0164）。ELFS = 0 は基底状態。励起状態は、評価が
/// 生成を与えた準安定状態を指すが、励起エネルギーが decay と食い違う評価がある（Am-243 (n,γ) の Am-244m が 72.9 keV と
/// 86.1 keV）ので、エネルギーの一致では引かず、decay の準安定状態に、相対の差の小さい組から 1 対 1 で対応させる。対応の
/// ない状態（準安定状態がない、または数が足りない）は、即発の γ で落ちる準位として基底状態にする。返すのは状態ごとの
/// 生成物の (Z, A, state) と、警告に出す対応の注記
fn assign_states(
    parent: &str,
    mt: u32,
    entries: &[FinalStateEntry],
    product: (u32, u32),
    evaluations: &[&DecayEvaluation],
) -> Result<(Vec<Key>, Vec<String>)> {
    let za = 1000 * product.0 + product.1;
    if let Some(e) = entries.iter().find(|e| e.za != za) {
        bail!(
            "{parent}: MF8/MT{mt} names product ZA {}, expected {za}",
            e.za
        );
    }
    let isomers: Vec<&Material> = evaluations
        .iter()
        .map(|e| &e.material)
        .filter(|m| (m.z, m.a) == product && m.state > 0)
        .collect();
    let relative = |e: &FinalStateEntry, m: &Material| {
        (e.excitation - m.excitation).abs() / e.excitation.max(m.excitation)
    };
    let mut pairs: Vec<(f64, usize, usize)> = Vec::new();
    for (i, e) in entries
        .iter()
        .enumerate()
        .filter(|(_, e)| e.excitation > 0.0)
    {
        for (j, m) in isomers.iter().enumerate() {
            pairs.push((relative(e, m), i, j));
        }
    }
    pairs.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)).then(x.2.cmp(&y.2)));
    let mut assigned: Vec<Option<usize>> = vec![None; entries.len()];
    let mut used = vec![false; isomers.len()];
    for (_, i, j) in pairs {
        if assigned[i].is_none() && !used[j] {
            assigned[i] = Some(j);
            used[j] = true;
        }
    }
    let mut keys = Vec::with_capacity(entries.len());
    let mut notes = Vec::new();
    for (e, a) in entries.iter().zip(assigned) {
        let ground = (product.0, product.1, 0);
        match a {
            Some(j) => {
                let m = isomers[j];
                let k = key(m);
                if relative(e, m) > EXCITATION_TOLERANCE {
                    notes.push(format!(
                        "LFS {} at {} eV to {} at {} eV",
                        e.lfs,
                        e.excitation,
                        name(k)?,
                        m.excitation
                    ));
                }
                keys.push(k);
            }
            None => {
                if e.excitation > 0.0 {
                    notes.push(format!(
                        "LFS {} at {} eV to the ground state (no isomer left)",
                        e.lfs, e.excitation
                    ));
                }
                keys.push(ground);
            }
        }
    }
    Ok((keys, notes))
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
        if !evaluation.mts.iter().any(|&mt| spec.reaction.includes(mt)) {
            continue;
        }
        let (z, a) = shifted(m.z, m.a, spec.dz, spec.da)
            .ok_or_else(|| anyhow!("{parent}: MT{} leaves no nucleus", spec.reaction.mt))?;
        let light: Vec<Product> = match spec.light {
            None => Vec::new(),
            Some((lz, la)) => {
                let nuclide = index.get((lz, la, 0)).ok_or_else(|| {
                    anyhow!(
                        "the light product of MT{} is not in the decay data",
                        spec.reaction.mt
                    )
                })?;
                vec![Product { nuclide, count: 1 }]
            }
        };
        let states = evaluation
            .final_states
            .iter()
            .find(|(mt, _)| *mt == spec.reaction.mt)
            .map(|(_, entries)| entries);
        let channels = match states {
            Some(entries) => {
                let (keys, notes) =
                    assign_states(parent, spec.reaction.mt, entries, (z, a), evaluations)?;
                if !notes.is_empty() {
                    warnings.push(format!(
                        "{parent}: MF8/MT{} final states of {}: {}",
                        spec.reaction.mt,
                        name((z, a, 0))?,
                        notes.join(", ")
                    ));
                }
                let mut targets: Vec<(FinalState, Key)> = entries
                    .iter()
                    .zip(keys)
                    .map(|(e, k)| {
                        (
                            FinalState {
                                mf: e.mf,
                                lfs: Some(e.lfs),
                            },
                            k,
                        )
                    })
                    .collect();
                // 基底状態を挙げない反応は、挙げた状態を除いた残りが基底状態へ行く
                if entries.iter().all(|e| e.excitation > 0.0) {
                    let mf = entries[0].mf;
                    if entries.iter().any(|e| e.mf != mf) {
                        bail!("{parent}: MF8/MT{} mixes MF9 and MF10", spec.reaction.mt);
                    }
                    targets.push((FinalState { mf, lfs: None }, (z, a, 0)));
                }
                let mut channels = Vec::new();
                for (final_state, k) in targets {
                    let mut products = Vec::new();
                    match index.get(k) {
                        Some(i) => products.push(Product { nuclide: i, count: 1 }),
                        None if z == 0 => {}
                        None => warnings.push(format!(
                            "{parent}: the product {} of MT{} is not in the decay data; the channel produces no residual",
                            name(k)?,
                            spec.reaction.mt
                        )),
                    }
                    products.extend(light.iter().copied());
                    channels.push(Channel {
                        final_state: Some(final_state),
                        products: merged(products),
                    });
                }
                channels
            }
            None => {
                let mut products = Vec::new();
                match index.get((z, a, 0)) {
                    Some(i) => products.push(Product { nuclide: i, count: 1 }),
                    // 残りが中性子だけ（Z = 0。H-2 の (n,p)）なら、β−n の中性子と同じく追わない
                    None if z == 0 => {}
                    None => warnings.push(format!(
                        "{parent}: the product {} of MT{} is not in the decay data; the reaction produces no residual",
                        name((z, a, 0))?,
                        spec.reaction.mt
                    )),
                }
                products.extend(light.iter().copied());
                vec![Channel {
                    final_state: None,
                    products: merged(products),
                }]
            }
        };
        reactions.push(Reaction {
            mt: spec.reaction.mt,
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
    if let Some(w) = yield_evaluations
        .windows(2)
        .find(|w| key(&w[0].material) == key(&w[1].material))
    {
        bail!(
            "{} appears twice in the fission yield data",
            name(key(&w[0].material))?
        );
    }
    for y in &yield_evaluations {
        if index.get(key(&y.material)).is_none() {
            bail!(
                "fission yields for {}, which is not in the decay data",
                name(key(&y.material))?
            );
        }
    }
    let mut fissionable = Vec::new();
    let mut seen = std::collections::HashSet::new();
    // 中性子そのものの評価（n + n）もチェーンに入れない
    for evaluation in neutrons.iter().filter(|n| n.material.z > 0) {
        let k = key(&evaluation.material);
        let i = index.get(k).ok_or_else(|| {
            anyhow!(
                "neutron data for {}, which is not in the decay data",
                name(k).unwrap_or_default()
            )
        })?;
        if !seen.insert(i) {
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
    // nfy の評価は核分裂する核種のもの（自分の収率になり、借りる先にもなる）
    for y in &yield_evaluations {
        let k = key(&y.material);
        if !fissionable.iter().any(|&i| index.get(k) == Some(i)) {
            bail!(
                "fission yields for {}, which has no fission in the neutron data",
                name(k)?
            );
        }
    }
    // 収率の表は評価ごとに 1 度だけ作り（警告も 1 度）、借りる核種にはコピーする
    let mut tables: HashMap<usize, Fission> = HashMap::new();
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
        let fission = match tables.get(&source) {
            Some(f) => f.clone(),
            None => {
                let f = fission_yields(&source_name, evaluation, &index, &nuclides, &mut warnings)?;
                tables.insert(source, f.clone());
                f
            }
        };
        nuclides[i].fission = Some(fission);
    }

    let chain = Chain { nuclides };
    chain.check()?;
    Ok(Built { chain, warnings })
}

#[cfg(test)]
mod tests;
