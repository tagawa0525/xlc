//! 格子（ピンセル）の断面積を、指定した群で OpenMC の実効断面積に置き換えて解き、k の変化を測る（同じ組成の OpenMC との k の
//! 差を、核種、物質、エネルギーの範囲ごとの断面積の差に分ける。#0133）。
//!
//! 使い方: swap_absorption <格子の入力 (TOML)> <OpenMC の openmc-rates.msgpack> <置き換えの一覧（テキスト）>
//!     [--materials <物質の吸収のある OpenMC の出力>] [--transfer <OpenMC の移行の出力>] [--chi <OpenMC の χ>]
//!
//! 一覧は 1 行に 1 つの置き換え（空行と `#` からの行は読まない）。置き換えは次の部分を `,` でつないだもので、置き換えごとに別に
//! 解き、置き換えない k との差を書く。群は 1 始まりで、高いエネルギーから数え、`<最初の群>-<最後の群>` は両端を含む。
//!
//! - `<核種>:<最初の群>-<最後の群>`: 燃料の核種の吸収と核分裂。核種を `*` にすると、OpenMC の反応率のある燃料の全核種
//! - `@<物質>=<OpenMC の物質>:<最初の群>-<最後の群>`: 物質（被覆管、減速材など。xlc の入力の名前）の吸収（`--materials`。
//!   openmc_pin_rates.py の出力の `materials` の `flux` と `absorption`）
//! - `~<物質>:<最初の群>-<最後の群>`: 物質の散乱の移行の行（`--transfer`。openmc_pin_transfer.py の出力の `materials`。物質の
//!   名前は xlc と同じ）
//! - `chi`: 燃料の核分裂の項の χ（`--chi`。openmc_chi.py の出力の `chi` にある核種の項だけ）
//!
//! 燃料の核種の OpenMC の実効断面積は、燃料の群ごとの反応率（`(n,gamma)`、`fission`、`nu-fission`）/ 燃料の中性子束 / 数密度
//! （openmc_pin_rates.py などの出力。ピンの燃料の平均）。xlc の側は、最初の燃料のピンの、ピンの燃料で平均した物理的な実効断面積。
//! 置き換えるのは捕獲（MT102）+ 核分裂と、核分裂と νσf で、両方の側で同じ反応にそろえる（OpenMC の `absorption` は (n,p)、(n,α)
//! などを含み、xlc の捕獲は含まないので使わない）。非弾性散乱、(n,2n)、(n,p)、(n,α) などは置き換えず、高速の群でもそのまま残る。
//! 置き換えは比 R = OpenMC / xlc を、燃料のリングごとの核種の値に掛ける形にし（リングの間の違いを残す）、リングの物質の全断面積、
//! 核分裂、νΣf と、その核種の核分裂の項を、SPH 因子（物質の全断面積 / 核種の全断面積の和。SPH 因子を掛けた物質の全断面積は核種の
//! 和の μ 倍）を掛けて動かす。散乱は変えない（吸収の差だけを測る）。
//!
//! 物質の吸収は、OpenMC の `absorption`（消滅の反応のすべて: 捕獲、(n,p)、(n,α) など）/ 物質の中性子束と、xlc の物質全体を物理的な
//! 中性子束で平均した吸収の比 R を、物質のすべての領域の吸収に掛ける。xlc の物質の吸収 Σt − Σ_g' Σs(g→g') は、散乱行列が
//! (n,2n) などの多重度を数えるので、消滅の吸収から中性子の過剰（(n,2n) の σ、(n,3n) の 2σ、(n,4n) の 3σ。ライブラリの燃焼チェーンの
//! 反応の無限希釈の群の断面積）を引いた値になる（被覆管の Zr-91 は 6〜10 MeV の群で (n,2n) が消滅の吸収の数倍で、正味の吸収が
//! 負になりうる）。比べるのは過剰を足し戻した値で、全断面積を (R − 1) ×（正味の吸収 + 過剰）だけ動かし、散乱は変えない。ライブラリが
//! 群の断面積を持たない多重度の反応（MT11、24、25、30、41、42）は過剰に入らない。
//!
//! 移行は `swap_transfer` の「SPH kept」と同じ置き換え（物質全体の移行の断面積の要素ごとの比 R を、物質のすべての領域の行に掛け、
//! xlc が 0 で OpenMC が正の要素は OpenMC の値を置く。P1 以上の行列に同じ比を掛け、全断面積は行の和の変化だけ動かす）。χ は、燃料の
//! リングの物質の核分裂の項のうち `chi` にある核種の χ を置き換える（νΣf は変えない）。
//!
//! 物質の吸収と移行の比は、置き換えない xlc の解の反応率で求め、同じ行の別の部分の置き換えで変えない（置き換えどうしの相互作用を、
//! 別々に解いた和と比べて測るため）。物質の吸収と移行は、SPH 因子を外した物理的な物質を置き換えてから、同じ因子を掛け直す。
//! 外周が鏡面反射の格子だけを受け付ける（ピンセルの比較のため）。燃料のピンが複数あれば、どのピンの燃料にも同じ比を掛ける。
//!
//! SPH 因子は置き換える前のものをそのまま使い、置き換えた断面積で求め直さない（`swap_transfer` の (A) と同じ）。UFG の範囲の
//! 群では、置き換えでリングの間の中性子束の形が変わる分を SPH 因子が追わないので、その分は k の変化に入らない。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

use std::collections::BTreeMap;
use std::ops::Range;
use std::path::Path;

use anyhow::{Context, anyhow};
use serde::Deserialize;
use serde::de::DeserializeOwned;
use xlc_format::library::Library;
use xlc_format::msgpack::{from_value, read_value};
use xlc_lattice::balance::by_material;
use xlc_lattice::input::{CrossSections, LatticeInput, MaterialInput};
use xlc_lattice::material::Material;
use xlc_lattice::moc::{Options, Solution, available_order, solve};
use xlc_lattice::shielding::{FuelShielding, shield_fuel};
use xlc_lattice::tracking::{Boundaries, Tracking};

#[derive(Deserialize)]
struct Score {
    mean: Vec<f64>,
}

#[derive(Deserialize)]
struct Reference {
    /// 群の境界 (eV)、降順
    bounds: Vec<f64>,
    flux: Vec<f64>,
    densities: BTreeMap<String, f64>,
    rates: BTreeMap<String, BTreeMap<String, Score>>,
}

/// OpenMC の物質の群ごとの中性子束と吸収（openmc_pin_rates.py の `materials`）
#[derive(Deserialize)]
struct MaterialRates {
    flux: Vec<f64>,
    absorption: Vec<f64>,
}

#[derive(Deserialize)]
struct MaterialReference {
    bounds: Vec<f64>,
    materials: BTreeMap<String, MaterialRates>,
}

/// OpenMC の物質の群ごとの中性子束と移行（openmc_pin_transfer.py の `materials`）
#[derive(Deserialize)]
struct TransferRates {
    flux: Vec<f64>,
    transfer: Vec<Vec<f64>>,
}

#[derive(Deserialize)]
struct TransferReference {
    bounds: Vec<f64>,
    materials: BTreeMap<String, TransferRates>,
}

/// OpenMC の核種ごとの χ（openmc_chi.py の出力）
#[derive(Deserialize)]
struct ChiReference {
    bounds: Vec<f64>,
    chi: BTreeMap<String, Vec<f64>>,
}

/// 1 つの核種の、群ごとの吸収と核分裂の比 OpenMC / xlc（置き換えない群は 1）
struct Ratio {
    absorption: Vec<f64>,
    fission: Vec<f64>,
    nu_fission: Vec<f64>,
}

fn ratio(
    reference: &Reference,
    name: &str,
    xlc: &xlc_format::library::Shielded,
    groups: &Range<usize>,
    ng: usize,
) -> anyhow::Result<Ratio> {
    let rates = reference
        .rates
        .get(name)
        .ok_or_else(|| anyhow!("{name} has no OpenMC rates"))?;
    let score = |s: &str| {
        rates
            .get(s)
            .map(|v| &v.mean)
            .ok_or_else(|| anyhow!("{name} has no OpenMC `{s}`"))
    };
    let (capture, fission, nu_fission) =
        (score("(n,gamma)")?, score("fission")?, score("nu-fission")?);
    let density = *reference
        .densities
        .get(name)
        .ok_or_else(|| anyhow!("{name} has no OpenMC density"))?;
    if [capture, fission, nu_fission, &reference.flux]
        .iter()
        .any(|v| v.len() != ng)
    {
        return Err(anyhow!(
            "the OpenMC rates of {name} do not have {ng} groups"
        ));
    }
    let mut r = Ratio {
        absorption: vec![1.0; ng],
        fission: vec![1.0; ng],
        nu_fission: vec![1.0; ng],
    };
    for g in groups.clone() {
        let phi = reference.flux[g] * density;
        let a_x = xlc.capture[g] + xlc.fission[g];
        if !(phi > 0.0 && a_x > 0.0) {
            return Err(anyhow!(
                "{name} group {}: no flux or no absorption to compare",
                g + 1
            ));
        }
        r.absorption[g] = (capture[g] + fission[g]) / phi / a_x;
        if xlc.fission[g] > 0.0 {
            if !xlc.nu_fission[g].is_finite() || xlc.nu_fission[g] <= 0.0 {
                return Err(anyhow!("{name} group {}: fission without ν", g + 1));
            }
            r.fission[g] = fission[g] / phi / xlc.fission[g];
            r.nu_fission[g] = nu_fission[g] / phi / xlc.nu_fission[g];
        } else if fission[g] > 0.0 {
            return Err(anyhow!(
                "{name} group {}: OpenMC fissions but xlc does not",
                g + 1
            ));
        }
    }
    Ok(r)
}

/// 物質の吸収の比 R = OpenMC / xlc（`groups` の外は 1）。`xlc_absorption` は xlc の物質全体の、多重度を差し引かない吸収
/// （正味の吸収 + 中性子の過剰）の物理的な値
fn material_ratio(
    openmc_absorption: &[f64],
    openmc_flux: &[f64],
    xlc_absorption: &[f64],
    groups: &Range<usize>,
) -> anyhow::Result<Vec<f64>> {
    let ng = xlc_absorption.len();
    if openmc_absorption.len() != ng || openmc_flux.len() != ng || groups.end > ng {
        return Err(anyhow!("the OpenMC absorption does not have {ng} groups"));
    }
    let mut r = vec![1.0; ng];
    for g in groups.clone() {
        if !(openmc_flux[g] > 0.0 && xlc_absorption[g] > 0.0) {
            return Err(anyhow!(
                "group {}: no flux or no absorption to compare",
                g + 1
            ));
        }
        r[g] = openmc_absorption[g] / openmc_flux[g] / xlc_absorption[g];
    }
    Ok(r)
}

/// 物理的な物質 `m` の、多重度を差し引かない吸収（正味の吸収 + `excess`）に比 `ratio` を掛けた物質（全断面積だけを動かす）
fn with_absorption(m: &Material, ratio: &[f64], excess: &[f64]) -> Material {
    let mut m = m.clone();
    let absorption = m.absorption();
    for (g, total) in m.total.iter_mut().enumerate() {
        *total += (ratio[g] - 1.0) * (absorption[g] + excess[g]);
    }
    m
}

/// 物質 `m` の核分裂の項（`terms` は項ごとの核種の名前）のうち、`chi` にある核種の χ を置き換える
fn replace_chi(
    m: &mut Material,
    terms: &[String],
    chi: &BTreeMap<String, Vec<f64>>,
) -> anyhow::Result<()> {
    if terms.len() != m.fission_terms.len() {
        return Err(anyhow!(
            "{} nuclide names for {} fission terms",
            terms.len(),
            m.fission_terms.len()
        ));
    }
    for (name, values) in chi {
        if values.len() != m.groups() {
            return Err(anyhow!(
                "the chi of {name} does not have {} groups",
                m.groups()
            ));
        }
        let t = terms
            .iter()
            .position(|n| n == name)
            .ok_or_else(|| anyhow!("{name} has no fission term"))?;
        m.fission_terms[t].chi = values.clone();
    }
    m.validate()?;
    Ok(())
}

/// 物質の中性子の過剰 Σ_核種 N (σ(n,2n) + 2σ(n,3n) + 3σ(n,4n)) (1/cm)。ライブラリの燃焼チェーンの反応の無限希釈の群の断面積
fn neutron_excess(
    library: &Library,
    material: &MaterialInput,
    ng: usize,
) -> anyhow::Result<Vec<f64>> {
    let mut excess = vec![0.0; ng];
    for (name, density) in &material.composition {
        let (_, table) = library.lookup(name, material.temperature)?;
        for reaction in &table.reactions {
            let extra = match reaction.mt {
                16 => 1.0,
                17 => 2.0,
                37 => 3.0,
                _ => continue,
            };
            for channel in &reaction.channels {
                for (e, s) in excess.iter_mut().zip(&channel.cross_section) {
                    *e += extra * density * s;
                }
            }
        }
    }
    Ok(excess)
}

/// 物質全体の移行の断面積 [入射][出射] (1/cm)
fn transfer_xs(flux: &[f64], transfer: &[Vec<f64>]) -> anyhow::Result<Vec<Vec<f64>>> {
    transfer
        .iter()
        .zip(flux)
        .map(|(row, &phi)| {
            if phi > 0.0 {
                Ok(row.iter().map(|t| t / phi).collect())
            } else {
                Err(anyhow!("no flux"))
            }
        })
        .collect()
}

/// 物理的な物質 `m` の行 `rows` に比 `ratio` を掛ける（xlc が 0 の要素は `reference` の値を置く）。全断面積は行の和の
/// 変化だけ動かす（`swap_transfer` の `swap` と同じ）
fn swap_rows(
    m: &Material,
    rows: &[bool],
    ratio: &[Vec<Option<f64>>],
    reference: &[Vec<f64>],
) -> Material {
    let mut m = m.clone();
    for g in (0..rows.len()).filter(|&g| rows[g]) {
        let before: f64 = m.scatter[g].iter().sum();
        for (h, r) in ratio[g].iter().enumerate() {
            match r {
                Some(r) => {
                    m.scatter[g][h] *= r;
                    for moment in &mut m.anisotropic {
                        moment[g][h] *= r;
                    }
                }
                None => m.scatter[g][h] = reference[g][h],
            }
        }
        m.total[g] += m.scatter[g].iter().sum::<f64>() - before;
    }
    m
}

/// 1 つの部分の置き換え
enum Swap {
    /// 燃料の核種の吸収と核分裂
    Nuclide(String, Range<usize>),
    /// 物質（xlc の名前）の吸収を、OpenMC の物質 `openmc` の値に
    Absorption {
        material: String,
        openmc: String,
        groups: Range<usize>,
    },
    /// 物質の散乱の移行の行
    Transfer {
        material: String,
        groups: Range<usize>,
    },
    /// 燃料の核分裂の項の χ
    Chi,
}

/// `<最初>-<最後>`（1 始まり、両端を含む）を、0 始まりで終わりを含まない範囲にする
fn groups_of(part: &str, range: &str, ng: usize) -> anyhow::Result<Range<usize>> {
    let (first, last) = range
        .split_once('-')
        .ok_or_else(|| anyhow!("`{part}`: expected <first>-<last>"))?;
    let (first, last): (usize, usize) = (first.parse()?, last.parse()?);
    if first == 0 || last < first || last > ng {
        return Err(anyhow!("`{part}`: groups count from 1 to {ng}"));
    }
    Ok(first - 1..last)
}

/// 置き換えの 1 行を部分の並びにする（燃料の核種の `*` は、OpenMC の反応率のある燃料の全核種に開く）
fn parse(
    spec: &str,
    fuel: &[String],
    reference: &Reference,
    ng: usize,
) -> anyhow::Result<Vec<Swap>> {
    let mut out = Vec::new();
    for part in spec.split(',') {
        if part == "chi" {
            out.push(Swap::Chi);
            continue;
        }
        let (name, range) = part
            .split_once(':')
            .ok_or_else(|| anyhow!("`{part}`: expected <name>:<first>-<last> or chi"))?;
        let groups = groups_of(part, range, ng)?;
        if let Some(pair) = name.strip_prefix('@') {
            let (material, openmc) = pair
                .split_once('=')
                .ok_or_else(|| anyhow!("`{part}`: expected @<material>=<OpenMC material>"))?;
            out.push(Swap::Absorption {
                material: material.to_string(),
                openmc: openmc.to_string(),
                groups,
            });
        } else if let Some(material) = name.strip_prefix('~') {
            out.push(Swap::Transfer {
                material: material.to_string(),
                groups,
            });
        } else if name == "*" {
            out.extend(
                fuel.iter()
                    .filter(|n| reference.rates.contains_key(*n))
                    .map(|n| Swap::Nuclide(n.clone(), groups.clone())),
            );
        } else if fuel.iter().any(|n| n == name) {
            out.push(Swap::Nuclide(name.to_string(), groups));
        } else {
            return Err(anyhow!("{name} is not in the fuel"));
        }
    }
    Ok(out)
}

/// 燃料のリングの物質（核分裂のある物質）の番号と、その核分裂の項の核種の名前（項の順）
fn fuel_terms(
    shielded: &FuelShielding,
    input: &LatticeInput,
    names: &[String],
) -> Vec<(usize, Vec<String>)> {
    let mut out: Vec<(usize, Vec<String>)> = Vec::new();
    for (region, &index) in shielded.region_material.iter().enumerate() {
        let Some(input_material) = input.materials.iter().find(|m| m.name == names[region]) else {
            continue;
        };
        let Some(nuclides) = shielded.region_cross_sections(region) else {
            continue;
        };
        let fissile = shielded.materials[index]
            .nu_fission
            .iter()
            .any(|&v| v > 0.0);
        if !fissile || out.iter().any(|(i, _)| *i == index) {
            continue;
        }
        let terms = input_material
            .composition
            .iter()
            .zip(nuclides)
            .filter(|(_, s)| s.nu_fission.iter().any(|&v| v > 0.0))
            .map(|((name, _), _)| name.clone())
            .collect();
        out.push((index, terms));
    }
    out
}

/// 燃料のリングの物質に、核種ごとの比を入れた物質の一覧
fn swapped(
    shielded: &FuelShielding,
    input: &LatticeInput,
    names: &[String],
    ratios: &BTreeMap<String, Ratio>,
) -> anyhow::Result<Vec<Material>> {
    let mut materials = shielded.materials.clone();
    let mut done = vec![false; materials.len()];
    for (region, &index) in shielded.region_material.iter().enumerate() {
        let Some(input_material) = input.materials.iter().find(|m| m.name == names[region]) else {
            continue;
        };
        let Some(nuclides) = shielded.region_cross_sections(region) else {
            continue;
        };
        let fissile = materials[index].nu_fission.iter().any(|&v| v > 0.0);
        if !fissile || done[index] {
            continue;
        }
        done[index] = true;
        let m = &mut materials[index];
        // SPH 因子: 置き換える前の物質の全断面積 / 核種の全断面積の和
        let sph: Vec<f64> = m
            .total
            .iter()
            .enumerate()
            .map(|(g, total)| {
                total
                    / input_material
                        .composition
                        .iter()
                        .zip(nuclides)
                        .map(|((_, n), s)| n * s.total[g])
                        .sum::<f64>()
            })
            .collect();
        let mut term = 0;
        for ((name, density), s) in input_material.composition.iter().zip(nuclides) {
            let has_term = s.nu_fission.iter().any(|&v| v > 0.0);
            if let Some(r) = ratios.get(name) {
                for (g, &mu) in sph.iter().enumerate() {
                    let da = mu * density * (r.absorption[g] - 1.0) * (s.capture[g] + s.fission[g]);
                    let df = mu * density * (r.fission[g] - 1.0) * s.fission[g];
                    let dnf = mu * density * (r.nu_fission[g] - 1.0) * s.nu_fission[g];
                    m.total[g] += da;
                    m.fission[g] += df;
                    m.nu_fission[g] += dnf;
                    if has_term {
                        m.fission_terms[term].nu_fission[g] += dnf;
                    }
                }
            }
            if has_term {
                term += 1;
            }
        }
        if term != m.fission_terms.len() {
            return Err(anyhow!(
                "material {index}: {term} fissile nuclides but {} fission terms",
                m.fission_terms.len()
            ));
        }
    }
    Ok(materials)
}

/// [ライブラリの群] の因子（`groups` の外は 1）
fn spread(mu: &[f64], groups: &Range<usize>, ng: usize) -> Vec<f64> {
    let mut all = vec![1.0; ng];
    all[groups.clone()].copy_from_slice(mu);
    all
}

/// 物質ごとの SPH 因子（[ライブラリの群]。因子のない物質は `None`）。燃料のピンの等価なセルのリングの物質だけが持つ
fn sph_of_materials(
    shielded: &FuelShielding,
    lattice: &xlc_lattice::geometry::Lattice,
    ng: usize,
) -> anyhow::Result<Vec<Option<Vec<f64>>>> {
    let mut factor_of: Vec<Option<Vec<f64>>> = vec![None; shielded.materials.len()];
    for pin in &shielded.pins {
        let Some(sph) = &pin.sph else { continue };
        if lattice
            .cell_type(pin.cell)
            .pin()
            .is_none_or(|p| p.frame().is_some())
        {
            return Err(anyhow!("cell {}: only pins without a frame", pin.cell));
        }
        let rings = pin
            .ring_regions(lattice)
            .ok_or_else(|| anyhow!("cell {}: no rings", pin.cell))?;
        for (ring, mu) in rings.iter().zip(&sph.factors) {
            let i = shielded.region_material[ring[0]];
            // リングの領域がすべて同じ物質（リングの因子を掛けたもの）であること
            if ring.iter().any(|&r| shielded.region_material[r] != i) {
                return Err(anyhow!("cell {}: a ring has several materials", pin.cell));
            }
            let all = spread(mu, &sph.groups, ng);
            match &factor_of[i] {
                Some(known) if known != &all => {
                    return Err(anyhow!(
                        "material {i} has different SPH factors in two pins"
                    ));
                }
                _ => factor_of[i] = Some(all),
            }
        }
    }
    Ok(factor_of)
}

/// 群の境界がライブラリと同じこと。いちばん上の境界は、OpenMC の tally が群 1 に上の中性子を数えるために高くしてよい
fn check_bounds(what: &str, given: &[f64], library: &[f64]) -> anyhow::Result<()> {
    let same = |a: f64, b: f64| (a - b).abs() <= 1e-9 * b.abs();
    if given.len() != library.len()
        || given[0] < library[0]
        || !given[1..]
            .iter()
            .zip(&library[1..])
            .all(|(&a, &b)| same(a, b))
    {
        return Err(anyhow!(
            "the groups of {what} {given:?} are not those of the library"
        ));
    }
    Ok(())
}

fn read_reference<T: DeserializeOwned>(path: &str) -> anyhow::Result<T> {
    Ok(from_value(read_value(
        &std::fs::read(path).with_context(|| path.to_string())?,
    )?)?)
}

/// 1 つの物質（xlc の名前）の、物理的な物質への置き換え
struct MaterialChange {
    /// 吸収の比（置き換えない群は 1）と中性子の過剰
    absorption: Option<(Vec<f64>, Vec<f64>)>,
    /// 置き換える移行の行
    rows: Vec<bool>,
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let usage = || {
        anyhow!(
            "usage: swap_absorption <lattice input> <openmc-rates.msgpack> <list of swaps> \
             [--materials <OpenMC rates with materials>] [--transfer <OpenMC transfers>] \
             [--chi <OpenMC chi>]"
        )
    };
    let [_, input, reference, list, options @ ..] = args.as_slice() else {
        return Err(usage());
    };
    let (mut materials_path, mut transfer_path, mut chi_path) = (None, None, None);
    for pair in options.chunks(2) {
        let [flag, value] = pair else {
            return Err(usage());
        };
        let slot = match flag.as_str() {
            "--materials" => &mut materials_path,
            "--transfer" => &mut transfer_path,
            "--chi" => &mut chi_path,
            _ => return Err(usage()),
        };
        if slot.replace(value.as_str()).is_some() {
            return Err(anyhow!("{flag} is given twice"));
        }
    }
    let specs: Vec<String> = std::fs::read_to_string(list)
        .with_context(|| list.clone())?
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(String::from)
        .collect();
    if specs.is_empty() {
        return Err(anyhow!("{list}: no swaps"));
    }
    let path = Path::new(input);
    let input = LatticeInput::parse(&std::fs::read_to_string(path)?).context("input")?;
    let CrossSections::Library(file) = &input.cross_sections else {
        return Err(anyhow!("the input must use a library"));
    };
    if !input.self_shielding {
        return Err(anyhow!("the input must be self-shielded"));
    }
    let base = path.parent().unwrap_or(Path::new("."));
    let library = Library::from_msgpack(&std::fs::read(base.join(file))?)?;
    let reference: Reference = read_reference(reference)?;
    let ng = library.groups.len();
    let bounds = library.groups.bounds();
    check_bounds("the OpenMC rates", &reference.bounds, bounds)?;
    let material_reference: Option<MaterialReference> =
        materials_path.map(read_reference).transpose()?;
    let transfer_reference: Option<TransferReference> =
        transfer_path.map(read_reference).transpose()?;
    let chi_reference: Option<ChiReference> = chi_path.map(read_reference).transpose()?;
    if let Some(r) = &material_reference {
        check_bounds("the OpenMC materials", &r.bounds, bounds)?;
    }
    if let Some(r) = &transfer_reference {
        check_bounds("the OpenMC transfers", &r.bounds, bounds)?;
    }
    if let Some(r) = &chi_reference {
        check_bounds("the OpenMC chi", &r.bounds, bounds)?;
    }
    let (lattice, names) = input.geometry()?;
    let tracking =
        Tracking::with_boundaries(&lattice, input.azimuthal, input.spacing, input.boundaries)?;
    if tracking.boundaries() != Boundaries::REFLECTIVE {
        return Err(anyhow!("the lattice must have reflective sides"));
    }
    let shielded = shield_fuel(
        &library,
        &tracking,
        &input.materials,
        &names,
        &input.shielding_options(),
    )?;
    let pin = shielded
        .pins
        .first()
        .ok_or_else(|| anyhow!("no fuel pin"))?;
    let fuel: Vec<String> = pin.nuclides.iter().map(|(n, _)| n.name.clone()).collect();
    let order = input.scattering_order;
    let solve_regions = |regions: &[&Material]| -> anyhow::Result<Solution> {
        let options = Options {
            scattering_order: order.unwrap_or_else(|| available_order(regions)),
            ..Options::default()
        };
        Ok(solve(&tracking, regions, &options)?)
    };
    let regions: Vec<&Material> = shielded
        .region_material
        .iter()
        .map(|&i| &shielded.materials[i])
        .collect();
    let solution = solve_regions(&regions)?;
    let k0 = solution.k;
    println!("k {k0:.6}");
    // 物質の吸収と移行の比に使う、置き換えない解の物質ごとの反応率（物理的な中性子束で平均）
    let physical = shielded.physical_flux(&tracking.lattice, &solution.flux)?;
    let rates = by_material(&solution, &regions, &names, &physical)?;
    let factor_of = sph_of_materials(&shielded, &lattice, ng)?;
    let terms = fuel_terms(&shielded, &input, &names);
    // 物質（xlc の名前）の番号。ほかの名前の領域と共有する物質は置き換えられない
    let indices_of = |material: &str| -> anyhow::Result<Vec<usize>> {
        let mut chosen: Vec<usize> = Vec::new();
        for (region, &i) in shielded.region_material.iter().enumerate() {
            if names[region] == material && !chosen.contains(&i) {
                chosen.push(i);
            }
        }
        if chosen.is_empty() {
            return Err(anyhow!("`{material}` is not in the lattice"));
        }
        for (region, &i) in shielded.region_material.iter().enumerate() {
            if chosen.contains(&i) && names[region] != material {
                return Err(anyhow!(
                    "material {i} is shared by `{material}` and `{}`",
                    names[region]
                ));
            }
        }
        Ok(chosen)
    };
    let rates_of = |material: &str| {
        rates
            .get(material)
            .ok_or_else(|| anyhow!("no xlc rates of `{material}`"))
    };
    for spec in &specs {
        let mut ratios: BTreeMap<String, Ratio> = BTreeMap::new();
        let mut changes: BTreeMap<String, MaterialChange> = BTreeMap::new();
        let mut chi = false;
        for part in parse(spec, &fuel, &reference, ng)? {
            match part {
                Swap::Nuclide(name, groups) => {
                    let k = fuel
                        .iter()
                        .position(|n| *n == name)
                        .ok_or_else(|| anyhow!("{name} is not in the fuel"))?;
                    let r = ratio(&reference, &name, &pin.effective[k], &groups, ng)?;
                    match ratios.get_mut(&name) {
                        // 同じ核種の別の範囲は、置き換える群を合わせる
                        Some(existing) => {
                            for g in groups {
                                existing.absorption[g] = r.absorption[g];
                                existing.fission[g] = r.fission[g];
                                existing.nu_fission[g] = r.nu_fission[g];
                            }
                        }
                        None => {
                            ratios.insert(name, r);
                        }
                    }
                }
                Swap::Absorption {
                    material,
                    openmc,
                    groups,
                } => {
                    let o = material_reference
                        .as_ref()
                        .ok_or_else(|| anyhow!("`{spec}` needs --materials"))?
                        .materials
                        .get(&openmc)
                        .ok_or_else(|| anyhow!("no OpenMC material `{openmc}`"))?;
                    let input_material = input
                        .materials
                        .iter()
                        .find(|m| m.name == material)
                        .ok_or_else(|| anyhow!("no input material `{material}`"))?;
                    let excess = neutron_excess(&library, input_material, ng)?;
                    let x = rates_of(&material)?;
                    let gross: Vec<f64> = x
                        .absorption
                        .iter()
                        .zip(&x.flux)
                        .zip(&excess)
                        .map(|((a, f), e)| a / f + e)
                        .collect();
                    let r = material_ratio(&o.absorption, &o.flux, &gross, &groups)
                        .with_context(|| format!("`{material}` and OpenMC `{openmc}`"))?;
                    let change = changes.entry(material).or_insert_with(|| MaterialChange {
                        absorption: None,
                        rows: vec![false; ng],
                    });
                    match &mut change.absorption {
                        // 同じ物質の別の範囲は、置き換える群を合わせる
                        Some((existing, _)) => existing[groups.clone()].copy_from_slice(&r[groups]),
                        None => change.absorption = Some((r, excess)),
                    }
                }
                Swap::Transfer { material, groups } => {
                    if transfer_reference.is_none() {
                        return Err(anyhow!("`{spec}` needs --transfer"));
                    }
                    let change = changes.entry(material).or_insert_with(|| MaterialChange {
                        absorption: None,
                        rows: vec![false; ng],
                    });
                    for g in groups {
                        change.rows[g] = true;
                    }
                }
                Swap::Chi => {
                    if chi_reference.is_none() {
                        return Err(anyhow!("`{spec}` needs --chi"));
                    }
                    chi = true;
                }
            }
        }
        let mut materials = swapped(&shielded, &input, &names, &ratios)?;
        for (material, change) in &changes {
            // 移行の比 R = OpenMC / xlc（xlc が 0 なら None）と OpenMC の移行の断面積
            let transfer = match &transfer_reference {
                Some(reference) if change.rows.iter().any(|&r| r) => {
                    let o = reference
                        .materials
                        .get(material)
                        .ok_or_else(|| anyhow!("no OpenMC transfers of `{material}`"))?;
                    let openmc = transfer_xs(&o.flux, &o.transfer)?;
                    let x = rates_of(material)?;
                    let xlc = transfer_xs(&x.flux, &x.transfer)?;
                    if openmc.len() != ng || openmc.iter().any(|r| r.len() != ng) {
                        return Err(anyhow!("the OpenMC transfers do not have {ng} groups"));
                    }
                    let ratio: Vec<Vec<Option<f64>>> = xlc
                        .iter()
                        .zip(&openmc)
                        .map(|(x, o)| {
                            x.iter()
                                .zip(o)
                                .map(|(&x, &o)| (x > 0.0).then(|| o / x))
                                .collect()
                        })
                        .collect();
                    Some((ratio, openmc))
                }
                _ => None,
            };
            for i in indices_of(material)? {
                let mu = &factor_of[i];
                let mut m = match mu {
                    Some(mu) => {
                        materials[i].scaled(&mu.iter().map(|m| 1.0 / m).collect::<Vec<_>>())?
                    }
                    None => materials[i].clone(),
                };
                if let Some((ratio, excess)) = &change.absorption {
                    m = with_absorption(&m, ratio, excess);
                }
                if let Some((ratio, openmc)) = &transfer {
                    m = swap_rows(&m, &change.rows, ratio, openmc);
                }
                materials[i] = match mu {
                    Some(mu) => m.scaled(mu)?,
                    None => m,
                };
            }
        }
        if let (true, Some(reference)) = (chi, &chi_reference) {
            for (index, names) in &terms {
                replace_chi(&mut materials[*index], names, &reference.chi)?;
            }
        }
        let regions: Vec<&Material> = shielded
            .region_material
            .iter()
            .map(|&i| &materials[i])
            .collect();
        let k = solve_regions(&regions)?.k;
        println!("{spec}: k {k:.6} ({:+.1} pcm)", (k - k0) * 1e5);
    }
    Ok(())
}
#[cfg(test)]
#[path = "swap_absorption/tests.rs"]
mod tests;
