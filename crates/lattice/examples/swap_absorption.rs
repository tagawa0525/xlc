//! 格子（ピンセル）の燃料の核種の吸収と核分裂の断面積を、指定した群で OpenMC の実効断面積に置き換えて解き、k の変化を測る
//! （同じ組成の OpenMC との k の差を、核種とエネルギーの範囲ごとの断面積の差に分ける。#0133）。
//!
//! 使い方: swap_absorption <格子の入力 (TOML)> <OpenMC の openmc-rates.msgpack> <置き換えの一覧（テキスト）>
//!
//! 一覧は 1 行に 1 つの置き換え（空行と `#` からの行は読まない）。置き換えは `<核種>:<最初の群>-<最後の群>` を `,` でつないだもの
//! （群は 1 始まりで、高いエネルギーから数える）。核種を `*` にすると、OpenMC の反応率のある燃料の全核種。置き換えごとに別に解き、
//! 置き換えない k との差を書く。
//!
//! OpenMC の実効断面積は、燃料の群ごとの反応率（`(n,gamma)`、`fission`、`nu-fission`）/ 燃料の中性子束 / 数密度（openmc_pin_rates.py
//! などの出力。ピンの燃料の平均）。xlc の側は、最初の燃料のピンの、ピンの燃料で平均した物理的な実効断面積。置き換えるのは捕獲
//! （MT102）+ 核分裂と、核分裂と νσf で、両方の側で同じ反応にそろえる（OpenMC の `absorption` は (n,p)、(n,α) などを含み、xlc の
//! 捕獲は含まないので使わない）。非弾性散乱、(n,2n)、(n,p)、(n,α) などは置き換えず、高速の群でもそのまま残る。
//! 置き換えは比 R = OpenMC / xlc を、燃料のリングごとの核種の値に掛ける形にし（リングの間の違いを残す）、リングの物質の全断面積、
//! 核分裂、νΣf と、その核種の核分裂の項を、SPH 因子（物質の全断面積 / 核種の全断面積の和。SPH 因子を掛けた物質の全断面積は核種の
//! 和の μ 倍）を掛けて動かす。散乱は変えない（吸収の差だけを測る）。
//! 外周が鏡面反射の格子だけを受け付ける（ピンセルの比較のため）。燃料のピンが複数あれば、どのピンの燃料にも同じ比を掛ける。
//!
//! SPH 因子は置き換える前のものをそのまま使い、置き換えた断面積で求め直さない（`swap_transfer` の (A) と同じ）。UFG の範囲の
//! 群では、置き換えでリングの間の中性子束の形が変わる分を SPH 因子が追わないので、その分は k の変化に入らない。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, anyhow};
use serde::Deserialize;
use xlc_format::library::Library;
use xlc_format::msgpack::{from_value, read_value};
use xlc_lattice::input::{CrossSections, LatticeInput};
use xlc_lattice::material::Material;
use xlc_lattice::moc::{Options, available_order, solve};
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
    groups: &std::ops::Range<usize>,
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

/// `spec`（`<核種>:<最初>-<最後>` を `,` でつないだもの）を、[核種] の群の範囲（0 始まり、終わりを含まない）にする
fn parse(
    spec: &str,
    fuel: &[String],
    reference: &Reference,
    ng: usize,
) -> anyhow::Result<Vec<(String, std::ops::Range<usize>)>> {
    let mut out = Vec::new();
    for part in spec.split(',') {
        let (name, range) = part
            .split_once(':')
            .ok_or_else(|| anyhow!("`{part}`: expected <nuclide>:<first>-<last>"))?;
        let (first, last) = range
            .split_once('-')
            .ok_or_else(|| anyhow!("`{part}`: expected <first>-<last>"))?;
        let (first, last): (usize, usize) = (first.parse()?, last.parse()?);
        if first == 0 || last < first || last > ng {
            return Err(anyhow!("`{part}`: groups count from 1 to {ng}"));
        }
        let names: Vec<String> = if name == "*" {
            fuel.iter()
                .filter(|n| reference.rates.contains_key(*n))
                .cloned()
                .collect()
        } else if fuel.iter().any(|n| n == name) {
            vec![name.to_string()]
        } else {
            return Err(anyhow!("{name} is not in the fuel"));
        };
        out.extend(names.into_iter().map(|n| (n, first - 1..last)));
    }
    Ok(out)
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

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let [_, input, reference, list] = args.as_slice() else {
        return Err(anyhow!(
            "usage: swap_absorption <lattice input> <openmc-rates.msgpack> <list of swaps>"
        ));
    };
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
    let reference: Reference = from_value(read_value(&std::fs::read(reference)?)?)?;
    let ng = library.groups.len();
    // 群の境界がライブラリと同じこと。いちばん上の境界は、OpenMC の tally が群 1 に上の中性子を数えるために高くしてよい
    let bounds = library.groups.bounds();
    let same = |a: f64, b: f64| (a - b).abs() <= 1e-9 * b.abs();
    if reference.bounds.len() != bounds.len()
        || reference.bounds[0] < bounds[0]
        || !reference.bounds[1..]
            .iter()
            .zip(&bounds[1..])
            .all(|(&a, &b)| same(a, b))
    {
        return Err(anyhow!(
            "the OpenMC groups {:?} are not those of the library",
            reference.bounds
        ));
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
    let solve_k = |materials: &[Material]| -> anyhow::Result<f64> {
        let regions: Vec<&Material> = shielded
            .region_material
            .iter()
            .map(|&i| &materials[i])
            .collect();
        let options = Options {
            scattering_order: input
                .scattering_order
                .unwrap_or_else(|| available_order(&regions)),
            ..Options::default()
        };
        Ok(solve(&tracking, &regions, &options)?.k)
    };
    let k0 = solve_k(&shielded.materials)?;
    println!("k {k0:.6}");
    for spec in &specs {
        let mut ratios: BTreeMap<String, Ratio> = BTreeMap::new();
        for (name, groups) in parse(spec, &fuel, &reference, ng)? {
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
        let materials = swapped(&shielded, &input, &names, &ratios)?;
        let k = solve_k(&materials)?;
        println!("{spec}: k {k:.6} ({:+.1} pcm)", (k - k0) * 1e5);
    }
    Ok(())
}

#[cfg(test)]
#[path = "swap_absorption/tests.rs"]
mod tests;
