//! 格子（ピンセルか集合体）の 1 つの物質の散乱行列のうち、指定した入射の群の行を OpenMC の行列に置き換えて解き、k の変化を
//! 測る（OpenMC との k の差のうち、群の間の散乱による移行の差による分を見積もる。examples/vera/verification.md「k の差の
//! 内訳」）。
//!
//! 使い方: swap_transfer <格子の入力 (TOML)> <OpenMC の balance.msgpack> <物質> <最初の群> <最後の群>
//!
//! 群は 1 始まりで、高いエネルギーから数える。OpenMC の行列は `openmc_balance.py` の物質の移行の率 / 中性子束（物質
//! 全体の平均）。ピンごとの等価なセルのリングと、それ以外の領域（案内管の中の水など）に物理的な（SPH 因子を掛けていない）
//! 行列があるので、置き換えは要素ごとの比 R = OpenMC / xlc（xlc の側も物質全体を物理的な中性子束で平均した値）を、その
//! 物質のすべての領域の行に掛ける形にし、リングとピンの間の違いを残す。xlc が 0 で OpenMC が正の要素（自由ガスの上方
//! 散乱）は、OpenMC の値をどの領域にも置く。P1 以上の行列は同じ比を掛ける（移行ごとの平均の余弦を保つ）。全断面積は行の和の変化だけ動かし、正味の吸収と生成は変えない。
//!
//! SPH 因子は 2 通り: (A) xlc 自身の行列で求めた因子をそのまま使う、(B) 置き換えた行列で、ピンごとに等価なセルの因子を
//! 求め直す（対称操作で移り合うピンは物質と因子を共有するので、同じ因子になることを確かめる）。
//! (B) の参照の中性子束（UFG）は保存していないので、元の因子を掛けた物質で等価なセルを多群で解いた中性子束に因子を
//! 掛けて戻す（収束した SPH 因子では φ_ref ∝ μ φ̃）。元の物質でその参照から因子を求め直し、元の因子と合うことを確かめる。

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, anyhow};
use serde::Deserialize;
use xlc_format::library::Library;
use xlc_format::msgpack::{from_value, read_value};
use xlc_lattice::balance::by_material;
use xlc_lattice::cp::Cylinder;
use xlc_lattice::input::{CrossSections, LatticeInput};
use xlc_lattice::material::Material;
use xlc_lattice::moc::{Options, available_order, solve};
use xlc_lattice::shielding::{FuelShielding, PinShielding, ShieldingOptions, Sph, shield_fuel};
use xlc_lattice::sph::{multigroup_flux, sph_factors};
use xlc_lattice::tracking::{Boundaries, Tracking};

#[derive(Deserialize)]
struct Rates {
    flux: Vec<f64>,
    transfer: Vec<Vec<f64>>,
}

#[derive(Deserialize)]
struct Reference {
    k: f64,
    materials: BTreeMap<String, Rates>,
}

/// 物質全体の移行の断面積 [入射][出射] (1/cm)
fn transfer_xs(rates: &Rates) -> anyhow::Result<Vec<Vec<f64>>> {
    rates
        .transfer
        .iter()
        .zip(&rates.flux)
        .map(|(row, &phi)| {
            if phi > 0.0 {
                Ok(row.iter().map(|t| t / phi).collect())
            } else {
                Err(anyhow!("no flux"))
            }
        })
        .collect()
}

struct Solved {
    k: f64,
    absorption: f64,
    /// 置き換える物質の移行の断面積（物理的な値）
    swapped: Vec<Vec<f64>>,
}

fn solve_with(
    tracking: &Tracking,
    shielded: &FuelShielding,
    names: &[String],
    order: Option<usize>,
    swapped: &str,
) -> anyhow::Result<Solved> {
    let materials: Vec<&Material> = shielded
        .region_material
        .iter()
        .map(|&i| &shielded.materials[i])
        .collect();
    let options = Options {
        scattering_order: order.unwrap_or_else(|| available_order(&materials)),
        ..Options::default()
    };
    let solution = solve(tracking, &materials, &options)?;
    let physical = shielded.physical_flux(&tracking.lattice, &solution.flux)?;
    let rates = by_material(&solution, &materials, names, &physical)?;
    let rates_of = rates
        .get(swapped)
        .ok_or_else(|| anyhow!("no `{swapped}`"))?;
    Ok(Solved {
        k: solution.k,
        absorption: rates.values().flat_map(|r| &r.absorption).sum(),
        swapped: transfer_xs(&Rates {
            flux: rates_of.flux.clone(),
            transfer: rates_of.transfer.clone(),
        })?,
    })
}

/// 物理的な物質 `m` の行 `rows` に比 `ratio` を掛ける（xlc が 0 の要素は `reference` の値を置く）。全断面積は行の和の
/// 変化だけ動かす
fn swap(
    m: &Material,
    rows: &std::ops::Range<usize>,
    ratio: &[Vec<Option<f64>>],
    reference: &[Vec<f64>],
) -> Material {
    let mut m = m.clone();
    for g in rows.clone() {
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

/// [ng] の因子（`groups` の外は 1）
fn spread(mu: &[f64], groups: &std::ops::Range<usize>, ng: usize) -> Vec<f64> {
    let mut all = vec![1.0; ng];
    all[groups.clone()].copy_from_slice(mu);
    all
}

/// 1 本のピン（対称操作で移り合う組の 1 つ）の、等価なセルのリングの物質と SPH 因子
struct PinRings<'a> {
    pin: &'a PinShielding,
    sph: &'a Sph,
    /// [リング] の物質の番号
    material: Vec<usize>,
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let [_, input, reference, name, first, last] = args.as_slice() else {
        return Err(anyhow!(
            "usage: swap_transfer <lattice input> <OpenMC balance> <material> <first group> \
             <last group>"
        ));
    };
    let (first, last): (usize, usize) = (first.parse()?, last.parse()?);
    if first == 0 || last < first {
        return Err(anyhow!("groups {first}..={last}: count from 1"));
    }
    let rows = first - 1..last;
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
    let openmc = transfer_xs(
        reference
            .materials
            .get(name)
            .ok_or_else(|| anyhow!("no OpenMC `{name}`"))?,
    )?;
    if openmc.len() != ng || openmc.iter().any(|r| r.len() != ng) || rows.end > ng {
        return Err(anyhow!("the OpenMC matrix does not have {ng} groups"));
    }
    let (lattice, names) = input.geometry()?;
    let tracking =
        Tracking::with_boundaries(&lattice, input.azimuthal, input.spacing, input.boundaries)?;
    if tracking.boundaries() != Boundaries::REFLECTIVE {
        return Err(anyhow!(
            "1/k is the absorption per production only with reflective sides"
        ));
    }
    let shielded = shield_fuel(
        &library,
        &tracking,
        &input.materials,
        &names,
        &ShieldingOptions::default(),
    )?;
    if !names.iter().any(|n| n == name) {
        return Err(anyhow!("`{name}` is not in the lattice"));
    }
    // 置き換える物質の番号（`name` の領域の物質。燃料のピンの等価なセルのリングでは SPH 因子を掛けた物質）
    let mut chosen = vec![false; shielded.materials.len()];
    for (region, &i) in shielded.region_material.iter().enumerate() {
        if &names[region] == name {
            chosen[i] = true;
        }
    }
    for (region, &i) in shielded.region_material.iter().enumerate() {
        if chosen[i] && &names[region] != name {
            return Err(anyhow!(
                "material {i} is shared by `{name}` and `{}`",
                names[region]
            ));
        }
    }
    let mut pins = Vec::new();
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
            .ring_regions(&lattice)
            .ok_or_else(|| anyhow!("cell {}: no rings", pin.cell))?;
        let material: Vec<usize> = rings
            .iter()
            .map(|r| shielded.region_material[r[0]])
            .collect();
        // リングの領域がすべて同じ物質（リングの因子を掛けたもの）であること
        for (ring, &i) in rings.iter().zip(&material) {
            if ring.iter().any(|&r| shielded.region_material[r] != i) {
                return Err(anyhow!("cell {}: a ring has several materials", pin.cell));
            }
        }
        pins.push(PinRings { pin, sph, material });
    }
    // 物質ごとの SPH 因子（[ライブラリの群]。因子のない物質は持たない）。移り合うピンは物質と因子を共有する
    let mut factor_of: Vec<Option<Vec<f64>>> = vec![None; shielded.materials.len()];
    for p in &pins {
        for (&i, mu) in p.material.iter().zip(&p.sph.factors) {
            let all = spread(mu, &p.sph.groups, ng);
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
    let order = input.scattering_order;

    let base_solved = solve_with(&tracking, &shielded, &names, order, name)?;
    // 比 R = OpenMC / xlc（xlc が 0 なら None）
    let ratio: Vec<Vec<Option<f64>>> = base_solved
        .swapped
        .iter()
        .zip(&openmc)
        .map(|(x, o)| {
            x.iter()
                .zip(o)
                .map(|(&x, &o)| (x > 0.0).then(|| o / x))
                .collect()
        })
        .collect();
    // 物質ごとの物理的な物質（SPH 因子を外したもの）と、置き換えたもの
    let physical: Vec<Material> = shielded
        .materials
        .iter()
        .zip(&factor_of)
        .map(|(m, mu)| match mu {
            Some(mu) => m.scaled(&mu.iter().map(|m| 1.0 / m).collect::<Vec<_>>()),
            None => Ok(m.clone()),
        })
        .collect::<Result<_, _>>()?;
    let swapped: Vec<Material> = physical
        .iter()
        .zip(&chosen)
        .map(|(m, &chosen)| {
            if chosen {
                swap(m, &rows, &ratio, &openmc)
            } else {
                m.clone()
            }
        })
        .collect();

    // (A) 元の SPH 因子をそのまま掛ける
    let mut kept_shielding = shielded.clone();
    for (i, (m, mu)) in swapped.iter().zip(&factor_of).enumerate() {
        kept_shielding.materials[i] = match mu {
            Some(mu) => m.scaled(mu)?,
            None => m.clone(),
        };
    }

    // (B) ピンごとに、等価なセル（減速材のリング（最後のリング）の外径は UFG の計算の等価な半径）で因子を求め直す
    let bounds = library.groups.bounds();
    let mut refit_shielding = kept_shielding.clone();
    let mut refit_factor: Vec<Option<Vec<f64>>> = vec![None; shielded.materials.len()];
    let (mut deviation, mut change) = (0.0_f64, 0.0_f64);
    let mut factors_of_cell = BTreeMap::new();
    for p in &pins {
        let geometry = lattice
            .cell_type(p.pin.cell)
            .pin()
            .ok_or_else(|| anyhow!("cell {}: not a pin", p.pin.cell))?;
        let start = p.material.len() - 1;
        let cylinder = Cylinder::new(
            geometry.radii()[..start]
                .iter()
                .copied()
                .chain([p
                    .pin
                    .equivalent
                    .ok_or_else(|| anyhow!("cell {}: no equivalent radius", p.pin.cell))?])
                .collect(),
        )?;
        let groups = p.sph.groups.clone();
        let above: Vec<f64> = (0..groups.start)
            .map(|g| (bounds[g] / bounds[g + 1]).ln())
            .collect();
        let stored: Vec<Material> = p
            .material
            .iter()
            .map(|&i| shielded.materials[i].clone())
            .collect();
        let flux = multigroup_flux(&cylinder, &stored, groups.clone(), &above)?;
        let target: Vec<Vec<f64>> = flux
            .iter()
            .zip(&p.sph.factors)
            .map(|(f, mu)| f.iter().zip(mu).map(|(f, m)| f * m).collect())
            .collect();
        let fixed = vec![false; p.material.len()];
        let ring_physical: Vec<Material> =
            p.material.iter().map(|&i| physical[i].clone()).collect();
        let again = sph_factors(
            &cylinder,
            &ring_physical,
            &fixed,
            groups.clone(),
            &above,
            &target,
        )?;
        deviation = again
            .iter()
            .zip(&p.sph.factors)
            .flat_map(|(a, b)| a.iter().zip(b))
            .map(|(a, b)| (a - b).abs())
            .fold(deviation, f64::max);
        let ring_swapped: Vec<Material> = p.material.iter().map(|&i| swapped[i].clone()).collect();
        let recomputed = sph_factors(
            &cylinder,
            &ring_swapped,
            &fixed,
            groups.clone(),
            &above,
            &target,
        )?;
        change = recomputed
            .iter()
            .zip(&p.sph.factors)
            .flat_map(|(a, b)| a.iter().zip(b))
            .map(|(a, b)| (a / b - 1.0).abs())
            .fold(change, f64::max);
        for (&i, mu) in p.material.iter().zip(&recomputed) {
            let all = spread(mu, &groups, ng);
            match &refit_factor[i] {
                Some(known) if known != &all => {
                    return Err(anyhow!(
                        "material {i}: the recomputed SPH factors differ between pins"
                    ));
                }
                Some(_) => {}
                None => {
                    refit_shielding.materials[i] = swapped[i].scaled(&all)?;
                    refit_factor[i] = Some(all);
                }
            }
        }
        factors_of_cell.insert(p.pin.cell, recomputed);
    }
    for pin in &mut refit_shielding.pins {
        if let (Some(f), Some(recomputed)) = (&mut pin.sph, factors_of_cell.get(&pin.cell)) {
            f.factors = recomputed.clone();
        }
    }

    let kept = solve_with(&tracking, &kept_shielding, &names, order, name)?;
    let refitted = solve_with(&tracking, &refit_shielding, &names, order, name)?;
    println!(
        "groups {first}..={last}; {} pins; SPH factors recomputed from the reference: max |Δμ| \
         {deviation:.1e}; after the swap: max |μ'/μ − 1| {change:.2e}",
        pins.len()
    );
    println!("OpenMC k {:.6}", reference.k);
    for (label, s) in [
        ("xlc", &base_solved),
        ("swapped, SPH kept", &kept),
        ("swapped, SPH recomputed", &refitted),
    ] {
        // 置き換えた行の、群から出る移行の OpenMC との比（物質全体の平均）
        let out = |m: &[Vec<f64>], g: usize| m[g].iter().sum::<f64>() - m[g][g];
        let worst = rows
            .clone()
            .map(|g| (out(&s.swapped, g) / out(&openmc, g) - 1.0).abs())
            .fold(0.0, f64::max);
        println!(
            "{label:24} k {:.6} (1/A {:.6}), k − OpenMC {:+7.1} pcm; \
             max |out-scatter xs / OpenMC − 1| {worst:.2e}",
            s.k,
            1.0 / s.absorption,
            (s.k - reference.k) * 1e5
        );
    }
    Ok(())
}
