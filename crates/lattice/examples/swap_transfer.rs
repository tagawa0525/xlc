//! ピンセルの燃料の散乱行列のうち、指定した入射の群の行を OpenMC の行列に置き換えて解き、k の変化を測る（OpenMC との
//! k の差のうち、群の間の散乱による移行の差による分を見積もる。examples/vera/verification.md「k の差の内訳」）。
//!
//! 使い方: swap_transfer <格子の入力 (TOML)> <OpenMC の balance.msgpack> <最初の群> <最後の群>
//!
//! 群は 1 始まりで、高いエネルギーから数える。OpenMC の行列は `openmc_balance.py` の燃料の移行の率 / 中性子束（燃料の
//! 物質全体の平均）。燃料のリングごとに物理的な（SPH 因子を掛けていない）行列があるので、置き換えは要素ごとの比
//! R = OpenMC / xlc（xlc の側も燃料全体を物理的な中性子束で平均した値）を各リングの行に掛ける形にし、リングの間の違いを
//! 残す。xlc が 0 で OpenMC が正の要素（自由ガスの上方散乱）は、OpenMC の値をどのリングにも置く。P1 以上の行列は同じ比を
//! 掛ける（移行ごとの平均の余弦を保つ）。全断面積は行の和の変化だけ動かし、正味の吸収と生成は変えない。
//!
//! SPH 因子は 2 通り: (A) xlc 自身の行列で求めた因子をそのまま使う、(B) 置き換えた行列で等価なセルの因子を求め直す。
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
use xlc_lattice::shielding::{FuelShielding, ShieldingOptions, shield_fuel};
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

/// 燃料の物質全体の移行の断面積 [入射][出射] (1/cm)
fn transfer_xs(rates: &Rates) -> anyhow::Result<Vec<Vec<f64>>> {
    rates
        .transfer
        .iter()
        .zip(&rates.flux)
        .map(|(row, &phi)| {
            if phi > 0.0 {
                Ok(row.iter().map(|t| t / phi).collect())
            } else {
                Err(anyhow!("no fuel flux"))
            }
        })
        .collect()
}

struct Solved {
    k: f64,
    absorption: f64,
    /// 燃料の移行の断面積（物理的な値）
    fuel: Vec<Vec<f64>>,
}

fn solve_with(
    tracking: &Tracking,
    shielded: &FuelShielding,
    names: &[String],
    order: Option<usize>,
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
    let fuel = rates.get("fuel").ok_or_else(|| anyhow!("no `fuel`"))?;
    Ok(Solved {
        k: solution.k,
        absorption: rates.values().flat_map(|r| &r.absorption).sum(),
        fuel: transfer_xs(&Rates {
            flux: fuel.flux.clone(),
            transfer: fuel.transfer.clone(),
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

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let [_, input, reference, first, last] = args.as_slice() else {
        return Err(anyhow!(
            "usage: swap_transfer <lattice input> <OpenMC balance> <first group> <last group>"
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
            .get("fuel")
            .ok_or_else(|| anyhow!("no OpenMC `fuel`"))?,
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
    let [pin] = shielded.pins.as_slice() else {
        return Err(anyhow!("only a pin cell with one fuel pin"));
    };
    let sph = pin.sph.as_ref().ok_or_else(|| anyhow!("no SPH factors"))?;
    let geometry = lattice
        .cell_type(pin.cell)
        .pin()
        .ok_or_else(|| anyhow!("the fuel cell is not a pin"))?;
    if geometry.frame().is_some() {
        return Err(anyhow!("cells with a frame are not handled"));
    }
    let rings = pin
        .ring_regions(&lattice)
        .ok_or_else(|| anyhow!("no rings"))?;
    let ring_material: Vec<usize> = rings
        .iter()
        .map(|r| shielded.region_material[r[0]])
        .collect();
    let fuel_rings: Vec<bool> = rings.iter().map(|r| names[r[0]] == "fuel").collect();
    let order = input.scattering_order;

    let base_solved = solve_with(&tracking, &shielded, &names, order)?;
    // 比 R = OpenMC / xlc（xlc が 0 なら None）
    let ratio: Vec<Vec<Option<f64>>> = base_solved
        .fuel
        .iter()
        .zip(&openmc)
        .map(|(x, o)| {
            x.iter()
                .zip(o)
                .map(|(&x, &o)| (x > 0.0).then(|| o / x))
                .collect()
        })
        .collect();
    // リングの物理的な物質（元の SPH 因子を外す）
    let physical: Vec<Material> = ring_material
        .iter()
        .zip(&sph.factors)
        .map(|(&i, mu)| {
            let inverse: Vec<f64> = mu.iter().map(|m| 1.0 / m).collect();
            shielded.materials[i].scaled(&spread(&inverse, &sph.groups, ng))
        })
        .collect::<Result<_, _>>()?;
    let swapped: Vec<Material> = physical
        .iter()
        .zip(&fuel_rings)
        .map(|(m, &fuel)| {
            if fuel {
                swap(m, &rows, &ratio, &openmc)
            } else {
                m.clone()
            }
        })
        .collect();

    let with_factors = |factors: &[Vec<f64>]| -> anyhow::Result<FuelShielding> {
        let mut s = shielded.clone();
        for ((&i, m), mu) in ring_material.iter().zip(&swapped).zip(factors) {
            s.materials[i] = m.scaled(&spread(mu, &sph.groups, ng))?;
        }
        for p in &mut s.pins {
            if let Some(f) = &mut p.sph {
                f.factors = factors.to_vec();
            }
        }
        Ok(s)
    };

    // (B) 等価なセル: 減速材のリング（最後のリング）の外径は UFG の計算の等価な半径
    let start = rings.len() - 1;
    let cylinder = Cylinder::new(
        geometry.radii()[..start]
            .iter()
            .copied()
            .chain([pin
                .equivalent
                .ok_or_else(|| anyhow!("no equivalent radius"))?])
            .collect(),
    )?;
    let bounds = library.groups.bounds();
    let above: Vec<f64> = (0..sph.groups.start)
        .map(|g| (bounds[g] / bounds[g + 1]).ln())
        .collect();
    let stored: Vec<Material> = ring_material
        .iter()
        .map(|&i| shielded.materials[i].clone())
        .collect();
    let flux = multigroup_flux(&cylinder, &stored, sph.groups.clone(), &above)?;
    let target: Vec<Vec<f64>> = flux
        .iter()
        .zip(&sph.factors)
        .map(|(f, mu)| f.iter().zip(mu).map(|(f, m)| f * m).collect())
        .collect();
    let fixed = vec![false; rings.len()];
    let again = sph_factors(
        &cylinder,
        &physical,
        &fixed,
        sph.groups.clone(),
        &above,
        &target,
    )?;
    let deviation = again
        .iter()
        .zip(&sph.factors)
        .flat_map(|(a, b)| a.iter().zip(b))
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    let recomputed = sph_factors(
        &cylinder,
        &swapped,
        &fixed,
        sph.groups.clone(),
        &above,
        &target,
    )?;
    let change = recomputed
        .iter()
        .zip(&sph.factors)
        .flat_map(|(a, b)| a.iter().zip(b))
        .map(|(a, b)| (a / b - 1.0).abs())
        .fold(0.0, f64::max);

    let kept = solve_with(&tracking, &with_factors(&sph.factors)?, &names, order)?;
    let refitted = solve_with(&tracking, &with_factors(&recomputed)?, &names, order)?;
    println!(
        "groups {first}..={last}; SPH factors recomputed from the reference: max |Δμ| {deviation:.1e}; \
         after the swap: max |μ'/μ − 1| {change:.2e}"
    );
    println!("OpenMC k {:.6}", reference.k);
    for (label, s) in [
        ("xlc", &base_solved),
        ("swapped, SPH kept", &kept),
        ("swapped, SPH recomputed", &refitted),
    ] {
        // 置き換えた行の燃料の群から出る移行の、OpenMC との比（物質全体の平均）
        let out = |m: &[Vec<f64>], g: usize| m[g].iter().sum::<f64>() - m[g][g];
        let worst = rows
            .clone()
            .map(|g| (out(&s.fuel, g) / out(&openmc, g) - 1.0).abs())
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
