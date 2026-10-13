//! 格子の入力を自己遮蔽して MOC で解き、物質ごと・群ごとの物理的な中性子束、正味の吸収、生成と、群の間の散乱による移行を
//! 書く（OpenMC との k の差の内訳用。[`xlc_lattice::balance`]）。
//!
//! 使い方: balance <格子の入力 (TOML)> <出力 (msgpack)>
//!
//! 外周が鏡面反射の格子では、生成 1 あたりの正味の吸収の和が 1/k になる。中性子束は SPH 因子を戻した物理的な値で、
//! 反応率 / 中性子束が SPH 因子を掛ける前の物理的な実効断面積になる。群は高いエネルギーから並べる。対称性は使わない
//! （領域ごとの値をそのまま足すため）。

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, anyhow};
use serde::Serialize;
use xlc_format::library::Library;
use xlc_format::msgpack::to_named;
use xlc_lattice::balance::by_material;
use xlc_lattice::input::{CrossSections, LatticeInput};
use xlc_lattice::material::Material;
use xlc_lattice::moc::{Options, available_order, solve};
use xlc_lattice::shielding::shield_fuel;
use xlc_lattice::tracking::{Boundaries, Tracking};

#[derive(Serialize)]
struct Rates {
    flux: Vec<f64>,
    absorption: Vec<f64>,
    production: Vec<f64>,
    transfer: Vec<Vec<f64>>,
}

#[derive(Serialize)]
struct Balance {
    format: &'static str,
    k: f64,
    materials: BTreeMap<String, Rates>,
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let [_, input, output] = args.as_slice() else {
        return Err(anyhow!("usage: balance <lattice input> <output>"));
    };
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
        &input.shielding_options(),
    )?;
    let materials: Vec<&Material> = shielded
        .region_material
        .iter()
        .map(|&i| &shielded.materials[i])
        .collect();
    let options = Options {
        scattering_order: input
            .scattering_order
            .unwrap_or_else(|| available_order(&materials)),
        ..Options::default()
    };
    let solution = solve(&tracking, &materials, &options)?;
    let physical = shielded.physical_flux(&tracking.lattice, &solution.flux)?;
    let rates = by_material(&solution, &materials, &names, &physical)?;
    let absorption: f64 = rates.values().flat_map(|r| &r.absorption).sum();
    println!("k {:.6} 1/absorption {:.6}", solution.k, 1.0 / absorption);
    let materials = rates
        .into_iter()
        .map(|(name, r)| {
            (
                name,
                Rates {
                    flux: r.flux,
                    absorption: r.absorption,
                    production: r.production,
                    transfer: r.transfer,
                },
            )
        })
        .collect();
    std::fs::write(
        output,
        to_named(&Balance {
            format: "xlc-balance",
            k: solution.k,
            materials,
        })?,
    )?;
    Ok(())
}
