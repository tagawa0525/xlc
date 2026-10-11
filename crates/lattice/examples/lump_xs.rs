//! 格子の入力の、構造材の塊（[`xlc_lattice::shielding`] の `shield_lumps`）の自己遮蔽した核種ごとの実効断面積を書く
//! （OpenMC との比較用。examples/vera/openmc_balance.py の `nuclides` と比べる）。
//!
//! 使い方: lump_xs <格子の入力 (TOML)> <物質の名前> <出力 (msgpack)>
//!
//! 出力は、その物質の塊ごとに、塊のあるセル、Dancoff 係数、平均弦長と、核種ごとの群ごとの全断面積・弾性散乱・正味の
//! 吸収 (barn) を、高いエネルギーの群から並べる。正味の吸収は、全断面積から、放出の数を含む散乱行列の P0 の行の和を
//! 引いたもの（`balance` の例の正味の吸収と同じ量）。核種の名前は OpenMC と同じ（Ag-107 は Ag107）にする。

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, anyhow};
use serde::Serialize;
use xlc_format::library::Library;
use xlc_format::msgpack::to_named;
use xlc_lattice::input::{CrossSections, LatticeInput};
use xlc_lattice::shielding::{ShieldingOptions, shield_fuel, shielded_nuclide};
use xlc_lattice::tracking::Tracking;

#[derive(Serialize)]
struct Rates {
    total: Vec<f64>,
    elastic: Vec<f64>,
    absorption: Vec<f64>,
}

#[derive(Serialize)]
struct Lump {
    cell: usize,
    regions: usize,
    dancoff: Vec<f64>,
    chord: f64,
    nuclides: BTreeMap<String, Rates>,
}

#[derive(Serialize)]
struct LumpXs {
    format: &'static str,
    lumps: Vec<Lump>,
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let [_, input, name, output] = args.as_slice() else {
        return Err(anyhow!(
            "usage: lump_xs <lattice input> <material> <output>"
        ));
    };
    let path = Path::new(input);
    let input = LatticeInput::parse(&std::fs::read_to_string(path)?).context("input")?;
    let CrossSections::Library(file) = &input.cross_sections else {
        return Err(anyhow!("the input must use a library"));
    };
    let base = path.parent().unwrap_or(Path::new("."));
    let library = Library::from_msgpack(&std::fs::read(base.join(file))?)?;
    let (lattice, regions) = input.geometry()?;
    let tracking =
        Tracking::with_boundaries(&lattice, input.azimuthal, input.spacing, input.boundaries)?;
    let shielded = shield_fuel(
        &library,
        &tracking,
        &input.materials,
        &regions,
        &ShieldingOptions::default(),
    )?;
    let material = input
        .materials
        .iter()
        .find(|m| &m.name == name)
        .ok_or_else(|| anyhow!("material `{name}` is not defined"))?;
    let mut lumps = Vec::new();
    for lump in shielded.lumps.iter().filter(|l| &l.material == name) {
        let mut nuclides = BTreeMap::new();
        for ((nuclide, _), effective) in material.composition.iter().zip(&lump.effective) {
            let data = shielded_nuclide(&library, nuclide, material.temperature, effective)?;
            let absorption = effective
                .total
                .iter()
                .zip(&data.scatter)
                .map(|(t, row)| t - row.iter().sum::<f64>())
                .collect();
            nuclides.insert(
                nuclide.replace('-', ""),
                Rates {
                    total: effective.total.clone(),
                    elastic: effective.elastic.clone(),
                    absorption,
                },
            );
        }
        lumps.push(Lump {
            cell: lump.cell,
            regions: lump.regions.len(),
            dancoff: lump.dancoff.clone(),
            chord: lump.chord,
            nuclides,
        });
    }
    if lumps.is_empty() {
        return Err(anyhow!("material `{name}` has no shielded lump"));
    }
    for lump in &lumps {
        println!(
            "cell {} regions {} chord {:.5} cm, Dancoff factors {:.4}..{:.4}",
            lump.cell,
            lump.regions,
            lump.chord,
            lump.dancoff.iter().copied().fold(f64::INFINITY, f64::min),
            lump.dancoff
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max),
        );
    }
    std::fs::write(
        output,
        to_named(&LumpXs {
            format: "xlc-lump-xs",
            lumps,
        })?,
    )?;
    Ok(())
}
