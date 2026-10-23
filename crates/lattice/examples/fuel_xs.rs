//! 格子の入力の、最初の燃料のピンの自己遮蔽した核種ごとの実効断面積を書く（OpenMC との比較用）。
//!
//! 使い方: fuel_xs <格子の入力 (TOML)> <出力 (msgpack)>
//!
//! 出力は examples/vera/openmc_model.py の fuel-xs.msgpack と同じ形式に捕獲を足したもので、核種ごとに群ごとの
//! 全断面積・弾性散乱・核分裂・捕獲（MT102）(barn) を、高いエネルギーの群から並べる（散乱行列の行は放出される
//! 中性子の数を数えるので、吸収は行の和から求めずに、全断面積と弾性散乱で比べる）。核種の名前は OpenMC と同じ
//! （U-238 は U238）にする。

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, anyhow};
use serde::Serialize;
use xlc_format::library::Library;
use xlc_format::msgpack::to_named;
use xlc_lattice::input::{CrossSections, LatticeInput};
use xlc_lattice::shielding::shield_fuel;
use xlc_lattice::tracking::Tracking;

#[derive(Serialize)]
struct Rates {
    total: Vec<f64>,
    elastic: Vec<f64>,
    fission: Vec<f64>,
    capture: Vec<f64>,
}

#[derive(Serialize)]
struct FuelXs {
    format: &'static str,
    nuclides: BTreeMap<String, Rates>,
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let [_, input, output] = args.as_slice() else {
        return Err(anyhow!("usage: fuel_xs <lattice input> <output>"));
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
        &input.shielding_options(),
    )?;
    let pin = shielded
        .pins
        .first()
        .ok_or_else(|| anyhow!("no fuel pin"))?;
    let nuclides = pin
        .nuclides
        .iter()
        .zip(pin.effective.iter())
        .map(|((n, _), s)| {
            (
                n.name.replace('-', ""),
                Rates {
                    total: s.total.clone(),
                    elastic: s.elastic.clone(),
                    fission: s.fission.clone(),
                    capture: s.capture.clone(),
                },
            )
        })
        .collect();
    let bytes = to_named(&FuelXs {
        format: "xlc-fuel-xs",
        nuclides,
    })?;
    std::fs::write(output, bytes)?;
    println!("Dancoff factors: {:?}", pin.dancoff);
    Ok(())
}
