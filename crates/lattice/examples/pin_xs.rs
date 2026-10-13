//! 集合体の入力の、燃料ピンごとの自己遮蔽した核種ごとの実効断面積を書く（OpenMC との比較用）。
//!
//! 使い方: pin_xs <格子の入力 (TOML)> <出力 (msgpack)>
//!
//! 出力は examples/vera/openmc_model.py の pin-xs.msgpack と同じ形式で、核種ごとに全断面積・弾性散乱・核分裂
//! (barn) を [行（上から）][列][群（高いエネルギーから）] に並べる。燃料のないセルは 0。値は SPH 因子を掛けていない
//! 物理的な実効断面積（`PinShielding::effective`）。核種の名前は OpenMC と同じ（U-238 は U238）にする。

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, anyhow};
use serde::Serialize;
use xlc_format::library::Library;
use xlc_format::msgpack::to_named;
use xlc_lattice::input::{CrossSections, LatticeInput};
use xlc_lattice::shielding::shield_fuel;
use xlc_lattice::tracking::Tracking;

/// [行][列][群]
type Map = Vec<Vec<Vec<f64>>>;

#[derive(Serialize)]
struct Rates {
    total: Map,
    elastic: Map,
    fission: Map,
}

#[derive(Serialize)]
struct PinXs {
    format: &'static str,
    nuclides: BTreeMap<String, Rates>,
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let [_, input, output] = args.as_slice() else {
        return Err(anyhow!("usage: pin_xs <lattice input> <output>"));
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
    // ピンセルのある行と列（余白の行と列を除く）での位置
    let columns = lattice.columns();
    let (pin_rows, pin_columns) = lattice.pin_rows_and_columns();
    let n = pin_rows.len();
    if pin_columns.len() != n {
        anyhow::bail!("pins must fill a square block");
    }
    let ng = library.groups.len();
    let mut nuclides: BTreeMap<String, Rates> = BTreeMap::new();
    for pin in &shielded.pins {
        let position = |list: &[usize], k: usize| list.iter().position(|&x| x == k);
        let (Some(row), Some(column)) = (
            position(&pin_rows, pin.cell / columns),
            position(&pin_columns, pin.cell % columns),
        ) else {
            anyhow::bail!("a shielded pin outside the pin block");
        };
        for ((nuclide, _), effective) in pin.nuclides.iter().zip(&pin.effective) {
            let rates = nuclides
                .entry(nuclide.name.replace('-', ""))
                .or_insert_with(|| {
                    let zero = vec![vec![vec![0.0; ng]; n]; n];
                    Rates {
                        total: zero.clone(),
                        elastic: zero.clone(),
                        fission: zero,
                    }
                });
            rates.total[row][column].clone_from(&effective.total);
            rates.elastic[row][column].clone_from(&effective.elastic);
            rates.fission[row][column].clone_from(&effective.fission);
        }
    }
    let bytes = to_named(&PinXs {
        format: "xlc-pin-xs",
        nuclides,
    })?;
    std::fs::write(output, bytes)?;
    Ok(())
}
