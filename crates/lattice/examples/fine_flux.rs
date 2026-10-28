//! 格子の入力の最初の燃料のピンの等価なセル（自己遮蔽と同じ外半径と同じリング）を UFG の減速計算で解き、燃料と減速材の格子点ごとの
//! 中性子束を書く（群の中の重みを、参照解の細かいエネルギーの中性子束と比べる。#0133）。
//!
//! 使い方: fine_flux <格子の入力 (TOML)> <出力 (msgpack)>
//!
//! 等価なセルは、ピンのリングを減速材の始まりまでそのまま使い、外半径を自己遮蔽で求めた等価なセルの外半径にする（自己遮蔽の
//! `solve_equivalent_cell` と同じ形。支持格子のあるセルは受け付けない）。範囲、飛程の求め方、格子の間引き、散乱核の標的の扱いも自己遮蔽と同じ。
//! 燃料のリングの UFG の範囲の全断面積が自己遮蔽の結果と合わなければ止まる。出力は格子点のエネルギー（降順）と、燃料のリングと
//! 減速材のリングの体積で平均したレサジーあたりの中性子束（規格化は範囲の上端で 1/E の中性子束が 1）。

use std::path::Path;

use anyhow::{Context, anyhow};
use serde::Serialize;
use xlc_format::library::Library;
use xlc_format::msgpack::to_named;
use xlc_lattice::cp::Cylinder;
use xlc_lattice::input::{CrossSections, LatticeInput, MaterialInput, VOID};
use xlc_lattice::shielding::shield_fuel;
use xlc_lattice::tracking::Tracking;
use xlc_lattice::ufg::{cell_range, held_cells, hold_total, slowing_down_fine};

#[derive(Serialize)]
struct Fine {
    format: &'static str,
    range: [f64; 2],
    radius: f64,
    energies: Vec<f64>,
    fuel: Vec<f64>,
    moderator: Vec<f64>,
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let [_, input, output] = args.as_slice() else {
        return Err(anyhow!("usage: fine_flux <lattice input> <output>"));
    };
    let path = Path::new(input);
    let input = LatticeInput::parse(&std::fs::read_to_string(path)?).context("input")?;
    let CrossSections::Library(file) = &input.cross_sections else {
        return Err(anyhow!("the input must use a library"));
    };
    let base = path.parent().unwrap_or(Path::new("."));
    let library = Library::from_msgpack(&std::fs::read(base.join(file))?)?;
    let (lattice, names) = input.geometry()?;
    let tracking =
        Tracking::with_boundaries(&lattice, input.azimuthal, input.spacing, input.boundaries)?;
    let options = input.shielding_options();
    let shielded = shield_fuel(&library, &tracking, &input.materials, &names, &options)?;
    let pin = shielded
        .pins
        .first()
        .ok_or_else(|| anyhow!("no fuel pin"))?;
    let radius = pin
        .equivalent
        .ok_or_else(|| anyhow!("the fuel has no resolved range"))?;
    let geometry = lattice
        .cell_type(pin.cell)
        .pin()
        .ok_or_else(|| anyhow!("the fuel cell is not a pin"))?;
    if geometry.frame().is_some() {
        return Err(anyhow!("cells with a grid frame are not supported"));
    }
    let (radii, sectors) = (geometry.radii(), geometry.sectors());
    let first = lattice.cell_regions(pin.cell).start;
    // リングごとの物質の名前（最後はいちばん外の円の外）
    let ring_names: Vec<&str> = (0..=radii.len())
        .map(|r| names[first + r * sectors].as_str())
        .collect();
    let material = |name: &str| -> anyhow::Result<Option<&MaterialInput>> {
        if name == VOID {
            return Ok(None);
        }
        input
            .materials
            .iter()
            .find(|m| m.name == name)
            .map(Some)
            .ok_or_else(|| anyhow!("material `{name}` is not defined"))
    };
    let rings = pin.rings.len();
    let moderator = ring_names[radii.len()];
    let start = (rings..=radii.len())
        .find(|&r| ring_names[r..].iter().all(|&n| n == moderator))
        .unwrap_or(radii.len());
    let inputs: Vec<Option<&MaterialInput>> = ring_names[..start]
        .iter()
        .chain([&moderator])
        .map(|&n| material(n))
        .collect::<Result<_, _>>()?;
    let range = cell_range(&library, &inputs)?.ok_or_else(|| anyhow!("no UFG range"))?;
    let cylinder = Cylinder::new(radii[..start].iter().copied().chain([radius]).collect())?;
    let (spectrum, fine) = slowing_down_fine(
        &library,
        &cylinder,
        &inputs,
        range,
        options.collision,
        options.thinning,
        options.scattering,
        held_cells(hold_total(options.hold_memory)?, 1),
    )?;
    // 等価なセルの形を自己遮蔽の `solve_equivalent_cell` からコピーしているので、燃料のリングの UFG の範囲の群定数が自己遮蔽の結果と
    // 同じことを確かめる（自己遮蔽の側の作り方が変わって、ここが別のセルを解くようになれば止まる）
    for (ring, shielded) in pin.rings.iter().enumerate() {
        if shielded.len() != spectrum.nuclides[ring].len() {
            return Err(anyhow!(
                "ring {ring}: the shielding has {} nuclides, this cell {}",
                shielded.len(),
                spectrum.nuclides[ring].len()
            ));
        }
        for (k, (s, f)) in shielded.iter().zip(&spectrum.nuclides[ring]).enumerate() {
            for g in spectrum.groups.clone() {
                let (a, b) = (s.total[g], f.total[g - spectrum.groups.start]);
                if (a - b).abs() > 1e-12 * a.abs().max(b.abs()) {
                    return Err(anyhow!(
                        "ring {ring} nuclide {k} group {}: the shielding gave {a}, this cell {b}",
                        g + 1
                    ));
                }
            }
        }
    }
    let volumes = cylinder.volumes();
    let average = |ring_set: std::ops::Range<usize>| -> Vec<f64> {
        let v: f64 = volumes[ring_set.clone()].iter().sum();
        (0..fine.energies.len())
            .map(|m| {
                ring_set
                    .clone()
                    .map(|i| fine.flux[i][m] * volumes[i])
                    .sum::<f64>()
                    / v
            })
            .collect()
    };
    let out = Fine {
        format: "xlc-fine-flux",
        range,
        radius,
        fuel: average(0..rings),
        moderator: average(start..start + 1),
        energies: fine.energies,
    };
    std::fs::write(output, to_named(&out)?)?;
    println!(
        "range {range:?} eV, radius {radius} cm, {} points",
        out.energies.len()
    );
    Ok(())
}
