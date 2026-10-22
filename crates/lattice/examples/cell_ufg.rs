//! 格子の入力の 1 つのセルを、面積の等しい円筒（白色境界）にして超微細群の減速計算（[`xlc_lattice::ufg`]）で解き、
//! 指定した物質の核種ごとの実効断面積を書く（吸収体の棒の自己遮蔽を、等価原理によらずに求めて比べる）。
//!
//! 使い方: cell_ufg <格子の入力 (TOML)> <セルの名前> <物質の名前> <出力 (msgpack)> [--radius <外半径 (cm)>]
//!
//! セルの帯を入力のリングの数で等体積に分け、最後の帯はピッチの正方形と面積の等しい円（`--radius` を書けばその
//! 半径の円）まで延ばす（隣のセルの影は入れない。外周は白色境界なので、外半径が小さいほど同じ棒の格子の影が入る）。範囲は、ライブラリのポイントワイズの断面積の下限から、物質の核種の分離共鳴の上端のうち最も高いもの
//! （[`xlc_lattice::ufg::MAX_TOP`] で頭打ち）以下で最も高い群の境界まで。範囲の群は、物質のリングの値を中性子束と
//! 体積で平均し、範囲の外の群は無限希釈の値にする。出力は `lump_xs` の例と同じ形式で、塊は 1 つ（`cell` は 0、
//! Dancoff 係数は空、平均弦長は 0）。compare_lump_xs.py で OpenMC と比べる。

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, anyhow};
use serde::Serialize;
use xlc_format::library::Library;
use xlc_format::msgpack::to_named;
use xlc_lattice::cp::Cylinder;
use xlc_lattice::input::{CrossSections, LatticeInput, MaterialInput, VOID};
use xlc_lattice::ufg::{MAX_TOP, no_background, slowing_down_with};

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
    let usage =
        || anyhow!("usage: cell_ufg <lattice input> <cell> <material> <output> [--radius <cm>]");
    let (positional, radius) = match args.as_slice() {
        [rest @ .., flag, value] if flag == "--radius" => (rest, Some(value.parse::<f64>()?)),
        rest => (rest, None),
    };
    let [_, input, cell_name, name, output] = positional else {
        return Err(usage());
    };
    let path = Path::new(input);
    let input = LatticeInput::parse(&std::fs::read_to_string(path)?).context("input")?;
    let CrossSections::Library(file) = &input.cross_sections else {
        return Err(anyhow!("the input must use a library"));
    };
    let base = path.parent().unwrap_or(Path::new("."));
    let library = Library::from_msgpack(&std::fs::read(base.join(file))?)?;
    let cell = input
        .cells
        .iter()
        .find(|c| &c.name == cell_name)
        .ok_or_else(|| anyhow!("cell `{cell_name}` is not defined"))?;
    let material_of = |zone: &str| -> anyhow::Result<Option<&MaterialInput>> {
        if zone == VOID {
            return Ok(None);
        }
        input
            .materials
            .iter()
            .find(|m| m.name == zone)
            .map(Some)
            .ok_or_else(|| anyhow!("material `{zone}` is not defined"))
    };
    let edge = radius.unwrap_or(input.pitch / std::f64::consts::PI.sqrt());
    let mut radii = Vec::new();
    let mut rings = Vec::new();
    let mut inner: f64 = 0.0;
    for zone in &cell.zones {
        let outer = zone.outer.unwrap_or(edge).min(edge);
        if outer <= inner {
            continue;
        }
        for i in 1..=zone.rings {
            let fraction = i as f64 / zone.rings as f64;
            radii.push((inner * inner + fraction * (outer * outer - inner * inner)).sqrt());
            rings.push(material_of(&zone.material)?);
        }
        inner = outer;
    }
    // 外半径が帯の外半径より大きければ、最後の帯の物質で延ばす
    if inner < edge {
        let last = cell
            .zones
            .last()
            .ok_or_else(|| anyhow!("cell `{cell_name}` has no zones"))?;
        radii.push(edge);
        rings.push(material_of(&last.material)?);
    }
    let target = input
        .materials
        .iter()
        .find(|m| &m.name == name)
        .ok_or_else(|| anyhow!("material `{name}` is not defined"))?;
    let mut top: f64 = 0.0;
    for (nuclide, _) in &target.composition {
        if let Some(r) = &library.lookup(nuclide, target.temperature)?.0.resolved {
            top = top.max(r.top);
        }
    }
    let pointwise = library
        .pointwise_range
        .ok_or_else(|| anyhow!("the library has no pointwise data"))?;
    let top = top.min(MAX_TOP).min(pointwise[1]);
    let bounds = library.groups.bounds();
    let high = bounds
        .iter()
        .copied()
        .find(|&b| b <= top)
        .ok_or_else(|| anyhow!("no group bound below {top} eV"))?;
    let range = [pointwise[0], high];
    let first = bounds
        .iter()
        .position(|&b| b == high)
        .ok_or_else(|| anyhow!("no group bound at {high} eV"))?;
    let cylinder = Cylinder::new(radii.clone())?;
    let spectrum = slowing_down_with(
        &library,
        &cylinder,
        &rings,
        &no_background(&rings),
        range,
        input.collision,
    )?;
    let volumes = cylinder.volumes();
    let mine: Vec<usize> = (0..rings.len())
        .filter(|&j| rings[j].is_some_and(|m| &m.name == name))
        .collect();
    let ng = library.groups.len();
    let infinite = *library
        .sigma0
        .first()
        .ok_or_else(|| anyhow!("the library has no sigma0"))?;
    let mut nuclides = BTreeMap::new();
    for (k, (nuclide, _)) in target.composition.iter().enumerate() {
        let dilute = library.shielded(nuclide, target.temperature, &vec![infinite; ng])?;
        let data = library.nuclide(nuclide, target.temperature)?;
        let mut total = dilute.total.clone();
        let mut elastic = dilute.elastic.clone();
        for g in spectrum.groups.clone() {
            let weights: Vec<f64> = mine
                .iter()
                .map(|&j| spectrum.flux[j][g - first] * volumes[j])
                .collect();
            let sum: f64 = weights.iter().sum();
            let average = |of: fn(&xlc_format::library::Shielded) -> &Vec<f64>| -> f64 {
                mine.iter()
                    .zip(&weights)
                    .map(|(&j, w)| w * of(&spectrum.nuclides[j][k])[g - first])
                    .sum::<f64>()
                    / sum
            };
            total[g] = average(|s| &s.total);
            elastic[g] = average(|s| &s.elastic);
        }
        // 正味の吸収は、無限希釈の散乱行列の行の和を、弾性散乱の変化だけ動かして引く（`shielded_nuclide` と同じ）
        let absorption = (0..ng)
            .map(|g| {
                let row: f64 = data.scatter[g].iter().sum();
                total[g] - (row + elastic[g] - dilute.elastic[g])
            })
            .collect();
        nuclides.insert(
            nuclide.replace('-', ""),
            Rates {
                total,
                elastic,
                absorption,
            },
        );
    }
    println!(
        "range {:.4}..{:.1} eV, groups {:?}, rings {:?}",
        range[0], range[1], spectrum.groups, radii
    );
    std::fs::write(
        output,
        to_named(&LumpXs {
            format: "xlc-lump-xs",
            lumps: vec![Lump {
                cell: 0,
                regions: 1,
                dancoff: Vec::new(),
                chord: 0.0,
                nuclides,
            }],
        })?,
    )?;
    Ok(())
}
