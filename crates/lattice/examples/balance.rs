//! 格子の入力を自己遮蔽して MOC で解き、物質ごと・群ごとの物理的な中性子束、正味の吸収、生成を書く（OpenMC との k の差の
//! 内訳用）。
//!
//! 使い方: balance <格子の入力 (TOML)> <出力 (msgpack)>
//!
//! 外周が鏡面反射の格子では、生成 1 あたりの正味の吸収の和が 1/k になる。正味の吸収は、群の全断面積から、その群から
//! 出る散乱行列の行の和（(n,2n) などの放出の数を含む）を引いたもの。物質は入力の名前でまとめ、値は生成の和が 1 になる
//! ように規格化し、群は高いエネルギーから並べる。中性子束は、SPH 因子を掛けた領域と群で MOC の値に因子を掛けた物理的な
//! 値（反応率 / 中性子束が、SPH 因子を掛ける前の物理的な実効断面積になる）。対称性は使わない（領域ごとの値をそのまま
//! 足すため）。

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, anyhow};
use serde::Serialize;
use xlc_format::library::Library;
use xlc_format::msgpack::to_named;
use xlc_lattice::input::{CrossSections, LatticeInput};
use xlc_lattice::material::Material;
use xlc_lattice::moc::{Options, available_order, solve};
use xlc_lattice::shielding::{ShieldingOptions, shield_fuel};
use xlc_lattice::tracking::{Boundaries, Tracking};

#[derive(Serialize, Default)]
struct Rates {
    flux: Vec<f64>,
    absorption: Vec<f64>,
    production: Vec<f64>,
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
        &ShieldingOptions::default(),
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
    let groups = materials[0].total.len();
    // SPH 因子 μ を掛けた領域と群では、MOC の中性子束は物理的な中性子束の 1/μ（反応率 μσ × φ/μ は物理的な値のまま）。
    // 中性子束は μ を掛けて物理的な値に戻す
    let mut mu = vec![vec![1.0; groups]; names.len()];
    for pin in &shielded.pins {
        let (Some(sph), Some(rings)) = (pin.sph.as_ref(), pin.ring_regions(&tracking.lattice))
        else {
            continue;
        };
        for (factors, regions) in sph.factors.iter().zip(&rings) {
            for &region in regions {
                for (g, f) in sph.groups.clone().zip(factors) {
                    mu[region][g] = *f;
                }
            }
        }
    }
    let mut by_name: BTreeMap<String, Rates> = BTreeMap::new();
    for (region, (name, m)) in names.iter().zip(&materials).enumerate() {
        // 真空の領域には反応がない
        if m.is_void() {
            continue;
        }
        let rates = by_name.entry(name.clone()).or_insert_with(|| Rates {
            flux: vec![0.0; groups],
            absorption: vec![0.0; groups],
            production: vec![0.0; groups],
        });
        let volume = solution.volumes[region];
        for (g, factor) in mu[region].iter().enumerate() {
            let phi = solution.flux[region][g] * volume;
            let out: f64 = m.scatter[g].iter().sum();
            rates.flux[g] += phi * factor;
            rates.absorption[g] += (m.total[g] - out) * phi;
            rates.production[g] += m.nu_fission[g] * phi;
        }
    }
    let production: f64 = by_name.values().flat_map(|r| &r.production).sum();
    for rates in by_name.values_mut() {
        for v in rates
            .flux
            .iter_mut()
            .chain(&mut rates.absorption)
            .chain(&mut rates.production)
        {
            *v /= production;
        }
    }
    let absorption: f64 = by_name.values().flat_map(|r| &r.absorption).sum();
    println!("k {:.6} 1/absorption {:.6}", solution.k, 1.0 / absorption);
    std::fs::write(
        output,
        to_named(&Balance {
            format: "xlc-balance",
            k: solution.k,
            materials: by_name,
        })?,
    )?;
    Ok(())
}
