//! 均質な無限媒質（`examples/vera/homogeneous_medium.py` の入力）の k を、xlc のライブラリで求める（OpenMC の
//! `openmc_homogeneous.py` と比べる）。
//!
//! 使い方: homogeneous <媒質 (TOML)> <ライブラリ (msgpack)> [--balance <出力 (msgpack)>]
//!         [--swap <OpenMC の balance.msgpack> <核種>:<最初の群>-<最後の群> ...]
//!
//! 自己遮蔽は 3 通りで解く: (1) 無限希釈、(2) 背景断面積（同じ媒質の他の核種の全断面積の和 / 個数密度。逃げる
//! 項はない。分離共鳴のない核種は無限希釈）の反復、(3) (2) に、UFG の範囲の群を減速計算の値に置き換えたもの（1 つのリングの白色境界の円筒は、
//! 平坦な中性子束の無限媒質と同じ）。(3) が格子の燃料の扱いに対応する。散乱行列は、格子と同じく無限希釈の行の形
//! のまま、弾性散乱の変化だけ行の和を変える（`shielded_nuclide`）。背景断面積が表の上限より大きい群は、無限希釈の
//! 値を使う（表の上限が無限希釈）。
//!
//! k は、(Σt − Sᵀ) φ = (1/k) F φ（F は核種ごとの χ と νΣf の積の和）の最大の固有値を、べき乗法で求める。
//!
//! 釣り合いの出力は (3) の解の核種ごとの群ごとの中性子束、正味の吸収（全断面積 − 散乱行列の行の和）、生成、
//! 散乱による移行で、生成の和が 1 になるように規格化する（`balance` の例と同じ形式で、物質の代わりに核種。
//! `compare_balance.py` で OpenMC（`openmc_homogeneous.py`）と比べる）。
//!
//! `--swap` は、(3) の媒質の、指定した核種の入射の群（1 始まり、高いエネルギーから）の P0 の散乱行列の行を、OpenMC の
//! 移行の率 / 中性子束（`openmc_homogeneous.py` の核種ごとの `transfer`）に置き換えて解く。全断面積は行の和の変化だけ
//! 動かし、正味の吸収と生成は変えない（無限媒質の k は P0 の行列だけで決まる）。

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, anyhow};
use serde::{Deserialize, Serialize};
use xlc_format::library::{Library, LookupError, Shielded};
use xlc_format::msgpack::{from_value, read_value, to_named};
use xlc_lattice::cp::Cylinder;
use xlc_lattice::input::MaterialInput;
use xlc_lattice::material::Material;
use xlc_lattice::shielding::shielded_nuclide;
use xlc_lattice::ufg::{cell_range, slowing_down};
use xlc_numeric::lu::Lu;

#[derive(Deserialize)]
struct Medium {
    temperature: f64,
    nuclides: BTreeMap<String, f64>,
}

#[derive(Deserialize)]
struct Counted {
    flux: Vec<f64>,
    transfer: Vec<Vec<f64>>,
}

#[derive(Deserialize)]
struct Reference {
    materials: BTreeMap<String, Counted>,
}

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

/// 背景断面積の反復の打ち切り（全断面積の相対変化）と回数の上限
const TOLERANCE: f64 = 1.0e-10;
const ITERATIONS: usize = 1000;

fn bondarenko(library: &Library, medium: &MaterialInput) -> anyhow::Result<Vec<Shielded>> {
    let t = medium.temperature;
    let ng = library.groups.len();
    let (Some(&max), Some(&min)) = (library.sigma0.first(), library.sigma0.last()) else {
        return Err(LookupError::NoSigma0.into());
    };
    let mut tables: Vec<Shielded> = medium
        .composition
        .iter()
        .map(|(name, _)| library.shielded(name, t, &vec![max; ng]))
        .collect::<Result<_, _>>()?;
    for _ in 0..ITERATIONS {
        let next: Vec<Shielded> = medium
            .composition
            .iter()
            .enumerate()
            .map(|(k, (name, density))| {
                // 分離共鳴のない核種の表は σ0 でほとんど変わらず、薄い媒質では σ0 が表の下限を下回るので、
                // 無限希釈にする
                let resonant = library
                    .nuclides
                    .iter()
                    .find(|n| n.name == *name)
                    .ok_or_else(|| LookupError::Nuclide(name.clone()))?
                    .resolved
                    .is_some();
                if !resonant {
                    return Ok(library.shielded(name, t, &vec![max; ng])?);
                }
                let sigma0 = (0..ng)
                    .map(|g| {
                        let background: f64 = medium
                            .composition
                            .iter()
                            .zip(&tables)
                            .enumerate()
                            .filter(|&(j, _)| j != k)
                            .map(|(_, ((_, n), s))| n * s.total[g])
                            .sum();
                        let s = background / density;
                        if s < min {
                            return Err(anyhow!("{name}: σ0 {s} below the table in group {g}"));
                        }
                        Ok(s.min(max))
                    })
                    .collect::<anyhow::Result<Vec<_>>>()?;
                Ok(library.shielded(name, t, &sigma0)?)
            })
            .collect::<anyhow::Result<_>>()?;
        let change = next
            .iter()
            .zip(&tables)
            .flat_map(|(a, b)| a.total.iter().zip(&b.total))
            .map(|(a, b)| ((a - b) / b).abs())
            .fold(0.0, f64::max);
        tables = next;
        if change <= TOLERANCE {
            return Ok(tables);
        }
    }
    Err(anyhow!("the background cross sections did not converge"))
}

fn material(
    library: &Library,
    medium: &MaterialInput,
    tables: &[Shielded],
) -> anyhow::Result<Material> {
    let nuclides = medium
        .composition
        .iter()
        .zip(tables)
        .map(|((name, density), s)| {
            Ok((
                shielded_nuclide(library, name, medium.temperature, s)?,
                *density,
            ))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    Ok(Material::from_nuclides(&nuclides, library.groups.len())?)
}

/// 無限媒質の k と中性子束（べき乗法）
fn k_infinity(m: &Material) -> anyhow::Result<(f64, Vec<f64>)> {
    let ng = m.groups();
    // A[h][g] φ_g: 群 h の釣り合い（Σt,h φ_h − Σ_g Σs(g→h) φ_g）
    let mut a = vec![0.0; ng * ng];
    for h in 0..ng {
        a[h * ng + h] += m.total[h];
        for g in 0..ng {
            a[h * ng + g] -= m.scatter[g][h];
        }
    }
    let lu = Lu::new(a, ng).ok_or_else(|| anyhow!("singular balance"))?;
    let fission = m.fission_matrix();
    let mut flux = vec![1.0; ng];
    let mut k = 1.0;
    for _ in 0..100_000 {
        let mut source: Vec<f64> = (0..ng)
            .map(|h| (0..ng).map(|g| fission[g][h] * flux[g]).sum())
            .collect();
        let before: f64 = source.iter().sum();
        lu.solve(&mut source);
        let after: f64 = (0..ng)
            .map(|h| (0..ng).map(|g| fission[g][h] * source[g]).sum::<f64>())
            .sum();
        let next = after / before;
        flux = source;
        if (next - k).abs() < 1.0e-12 {
            return Ok((next, flux));
        }
        k = next;
    }
    Err(anyhow!("the power iteration did not converge"))
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let usage = || {
        anyhow!(
            "usage: homogeneous <medium> <library> [--balance <output>] \
             [--swap <OpenMC balance> <nuclide>:<first>-<last> ...]"
        )
    };
    let [_, medium_path, library_path, rest @ ..] = args.as_slice() else {
        return Err(usage());
    };
    let mut rest: &[String] = rest;
    let mut output = None;
    let mut swap = None;
    while let [flag, value, tail @ ..] = rest {
        match flag.as_str() {
            "--balance" => {
                output = Some(value);
                rest = tail;
            }
            "--swap" => {
                let tail: &[String] = tail;
                swap = Some((value, tail));
                rest = &[];
            }
            _ => return Err(usage()),
        }
    }
    if !rest.is_empty() {
        return Err(usage());
    }
    let medium: Medium =
        toml::from_str(&std::fs::read_to_string(medium_path)?).context("medium")?;
    let library = Library::from_msgpack(&std::fs::read(Path::new(library_path))?)?;
    let ng = library.groups.len();
    let input = MaterialInput {
        name: "mix".into(),
        temperature: medium.temperature,
        composition: medium.nuclides.into_iter().collect(),
    };

    let infinite = library.sigma0[0];
    let dilute: Vec<Shielded> = input
        .composition
        .iter()
        .map(|(name, _)| library.shielded(name, input.temperature, &vec![infinite; ng]))
        .collect::<Result<_, _>>()?;
    let background = bondarenko(&library, &input)?;
    let mut ufg = background.clone();
    let range = cell_range(&library, &[Some(&input)])?;
    if let Some(range) = range {
        let spectrum = slowing_down(&library, &Cylinder::new(vec![1.0])?, &[Some(&input)], range)?;
        let top = spectrum.groups.start;
        for (s, u) in ufg.iter_mut().zip(&spectrum.nuclides[0]) {
            for g in spectrum.groups.clone() {
                s.total[g] = u.total[g - top];
                s.elastic[g] = u.elastic[g - top];
                s.fission[g] = u.fission[g - top];
                s.nu_fission[g] = u.nu_fission[g - top];
                s.capture[g] = u.capture[g - top];
            }
        }
    }
    println!("{medium_path}");
    for (label, tables) in [
        ("infinite dilution", &dilute),
        ("background cross sections", &background),
        ("background + UFG", &ufg),
    ] {
        let (k, _) = k_infinity(&material(&library, &input, tables)?)?;
        println!("  {label:28} k {k:.6}");
    }
    // 核種ごとの物質
    let parts = input
        .composition
        .iter()
        .zip(&ufg)
        .map(|((name, density), s)| {
            let nuclide = shielded_nuclide(&library, name, input.temperature, s)?;
            Ok((
                name.clone(),
                Material::from_nuclides(&[(nuclide, *density)], ng)?,
            ))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    if let Some(output) = output {
        let (k, flux) = k_infinity(&material(&library, &input, &ufg)?)?;
        let produced: f64 = parts
            .iter()
            .flat_map(|(_, m)| m.nu_fission.iter().zip(&flux).map(|(n, f)| n * f))
            .sum();
        let phi: Vec<f64> = flux.iter().map(|f| f / produced).collect();
        let materials: BTreeMap<String, Rates> = parts
            .iter()
            .map(|(name, m)| {
                let rates = Rates {
                    flux: phi.clone(),
                    absorption: (0..ng)
                        .map(|g| (m.total[g] - m.scatter[g].iter().sum::<f64>()) * phi[g])
                        .collect(),
                    production: (0..ng).map(|g| m.nu_fission[g] * phi[g]).collect(),
                    transfer: (0..ng)
                        .map(|g| m.scatter[g].iter().map(|s| s * phi[g]).collect())
                        .collect(),
                };
                (name.clone(), rates)
            })
            .collect();
        let absorption: f64 = materials.values().flat_map(|r| &r.absorption).sum();
        println!("  balance: k {k:.6}, 1/absorption {:.6}", 1.0 / absorption);
        std::fs::write(
            output,
            to_named(&Balance {
                format: "xlc-balance",
                k,
                materials,
            })?,
        )?;
    }
    if let Some((reference, specs)) = swap {
        let reference: Reference = from_value(read_value(&std::fs::read(reference)?)?)?;
        let mut mix = material(&library, &input, &ufg)?;
        let (before, _) = k_infinity(&mix)?;
        for spec in specs {
            let error = || anyhow!("`{spec}`: expected <nuclide>:<first>-<last>");
            let (name, range) = spec.split_once(':').ok_or_else(error)?;
            let (first, last) = range.split_once('-').ok_or_else(error)?;
            let (first, last): (usize, usize) = (first.parse()?, last.parse()?);
            if first == 0 || last < first || last > ng {
                return Err(anyhow!("`{spec}`: groups 1..={ng}"));
            }
            let (_, mine) = parts
                .iter()
                .find(|(n, _)| n == name)
                .ok_or_else(|| anyhow!("`{name}` is not in the medium"))?;
            let theirs = reference
                .materials
                .get(name)
                .ok_or_else(|| anyhow!("`{name}` is not in the OpenMC balance"))?;
            if theirs.flux.len() != ng || theirs.transfer.len() != ng {
                return Err(anyhow!("the OpenMC balance does not have {ng} groups"));
            }
            for g in first - 1..last {
                if theirs.flux[g].is_nan() || theirs.flux[g] <= 0.0 {
                    return Err(anyhow!("no OpenMC flux in group {}", g + 1));
                }
                for h in 0..ng {
                    let change = theirs.transfer[g][h] / theirs.flux[g] - mine.scatter[g][h];
                    mix.scatter[g][h] += change;
                    mix.total[g] += change;
                }
            }
        }
        let (k, _) = k_infinity(&mix)?;
        println!(
            "  swapped {}: k {k:.6} ({:+.1} pcm)",
            specs.join(" "),
            (k - before) * 1e5
        );
    }
    if let Some(range) = range {
        println!("  UFG range {:.4e}〜{:.4e} eV", range[0], range[1]);
    }
    Ok(())
}
