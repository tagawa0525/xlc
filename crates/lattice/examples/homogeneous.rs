//! 均質な無限媒質（`examples/vera/homogeneous_medium.py` の入力）の k を、xlc のライブラリで求める（OpenMC の
//! `openmc_homogeneous.py` と比べる）。
//!
//! 使い方: homogeneous <媒質 (TOML)> <ライブラリ (msgpack)>
//!
//! 自己遮蔽は 3 通りで解く: (1) 無限希釈、(2) 背景断面積（同じ媒質の他の核種の全断面積の和 / 個数密度。逃げる
//! 項はない。分離共鳴のない核種は無限希釈）の反復、(3) (2) に、UFG の範囲の群を減速計算の値に置き換えたもの（1 つのリングの白色境界の円筒は、
//! 平坦な中性子束の無限媒質と同じ）。(3) が格子の燃料の扱いに対応する。散乱行列は、格子と同じく無限希釈の行の形
//! のまま、弾性散乱の変化だけ行の和を変える（`shielded_nuclide`）。背景断面積が表の上限より大きい群は、無限希釈の
//! 値を使う（表の上限が無限希釈）。
//!
//! k は、(Σt − Sᵀ) φ = (1/k) F φ（F は核種ごとの χ と νΣf の積の和）の最大の固有値を、べき乗法で求める。

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, anyhow};
use serde::Deserialize;
use xlc_format::library::{Library, LookupError, Shielded};
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

/// 無限媒質の k（べき乗法）
fn k_infinity(m: &Material) -> anyhow::Result<f64> {
    let ng = m.groups();
    // A[h][g] φ_g: 群 h のバランス（Σt,h φ_h − Σ_g Σs(g→h) φ_g）
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
            return Ok(next);
        }
        k = next;
    }
    Err(anyhow!("the power iteration did not converge"))
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let [_, medium_path, library_path] = args.as_slice() else {
        return Err(anyhow!("usage: homogeneous <medium> <library>"));
    };
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
            }
        }
    }
    println!("{medium_path}");
    for (label, tables) in [
        ("infinite dilution", &dilute),
        ("background cross sections", &background),
        ("background + UFG", &ufg),
    ] {
        let k = k_infinity(&material(&library, &input, tables)?)?;
        println!("  {label:28} k {k:.6}");
    }
    if let Some(range) = range {
        println!("  UFG range {:.4e}〜{:.4e} eV", range[0], range[1]);
    }
    Ok(())
}
