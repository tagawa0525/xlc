//! 無限希釈の群平均断面積を、NJOY (GROUPR) の GENDF と比べる。
//!
//! 使い方: compare_njoy_gendf <ENDF> <GENDF> <MAT> [MT...]
//!
//! 群構造は WIMS 69群、重み関数は iwt=4（軽水炉向けの値）とする。GENDF は同じ条件で
//! 0 K、無限希釈（σ0 = 1e10）として作っておくこと。MT452 は ν̄ を比べる。

use std::collections::BTreeMap;

use xs::endf::{Tape, parse_nubar};
use xs::multigroup::{GroupStructure, Weight, reaction_averages};
use xs::pointwise::{DEFAULT_TOLERANCE, reconstruct};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let [_, endf, gendf, mat, mts @ ..] = args.as_slice() else {
        return Err("usage: compare_njoy_gendf <ENDF> <GENDF> <MAT> [MT...]".into());
    };
    let mat: i32 = mat.parse()?;
    let tape = Tape::parse(&std::fs::read_to_string(endf)?)?;
    let tables = reconstruct(&tape, mat, DEFAULT_TOLERANCE)?;
    let nubar = match tape.section(mat, 1, 452) {
        Some(section) => Some(parse_nubar(&section)?),
        None => None,
    };
    let groups = GroupStructure::wims69();
    let avg = reaction_averages(&groups, &Weight::lwr(), &tables, nubar.as_ref());

    let reference = Tape::parse(&std::fs::read_to_string(gendf)?)?;
    let mts: Vec<u32> = if mts.is_empty() {
        reference.mts(mat, 3)
    } else {
        mts.iter().map(|m| m.parse()).collect::<Result<_, _>>()?
    };
    let ng = groups.len();
    for mt in mts {
        let Some(section) = reference.section(mat, 3, mt) else {
            println!("MT{mt}: not in GENDF");
            continue;
        };
        let njoy = read_gendf_mf3(&section, ng)?;
        let ours: Vec<f64> = if mt == 452 {
            let (Some(nu_fission), Some(fission)) = (&avg.nu_fission, avg.sigma.get(&18)) else {
                println!("MT452: no fission");
                continue;
            };
            nu_fission
                .iter()
                .zip(fission)
                .map(|(nf, f)| if *f > 0.0 { nf / f } else { 0.0 })
                .collect()
        } else {
            match avg.sigma.get(&mt) {
                Some(s) => s.clone(),
                None => {
                    println!("MT{mt}: not reconstructed");
                    continue;
                }
            }
        };
        let mut worst = (0.0_f64, 0);
        for (g, (&expected, &actual)) in njoy.iter().zip(&ours).enumerate() {
            if expected <= 0.0 {
                continue;
            }
            let rel = (actual - expected).abs() / expected;
            if rel > worst.0 {
                worst = (rel, g);
            }
        }
        let (lo, hi) = groups.range(worst.1);
        println!(
            "MT{mt:>3}: max rel diff {:.2e} in group {} [{lo:.4e}, {hi:.4e}] eV (ours {:.6e}, njoy {:.6e})",
            worst.0, worst.1, ours[worst.1], njoy[worst.1]
        );
    }
    Ok(())
}

/// GENDF の MF3 の群平均値（MT452 なら ν̄）を、高エネルギーの群から順に返す。
fn read_gendf_mf3(
    section: &xs::endf::Section<'_>,
    ng: usize,
) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    let mut r = section.reader();
    let head = r.cont()?;
    // 値は NL×NZ 個の中性子束のあとに断面積が並ぶ。先頭（P0、最初の σ0）を使う
    let (nl, nz) = (usize::try_from(head.l1)?, usize::try_from(head.l2)?);
    if usize::try_from(head.n2)? != ng {
        return Err(format!("GENDF has {} groups, expected {ng}", head.n2).into());
    }
    let mut values = BTreeMap::new();
    while !r.at_end() {
        let list = r.list()?;
        // GENDF の群番号 IG は低エネルギーから 1 始まり
        let ig = usize::try_from(list.n2)?;
        values.insert(ng - ig, list.values[nl * nz]);
    }
    Ok((0..ng)
        .map(|g| values.get(&g).copied().unwrap_or(0.0))
        .collect())
}
