//! 無限希釈の群平均断面積を、NJOY (GROUPR) の GENDF と比べる。
//!
//! 使い方: compare_njoy_gendf [--rows] <ENDF> <GENDF> <MAT> [MT...]
//!
//! 群構造は WIMS 69群、重み関数は iwt=4（軽水炉向けの値）とする。GENDF は同じ条件で
//! 0 K、無限希釈（σ0 = 1e10）として作っておくこと。MT452 は ν̄ を比べる。
//! MF5 MT18 があれば、GENDF の MF6 MT18 から求めた核分裂スペクトル χ も比べる。
//! NJOY は核分裂スペクトルを、θ や分布の表の2点目の入射エネルギー（econst）まで一定とみなす。
//! 解析式（LF=7）の U-236 では econst が 15 MeV になり、1e-5 eV の θ のスペクトルを全入射エネルギーに
//! 使うので、θ(E) の変化を積分するこちらの χ とずれる（高エネルギーの群で十数 %）。
//! GENDF の MF6（MT18 以外）は、散乱行列を Legendre の次数ごと・行ごとに比べる（行の差の絶対値の和 / P0 の行の和）。
//! 次数は GENDF にある分（GROUPR の lord）まで。
//! ENDF に MF6 がなく MF5 がある反応は、MF5 の出射エネルギー分布を使う。
//! 散乱行列はどちらも MF4 の角度分布を使う。`--rows` を付けると、全群の値を表示する。

use std::collections::BTreeMap;

use xlc_format::group::GroupStructure;
use xlc_xs::endf::{Mf3, Mf4, Mf5, Mf6, Tape, parse_nubar};
use xlc_xs::multigroup::{
    Emission, FissionNeutrons, Reactions, Weight, fission_spectrum, reaction_averages,
    scattering_moments, uses_angular,
};
use xlc_xs::pointwise::{DEFAULT_TOLERANCE, reconstruct};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (flags, args): (Vec<String>, Vec<String>) =
        std::env::args().partition(|arg| arg.starts_with("--"));
    let mut rows = false;
    for flag in &flags {
        match flag.as_str() {
            "--rows" => rows = true,
            _ => return Err(format!("unknown flag {flag}").into()),
        }
    }
    let [_, endf, gendf, mat, mts @ ..] = args.as_slice() else {
        return Err("usage: compare_njoy_gendf [--rows] <ENDF> <GENDF> <MAT> [MT...]".into());
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

    if let (Some(fission), Some(nu), Some(section), Some(gendf)) = (
        tables.get(&18),
        &nubar,
        tape.section(mat, 5, 18),
        reference.section(mat, 6, 18),
    ) {
        // GROUPR の MF6 MT18 は MF5 MT18 の分布だけを使うので、遅発中性子の分布（MF5 MT455）は混ぜない
        let ours = fission_spectrum(
            &groups,
            &Weight::lwr(),
            fission,
            &[FissionNeutrons {
                nubar: nu.clone(),
                spectrum: Mf5::parse(&section)?,
            }],
        );
        let njoy = read_gendf_chi(&gendf, ng)?;
        let mut worst = (0.0_f64, 0);
        for (g, (&expected, &actual)) in njoy.iter().zip(&ours).enumerate() {
            // 寄与の小さい群は相対差が大きくなるので、1e-6 以上の群だけを見る
            if expected < 1.0e-6 {
                continue;
            }
            let rel = (actual - expected).abs() / expected;
            if rel > worst.0 {
                worst = (rel, g);
            }
        }
        println!(
            "chi  : max rel diff {:.2e} in group {} (ours {:.6e}, njoy {:.6e}), sum of |diff| {:.2e}",
            worst.0,
            worst.1,
            ours[worst.1],
            njoy[worst.1],
            ours.iter()
                .zip(&njoy)
                .map(|(a, b)| (a - b).abs())
                .sum::<f64>()
        );
    }

    let awr = Mf3::parse(&tape.section(mat, 3, 2).ok_or("no MT2")?)?.awr;
    // MT18 は χ として上で比べた
    for mt in reference.mts(mat, 6).into_iter().filter(|&mt| mt != 18) {
        let (Some(gendf), Some(table)) = (reference.section(mat, 6, mt), tables.get(&mt)) else {
            continue;
        };
        let q_values = BTreeMap::from([(
            mt,
            Mf3::parse(&tape.section(mat, 3, mt).ok_or("no MF3")?)?.qi,
        )]);
        let mut emissions = BTreeMap::new();
        if let Some(section) = tape.section(mat, 6, mt) {
            match Mf6::parse(&section) {
                Ok(data) => {
                    emissions.insert(mt, Emission::Mf6(data));
                }
                Err(e) => {
                    println!("MF6 MT{mt}: skipped ({e})");
                    continue;
                }
            }
        } else if let Some(section) = tape.section(mat, 5, mt) {
            emissions.insert(mt, Emission::Mf5(Mf5::parse(&section)?));
        }
        // 角度分布は、弾性散乱・離散準位と、MF5 で与えられた反応の分（xs の処理と同じ）
        let mut angular = BTreeMap::new();
        if uses_angular(mt, emissions.get(&mt))
            && let Some(section) = tape.section(mat, 4, mt)
        {
            angular.insert(mt, Mf4::parse(&section)?);
        }
        let single = BTreeMap::from([(mt, table.clone())]);
        let njoy = read_gendf_matrix(&gendf, ng)?;
        let order = njoy.len() - 1;
        let ours = match scattering_moments(
            &groups,
            &Weight::lwr(),
            awr,
            &Reactions {
                tables: &single,
                q_values: &q_values,
                emissions: &emissions,
                angular: &angular,
            },
            None,
            order,
        ) {
            Ok(m) => m,
            Err(e) => {
                println!("MF6 MT{mt}: skipped ({e})");
                continue;
            }
        };
        // P_l の行列は、行の差の絶対値の和を P0 の行の和で割る
        let p0_sums: Vec<f64> = ours[0].iter().map(|row| row.iter().sum()).collect();
        for (l, (ours, njoy_l)) in ours.iter().zip(&njoy).enumerate() {
            let mut worst = (0.0_f64, 0);
            for g in 0..ng {
                // どちらか一方にだけ散乱がある群も差として数える
                let sum = njoy[0][g].iter().sum::<f64>().max(p0_sums[g]);
                if sum <= 0.0 {
                    continue;
                }
                let diff: f64 = ours[g]
                    .iter()
                    .zip(&njoy_l[g])
                    .map(|(a, b)| (a - b).abs())
                    .sum();
                if diff / sum > worst.0 {
                    worst = (diff / sum, g);
                }
                if rows {
                    println!("  MT{mt} P{l} group {g:2}: {:.2e}", diff / sum);
                }
            }
            let g = worst.1;
            println!(
                "MF6 MT{mt} P{l}: max row |diff| / P0 row sum {:.2e} in group {g} (row sum ours {:.6e}, njoy {:.6e}; self ours {:.6e}, njoy {:.6e})",
                worst.0,
                ours[g].iter().sum::<f64>(),
                njoy_l[g].iter().sum::<f64>(),
                ours[g][g],
                njoy_l[g][g]
            );
        }
    }
    Ok(())
}

/// GENDF の MF6 の行列を、Legendre の次数ごとに [l][入射の群][出射の群] で返す（最初の σ0）。
/// 値は出射の群ごとに (NL, NZ) の順（NL が速い）に並ぶ。
fn read_gendf_matrix(
    section: &xlc_xs::endf::Section<'_>,
    ng: usize,
) -> Result<Vec<Vec<Vec<f64>>>, Box<dyn std::error::Error>> {
    let mut r = section.reader();
    let head = r.cont()?;
    let (nl, nz) = (usize::try_from(head.l1)?, usize::try_from(head.l2)?);
    let mut matrix = vec![vec![vec![0.0; ng]; ng]; nl];
    while !r.at_end() {
        let list = r.list()?;
        let (ng2, ig2lo, ig) = (
            usize::try_from(list.l1)?,
            usize::try_from(list.l2)?,
            usize::try_from(list.n2)?,
        );
        if ig == 0 || ig2lo == 0 {
            continue;
        }
        for j in 0..ng2.saturating_sub(1) {
            for (l, m) in matrix.iter_mut().enumerate() {
                m[ng - ig][ng - (ig2lo + j)] = list.values[nl * nz * (j + 1) + l];
            }
        }
    }
    Ok(matrix)
}

/// GENDF の MF6 MT18 から、核分裂中性子の生成を出射の群ごとに集計し、和を 1 にして返す。
///
/// 低い群は「群ごとの νσf × 共通のスペクトル（IG=0 のレコード）」、高い群は生成行列の行で
/// 与えられる。どちらも群の中性子束を掛けて足し合わせる。
fn read_gendf_chi(
    section: &xlc_xs::endf::Section<'_>,
    ng: usize,
) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    let mut r = section.reader();
    let head = r.cont()?;
    let nlnz = usize::try_from(head.l1 * head.l2)?;
    let mut spectrum = vec![0.0; ng];
    let mut spectrum_source = 0.0;
    let mut production = vec![0.0; ng];
    while !r.at_end() {
        let list = r.list()?;
        let (ng2, ig2lo, ig) = (
            usize::try_from(list.l1)?,
            usize::try_from(list.l2)?,
            usize::try_from(list.n2)?,
        );
        if ig == 0 {
            for (j, &v) in list.values.iter().enumerate() {
                spectrum[ng - (ig2lo + j)] = v;
            }
        } else if ig2lo == 0 {
            spectrum_source += list.values[0] * list.values[nlnz];
        } else {
            let flux = list.values[0];
            for j in 0..ng2.saturating_sub(1) {
                production[ng - (ig2lo + j)] += flux * list.values[nlnz * (j + 1)];
            }
        }
    }
    let chi: Vec<f64> = production
        .iter()
        .zip(&spectrum)
        .map(|(p, s)| p + spectrum_source * s)
        .collect();
    let total: f64 = chi.iter().sum();
    Ok(chi.iter().map(|c| c / total).collect())
}

/// GENDF の MF3 の群平均値（MT452 なら ν̄）を、高エネルギーの群から順に返す。
fn read_gendf_mf3(
    section: &xlc_xs::endf::Section<'_>,
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
