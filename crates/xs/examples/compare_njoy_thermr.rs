//! S(α, β) による熱中性子の散乱行列を、NJOY（THERMR, GROUPR）の GENDF の MF6 MT222 と比べる。
//!
//! 使い方: compare_njoy_thermr <ENDF> <熱中性子散乱の ENDF> <GENDF> <MAT> <温度 K> <上限 eV>
//!
//! 群構造は WIMS 69群、重み関数は iwt=4。GENDF は同じ温度・上限で THERMR（iinc=2）を通して作っておくこと。
//! 上限より下の群ごとに、行の和（非弾性散乱の群断面積）の相対差と、行の差の絶対値の和 / 行の和を表示する。
//! 環境変数 XLC_ROWS を設定すると、全群の値を表示する。

use std::collections::BTreeMap;
use std::time::Instant;

use xs::endf::Tape;
use xs::multigroup::{GroupStructure, Thermal, ThermalSab, Weight, scattering_matrix};
use xs::pointwise::{DEFAULT_TOLERANCE, reconstruct_at};
use xs::process::read_thermal_scattering;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let [_, endf, tsl, gendf, mat, temperature, cutoff] = args.as_slice() else {
        return Err("usage: compare_njoy_thermr <ENDF> <TSL> <GENDF> <MAT> <K> <cutoff eV>".into());
    };
    let (mat, temperature, cutoff): (i32, f64, f64) =
        (mat.parse()?, temperature.parse()?, cutoff.parse()?);
    let read = |path: &str| -> Result<Tape, Box<dyn std::error::Error>> {
        Ok(Tape::parse(&std::fs::read_to_string(path)?)?)
    };
    let tape = read(endf)?;
    let tables = reconstruct_at(&tape, mat, DEFAULT_TOLERANCE, temperature)?;
    let elastic = BTreeMap::from([(2, tables.get(&2).ok_or("no MT2")?.clone())]);
    let awr = xs::endf::Mf3::parse(&tape.section(mat, 3, 2).ok_or("no MF3 MT2")?)?.awr;
    let data = read_thermal_scattering(&read(tsl)?)?;
    let sab = ThermalSab::new(&data, temperature, cutoff)?;

    let groups = GroupStructure::wims69();
    let ng = groups.len();
    let start = Instant::now();
    let ours = scattering_matrix(
        &groups,
        &Weight::lwr(),
        awr,
        &elastic,
        &BTreeMap::new(),
        &BTreeMap::new(),
        Some(&Thermal::Sab(sab)),
    )?;
    println!("matrix in {:.1} s", start.elapsed().as_secs_f64());

    let reference = read(gendf)?;
    let section = reference
        .section(mat, 6, 222)
        .ok_or("no MF6 MT222 in GENDF")?;
    let njoy = read_gendf_matrix(&section, ng)?;
    let (mut worst_sum, mut worst_shape) = ((0.0_f64, 0), (0.0_f64, 0));
    for g in (0..ng).filter(|&g| groups.range(g).1 <= cutoff) {
        let (ours_sum, njoy_sum) = (ours[g].iter().sum::<f64>(), njoy[g].iter().sum::<f64>());
        let sum_diff = (ours_sum / njoy_sum - 1.0).abs();
        let shape: f64 = ours[g]
            .iter()
            .zip(&njoy[g])
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>()
            / njoy_sum;
        if std::env::var("XLC_ROWS").is_ok() {
            println!(
                "group {g:2} [{:.3e}, {:.3e}]: row sum ours {ours_sum:.5e} njoy {njoy_sum:.5e} ({sum_diff:.1e}), |diff|/sum {shape:.1e}",
                groups.range(g).0,
                groups.range(g).1
            );
        }
        if sum_diff > worst_sum.0 {
            worst_sum = (sum_diff, g);
        }
        if shape > worst_shape.0 {
            worst_shape = (shape, g);
        }
    }
    println!(
        "row sum: max rel diff {:.2e} in group {}; row shape: max |diff|/sum {:.2e} in group {}",
        worst_sum.0, worst_sum.1, worst_shape.0, worst_shape.1
    );
    Ok(())
}

/// GENDF の MF6 の行列（P0）を [入射の群][出射の群] で返す。
fn read_gendf_matrix(
    section: &xs::endf::Section<'_>,
    ng: usize,
) -> Result<Vec<Vec<f64>>, Box<dyn std::error::Error>> {
    let mut r = section.reader();
    let head = r.cont()?;
    let (nl, nz) = (usize::try_from(head.l1)?, usize::try_from(head.l2)?);
    let mut matrix = vec![vec![0.0; ng]; ng];
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
            matrix[ng - ig][ng - (ig2lo + j)] = list.values[nl * nz * (j + 1)];
        }
    }
    Ok(matrix)
}
