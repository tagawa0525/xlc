//! ドップラー拡がりを、NJOY (BROADR) の出力と比べる。
//!
//! 使い方: compare_njoy_broadr <ENDF> <0 K の PENDF> <拡げた PENDF> <MAT> <温度 K> [MT...]
//!
//! NJOY RECONR が作った 0 K のポイントワイズの断面積（PENDF）の部分反応を、こちらの方法でまとめて拡げ、
//! 同じ 0 K の断面積を BROADR で拡げた結果と比べる。0 K の格子の違いを除いて、拡がりの計算だけを比べるためである。
//! 拡げた PENDF には、その温度の物質だけを入れておく（BROADR の出力は温度ごとに同じ MAT 番号の物質が並ぶ）。
//! 拡げる上限は ENDF の共鳴パラメータから決める。

use std::time::Instant;

use xs::doppler::{broaden, broadening_limit};
use xs::endf::{Mf2, Mf3, Tab1, Tape};
use xs::pointwise::DEFAULT_TOLERANCE;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let [_, endf, cold, hot, mat, temperature, mts @ ..] = args.as_slice() else {
        return Err(
            "usage: compare_njoy_broadr <ENDF> <PENDF 0 K> <PENDF hot> <MAT> <K> [MT...]".into(),
        );
    };
    let mat: i32 = mat.parse()?;
    let temperature: f64 = temperature.parse()?;
    let mts: Vec<u32> = if mts.is_empty() {
        vec![2, 18, 102]
    } else {
        mts.iter().map(|m| m.parse()).collect::<Result<_, _>>()?
    };
    let read = |path: &str| -> Result<Tape, Box<dyn std::error::Error>> {
        Ok(Tape::parse(&std::fs::read_to_string(path)?)?)
    };

    let evaluation = read(endf)?;
    let mf2 = match evaluation.section(mat, 2, 151) {
        Some(section) => Some(Mf2::parse(&section)?),
        None => None,
    };
    let cold = read(cold)?;
    let mut awr = 0.0;
    let mut present = Vec::new();
    let mut tables: Vec<Tab1> = Vec::new();
    for &mt in &mts {
        let Some(section) = cold.section(mat, 3, mt) else {
            println!("MT{mt}: missing in the 0 K PENDF");
            continue;
        };
        let mf3 = Mf3::parse(&section)?;
        awr = mf3.awr;
        present.push(mt);
        tables.push(mf3.sigma);
    }
    let e_max = broadening_limit(mf2.as_ref(), tables.iter());

    let start = Instant::now();
    let ours = broaden(&tables, awr, temperature, e_max, DEFAULT_TOLERANCE);
    println!(
        "broadened up to {e_max:.4e} eV in {:.1} s",
        start.elapsed().as_secs_f64()
    );

    let reference = read(hot)?;
    for (mt, table) in present.iter().zip(&ours) {
        let Some(section) = reference.section(mat, 3, *mt) else {
            println!("MT{mt}: missing in the broadened PENDF");
            continue;
        };
        let njoy = Mf3::parse(&section)?.sigma;
        let mut worst = (0.0_f64, 0.0_f64);
        let mut outside = 0;
        let mut compared = 0;
        for (&e, &expected) in njoy.x.iter().zip(&njoy.y) {
            if expected <= 0.0 || e >= e_max {
                continue;
            }
            let Some(actual) = table.eval(e) else {
                continue;
            };
            let rel = (actual - expected).abs() / expected;
            compared += 1;
            if rel > 2.0 * DEFAULT_TOLERANCE {
                outside += 1;
            }
            if rel > worst.0 {
                worst = (rel, e);
            }
        }
        let below = |t: &Tab1| t.x.iter().filter(|&&e| e < e_max).count();
        println!(
            "MT{mt}: points below limit ours {} njoy {}, compared {compared}, max rel diff {:.2e} at {:.6e} eV, > 2 tol: {outside}",
            below(table),
            below(&njoy),
            worst.0,
            worst.1
        );
    }
    Ok(())
}
