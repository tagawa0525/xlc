//! 0 K のポイントワイズの断面積を、NJOY (RECONR) で作った PENDF と比べる。
//!
//! 使い方: compare_njoy_pendf <ENDF> <PENDF> <MAT> [MT...]
//!
//! PENDF の各エネルギー点で、こちらの断面積との相対差を求め、反応ごとに最大値と、
//! 相対差が許容誤差の2倍を超える点の割合を表示する。

use std::time::Instant;

use xs::endf::{Mf3, Tape};
use xs::pointwise::{DEFAULT_TOLERANCE, reconstruct};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let [_, endf, pendf, mat, mts @ ..] = args.as_slice() else {
        return Err("usage: compare_njoy_pendf <ENDF> <PENDF> <MAT> [MT...]".into());
    };
    let mat: i32 = mat.parse()?;
    let mts: Vec<u32> = if mts.is_empty() {
        vec![1, 2, 18, 102]
    } else {
        mts.iter().map(|m| m.parse()).collect::<Result<_, _>>()?
    };

    let start = Instant::now();
    let tape = Tape::parse(&std::fs::read_to_string(endf)?)?;
    let ours = reconstruct(&tape, mat, DEFAULT_TOLERANCE)?;
    println!("reconstructed in {:.1} s", start.elapsed().as_secs_f64());

    let reference = Tape::parse(&std::fs::read_to_string(pendf)?)?;
    for mt in mts {
        let (Some(table), Some(section)) = (ours.get(&mt), reference.section(mat, 3, mt)) else {
            println!("MT{mt}: missing");
            continue;
        };
        let njoy = Mf3::parse(&section)?.sigma;
        let mut worst = (0.0_f64, 0.0_f64);
        let mut outside = 0;
        let mut compared = 0;
        for (&e, &expected) in njoy.x.iter().zip(&njoy.y) {
            if expected <= 0.0 {
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
        println!(
            "MT{mt}: points ours {} njoy {}, compared {compared}, max rel diff {:.2e} at {:.6e} eV, > 2 tol: {outside}",
            table.x.len(),
            njoy.x.len(),
            worst.0,
            worst.1
        );
    }
    Ok(())
}
