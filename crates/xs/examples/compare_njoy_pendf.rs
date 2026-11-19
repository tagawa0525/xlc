//! 0 K のポイントワイズの断面積を、NJOY (RECONR) で作った PENDF と比べる。
//!
//! 使い方: compare_njoy_pendf <ENDF> <PENDF> <MAT> [MT...]
//!
//! PENDF の各エネルギー点で、こちらの断面積との相対差を求め、反応ごとに最大値と、
//! 相対差が許容誤差の2倍を超える点の数を表示する。MF2 に共鳴の範囲があれば、範囲ごと（分離共鳴、非分離共鳴、
//! 範囲の外）にも同じ量を表示する。超える点があれば、その点の NJOY の値と差の絶対値の最大も表示する（弾性散乱が
//! 0 の近くで打ち消し合う点では、NJOY の有効数字の丸めが相対差を大きくする）。

use std::time::Instant;

use xlc_xs::endf::{Mf2, Mf3, RangeKind, Tab1, Tape};
use xlc_xs::pointwise::{DEFAULT_TOLERANCE, reconstruct};

/// 比べた点の集計
#[derive(Default)]
struct Stats {
    compared: usize,
    /// 相対差の最大と、そのエネルギー
    worst: (f64, f64),
    /// 許容誤差の 2 倍を超える点の数と、その点の NJOY の値と差の絶対値の最大
    outside: usize,
    outside_value: f64,
    outside_diff: f64,
}

impl Stats {
    fn add(&mut self, e: f64, actual: f64, expected: f64) {
        let rel = (actual - expected).abs() / expected;
        self.compared += 1;
        if rel > 2.0 * DEFAULT_TOLERANCE {
            self.outside += 1;
            self.outside_value = self.outside_value.max(expected);
            self.outside_diff = self.outside_diff.max((actual - expected).abs());
        }
        if rel > self.worst.0 {
            self.worst = (rel, e);
        }
    }

    fn summary(&self) -> String {
        let detail = if self.outside > 0 {
            format!(
                " (njoy up to {:.2e} b, |diff| up to {:.2e} b)",
                self.outside_value, self.outside_diff
            )
        } else {
            String::new()
        };
        format!(
            "compared {}, max rel diff {:.2e} at {:.6e} eV, > 2 tol: {}{detail}",
            self.compared, self.worst.0, self.worst.1, self.outside
        )
    }
}

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
    // 共鳴の範囲 [el, eh) と名前。散乱半径だけの範囲は分けない
    let ranges: Vec<(f64, f64, &str)> = match tape.section(mat, 2, 151) {
        Some(section) => Mf2::parse(&section)?
            .ranges
            .iter()
            .filter_map(|r| {
                let name = match r.kind {
                    RangeKind::ScatteringRadius { .. } => return None,
                    RangeKind::Unresolved(_) => "unresolved",
                    _ => "resolved",
                };
                Some((r.el, r.eh, name))
            })
            .collect(),
        None => Vec::new(),
    };

    let reference = Tape::parse(&std::fs::read_to_string(pendf)?)?;
    for mt in mts {
        let (Some(table), Some(section)) = (ours.get(&mt), reference.section(mat, 3, mt)) else {
            println!("MT{mt}: missing");
            continue;
        };
        let njoy = Mf3::parse(&section)?.sigma;
        let (all, by_range) = compare(table, &njoy, &ranges);
        println!(
            "MT{mt}: points ours {} njoy {}, {}",
            table.x.len(),
            njoy.x.len(),
            all.summary()
        );
        if ranges.is_empty() {
            continue;
        }
        let names = ranges
            .iter()
            .map(|&(el, eh, name)| format!("{name} {el:.6e}..{eh:.6e} eV"))
            .chain(["outside the resonance ranges".to_string()]);
        for (name, stats) in names.zip(&by_range) {
            if stats.compared > 0 {
                println!("  {name}: {}", stats.summary());
            }
        }
    }
    Ok(())
}

/// NJOY の値が正の点で比べた、全体と範囲ごと（最後は範囲の外）の集計。
fn compare(table: &Tab1, njoy: &Tab1, ranges: &[(f64, f64, &str)]) -> (Stats, Vec<Stats>) {
    let mut all = Stats::default();
    let mut by_range: Vec<Stats> = (0..=ranges.len()).map(|_| Stats::default()).collect();
    for (&e, &expected) in njoy.x.iter().zip(&njoy.y) {
        if expected <= 0.0 {
            continue;
        }
        let Some(actual) = table.eval(e) else {
            continue;
        };
        all.add(e, actual, expected);
        let index = ranges
            .iter()
            .position(|&(el, eh, _)| e >= el && e < eh)
            .unwrap_or(ranges.len());
        by_range[index].add(e, actual, expected);
    }
    (all, by_range)
}
