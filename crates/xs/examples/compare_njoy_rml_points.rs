//! R 行列の制限形（LRF=7）の分離共鳴の断面積を、線形化せずに NJOY (RECONR) の PENDF の点で比べる。
//!
//! 使い方: compare_njoy_rml_points <ENDF> <PENDF> <MAT>
//!
//! LRF=7 の範囲の中の PENDF の各エネルギー点で、`RMatrixLimitedXs` の値に MF3 のバックグラウンドを足したものと
//! NJOY の値の相対差を求め、反応（MT1、2、102）ごとに最大値を表示する。`compare_njoy_pendf` は線形化した表を比べるので、
//! 差が R 行列の式から来るか、線形化の格子から来るかをこれで分ける。

use xlc_xs::endf::{Mf2, Mf3, RangeKind, Tape};
use xlc_xs::resonance::RMatrixLimitedXs;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let [_, endf, pendf, mat] = args.as_slice() else {
        return Err("usage: compare_njoy_rml_points <ENDF> <PENDF> <MAT>".into());
    };
    let mat: i32 = mat.parse()?;
    let tape = Tape::parse(&std::fs::read_to_string(endf)?)?;
    let section = tape.section(mat, 2, 151).ok_or("no MF2 MT151")?;
    let mf2 = Mf2::parse(&section)?;
    let (range, rml) = mf2
        .ranges()
        .find_map(|r| match &r.kind {
            RangeKind::RMatrixLimited(rml) => Some((r, rml)),
            _ => None,
        })
        .ok_or("no LRF=7 range")?;
    let xs = RMatrixLimitedXs::new(rml)?;
    let reference = Tape::parse(&std::fs::read_to_string(pendf)?)?;
    // `cross_sections` の並び（全断面積、弾性散乱、捕獲）
    for (column, mt) in [1, 2, 102].into_iter().enumerate() {
        let background = match tape.section(mat, 3, mt) {
            Some(s) => Some(Mf3::parse(&s)?.sigma),
            None => None,
        };
        let njoy = Mf3::parse(&reference.section(mat, 3, mt).ok_or("no MT in PENDF")?)?.sigma;
        let (mut compared, mut worst, mut at) = (0, 0.0_f64, 0.0);
        for (&e, &expected) in njoy.x.iter().zip(&njoy.y) {
            if e < range.el || e >= range.eh || expected <= 0.0 {
                continue;
            }
            let bkg = background.as_ref().and_then(|b| b.eval(e)).unwrap_or(0.0);
            let actual = xs.cross_sections(e)?[column] + bkg;
            let rel = ((actual - expected) / expected).abs();
            compared += 1;
            if rel > worst {
                (worst, at) = (rel, e);
            }
        }
        println!("MT{mt}: compared {compared}, max rel diff {worst:.2e} at {at:.6e} eV");
    }
    Ok(())
}
