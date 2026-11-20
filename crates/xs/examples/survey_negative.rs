//! 評価ごとに、0 K の再構成（入力の宣言なし）で負の値が残るかを調べる（ADR-0267）。
//!
//! 使い方: survey_negative <ENDF>...
//!
//! ファイルごとに 1 行（タブ区切り）を書く: ファイル名と親のディレクトリ、結果（`OK`、負の値が残れば `NEG`、ほかの理由で
//! 止まれば `ERROR`、読めなければ `PARSE`）、`NEG` と `ERROR` は `ColdParts::read` の文。`NEG` の文は、負の値の残る反応ごとに
//! 点の数、エネルギーの範囲、最小の値を並べる。分離共鳴の弾性散乱の和は宣言なしで 0 にするので、数えない。

use xlc_xs::endf::Tape;
use xlc_xs::pointwise::{ColdParts, DEFAULT_TOLERANCE};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let paths: Vec<String> = std::env::args().skip(1).collect();
    if paths.is_empty() {
        return Err("usage: survey_negative <ENDF>...".into());
    }
    for path in &paths {
        let name: Vec<&str> = path.rsplit('/').take(2).collect();
        let name = name.join(" ");
        let tape = match Tape::parse(&std::fs::read_to_string(path)?) {
            Ok(tape) => tape,
            Err(error) => {
                println!("{name}\tPARSE\t{error}");
                continue;
            }
        };
        let Some(&mat) = tape.materials().first() else {
            println!("{name}\tPARSE\tno material");
            continue;
        };
        match ColdParts::read(&tape, mat, DEFAULT_TOLERANCE, &[]) {
            Ok(_) => println!("{name}\tOK"),
            // 負の値の文は、入力の宣言の項目の名前を示す
            Err(error) if error.message.contains("negative-to-zero") => {
                println!("{name}\tNEG\t{error}");
            }
            Err(error) => println!("{name}\tERROR\t{error}"),
        }
    }
    Ok(())
}
