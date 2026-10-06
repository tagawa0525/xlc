//! 格子計算を行い、少数群定数とセルごとの核分裂率を書き出すコマンド。
//!
//! 使い方: lattice <入力ファイル>
//!
//! 入力ファイル中のライブラリと出力のパスは、入力ファイルのあるディレクトリからの相対パス。

use std::path::Path;
use std::process::ExitCode;

use lattice::input::LatticeInput;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [input] = args.as_slice() else {
        eprintln!("usage: lattice <input>");
        return ExitCode::FAILURE;
    };
    match run(Path::new(input)) {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("error: {message}");
            ExitCode::FAILURE
        }
    }
}

fn run(path: &Path) -> Result<(), String> {
    let read = |p: &Path| std::fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()));
    let input =
        LatticeInput::parse(&read(path)?).map_err(|e| format!("{}: {e}", path.display()))?;
    todo!("{input:?}")
}
