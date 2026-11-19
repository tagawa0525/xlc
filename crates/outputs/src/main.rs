//! 出力を読むだけの入口（ADR-0258）。
//!
//! 使い方: read-outputs < <一覧>
//!
//! 標準入力の 1 行が 1 つのファイルで、`<種類>\t<パス>`（種類は `msgpack` か `xs-range`）。読めないファイルごとに
//! `UNREADABLE\t<パス>\t<理由>` を標準出力に 1 行で出す。終了コードは、すべて読めれば 0、読めないファイルがあれば 1、
//! 一覧の行が読めなければ 2（そこで止まる）。ファイルは 1 つずつ読む（大きいライブラリは 1 GB を超え、並べると記憶が足りない）。

use std::io::{BufRead, Write};
use std::path::Path;
use std::process::ExitCode;

use anyhow::{Context, anyhow};
use xlc_outputs::{Kind, read};

fn main() -> ExitCode {
    match run() {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::from(1),
        Err(e) => {
            let _ = writeln!(std::io::stderr(), "error: {e:#}");
            ExitCode::from(2)
        }
    }
}

/// すべて読めたか
fn run() -> anyhow::Result<bool> {
    let mut all = true;
    let mut out = std::io::stdout().lock();
    for (i, line) in std::io::stdin().lock().lines().enumerate() {
        let line = line.context("cannot read the list")?;
        let (kind, path) = line
            .split_once('\t')
            .ok_or_else(|| anyhow!("line {}: expected <kind>\\t<path>, got {line:?}", i + 1))?;
        let kind: Kind = kind.parse().with_context(|| format!("line {}", i + 1))?;
        if let Err(e) = read_file(kind, Path::new(path)) {
            all = false;
            let reason = format!("{e:#}").replace(['\n', '\t'], " ");
            writeln!(out, "UNREADABLE\t{path}\t{reason}")?;
        }
    }
    Ok(all)
}

fn read_file(kind: Kind, path: &Path) -> anyhow::Result<()> {
    let bytes = std::fs::read(path).with_context(|| format!("cannot open {}", path.display()))?;
    read(kind, &bytes)
}
