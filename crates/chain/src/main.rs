//! ENDF から全核種の燃焼鎖を作るコマンド（ADR-0164、ADR-0166）。
//!
//! 使い方: chain <入力ファイル>
//!
//! 入力ファイルは TOML（パスは入力ファイルのディレクトリからの相対パス）:
//!
//! ```toml
//! output = "chain.msgpack"                                         # 書く鎖（xlc_depletion::chain の形式）
//! decay-directory = "../../references/endf/ENDF-B-VIII.0_decay"    # decay の評価（*.endf をすべて読む）
//! nfy-directory = "../../references/endf/ENDF-B-VIII.0_nfy"        # 中性子の核分裂の収率
//! neutron-directory = "../../references/endf/ENDF-B-VIII.0_neutrons"  # 中性子の反応
//! ```
//!
//! データにないものを規則で補った箇所は警告（tracing の warn）として標準エラーに出す。

use std::env::VarError;
use std::io::IsTerminal;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{Context, Result, anyhow};
use rayon::prelude::*;
use serde::Deserialize;
use tracing_subscriber::EnvFilter;
use tracing_subscriber::filter::LevelFilter;
use tracing_subscriber::fmt::format::FmtSpan;
use xlc_chain::build::build;
use xlc_chain::endf::{read_decay, read_neutron, read_yields};
// MF8 の生成物の状態も、鎖の反応の MT で読む
use xlc_depletion::chain::REACTION_MTS;

const USAGE: &str = "usage: chain <input.toml>";

/// 入力ファイルの中身（パスは入力ファイルのディレクトリからの相対）
#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct Input {
    output: PathBuf,
    decay_directory: PathBuf,
    nfy_directory: PathBuf,
    neutron_directory: PathBuf,
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = init_logging().and_then(|()| match &args[..] {
        [path] => run(Path::new(path)),
        _ => Err(anyhow!(USAGE)),
    });
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error:#}");
            ExitCode::FAILURE
        }
    }
}

/// ログ（tracing）を標準エラーに出す。表示する量は環境変数 `RUST_LOG` で選ぶ（なければ info）。
#[allow(
    clippy::disallowed_methods,
    reason = "RUST_LOG は表示する量だけを選び、計算の結果を変えない（ADR-0137）"
)]
fn init_logging() -> Result<()> {
    // EnvFilter::from_env は UTF-8 でない値を黙って空とみなすので、自分で読む
    let directives = match std::env::var("RUST_LOG") {
        Ok(value) => value,
        Err(VarError::NotPresent) => String::new(),
        Err(error @ VarError::NotUnicode(_)) => return Err(anyhow!("RUST_LOG: {error}")),
    };
    let filter = EnvFilter::builder()
        .with_default_directive(LevelFilter::INFO.into())
        .parse(directives)
        .map_err(|error| anyhow!("RUST_LOG: {error}"))?;
    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_writer(std::io::stderr)
        .with_ansi(std::io::stderr().is_terminal())
        .with_span_events(FmtSpan::CLOSE)
        .init();
    Ok(())
}

/// ディレクトリの `*.endf` を名前の順に返す
fn endf_files(directory: &Path) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for entry in std::fs::read_dir(directory).with_context(|| format!("{}", directory.display()))? {
        let path = entry?.path();
        if path.extension().is_some_and(|e| e == "endf") {
            files.push(path);
        }
    }
    files.sort();
    if files.is_empty() {
        return Err(anyhow!("{}: no *.endf files", directory.display()));
    }
    Ok(files)
}

/// ファイルをすべて並列に読み、`read` で評価にする
fn read_all<T: Send>(directory: &Path, read: impl Fn(&str) -> Result<T> + Sync) -> Result<Vec<T>> {
    endf_files(directory)?
        .par_iter()
        .map(|path| {
            let text =
                std::fs::read_to_string(path).with_context(|| format!("{}", path.display()))?;
            read(&text).with_context(|| format!("{}", path.display()))
        })
        .collect()
}

fn run(path: &Path) -> Result<()> {
    let text = std::fs::read_to_string(path).with_context(|| format!("{}", path.display()))?;
    let input: Input = toml::from_str(&text).with_context(|| format!("{}", path.display()))?;
    let base = path.parent().unwrap_or(Path::new("."));
    let decays = read_all(&base.join(&input.decay_directory), read_decay)?;
    let yields = read_all(&base.join(&input.nfy_directory), read_yields)?;
    let neutrons = read_all(&base.join(&input.neutron_directory), |text| {
        read_neutron(text, &REACTION_MTS)
    })?;
    tracing::info!(
        decay = decays.len(),
        nfy = yields.len(),
        neutrons = neutrons.len(),
        "read the evaluations"
    );
    let built = build(&decays, &neutrons, &yields)?;
    for warning in &built.warnings {
        tracing::warn!("{warning}");
    }
    let chain = &built.chain;
    tracing::info!(
        nuclides = chain.nuclides.len(),
        unstable = chain.nuclides.iter().filter(|n| n.decay.is_some()).count(),
        with_reactions = chain
            .nuclides
            .iter()
            .filter(|n| !n.reactions.is_empty())
            .count(),
        fissionable = chain
            .nuclides
            .iter()
            .filter(|n| n.fission.is_some())
            .count(),
        warnings = built.warnings.len(),
        "built the chain"
    );
    let output = base.join(&input.output);
    std::fs::write(&output, chain.to_msgpack()?)
        .with_context(|| format!("{}", output.display()))?;
    Ok(())
}
