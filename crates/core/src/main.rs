//! 均質な直方体の炉心の固有値計算を行うコマンド。
//!
//! 使い方: core <入力ファイル>
//!
//! k_eff と漏れは標準出力に、経過（反復の回数、時間）はログとして標準エラーに出す
//! （`RUST_LOG=debug` で反復ごとの k と変化）。
//!
//! 入力ファイル中の少数群定数のパスは、入力ファイルのあるディレクトリからの相対パス。

use std::env::VarError;
use std::io::IsTerminal;
use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use anyhow::{Context, anyhow};
use tracing_subscriber::EnvFilter;
use tracing_subscriber::filter::LevelFilter;
use tracing_subscriber::fmt::format::FmtSpan;
use xlc_core::input::CoreInput;
use xlc_core::nem::{Mesh, Options, Problem, solve};
use xlc_lattice::fewgroup::FewGroup;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.as_slice() {
        [input] => init_logging().and_then(|()| run(Path::new(input))),
        _ => Err(anyhow!("usage: core <input>")),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error:#}");
            ExitCode::FAILURE
        }
    }
}

/// ログ（tracing）を標準エラーに出す。表示する量は環境変数 `RUST_LOG` で選ぶ（なければ info）。
/// 区間（span）が閉じるときに、その経過時間も出す。
fn init_logging() -> anyhow::Result<()> {
    // EnvFilter::from_env は UTF-8 でない値を黙って空とみなすので、自分で読む
    let directives = match std::env::var("RUST_LOG") {
        Ok(value) => value,
        Err(VarError::NotPresent) => String::new(),
        Err(error @ VarError::NotUnicode(_)) => return Err(anyhow!("RUST_LOG: {error}")),
    };
    // ParseError の表示は原因の文を含むので、原因の連鎖は付けない
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

fn run(path: &Path) -> anyhow::Result<()> {
    let start = Instant::now();
    let shown = |p: &Path| p.display().to_string();
    let text = std::fs::read_to_string(path).with_context(|| shown(path))?;
    let input = CoreInput::parse(&text).with_context(|| shown(path))?;
    let fewgroup_path = path
        .parent()
        .unwrap_or(Path::new("."))
        .join(&input.fewgroup);
    let bytes = std::fs::read(&fewgroup_path).with_context(|| shown(&fewgroup_path))?;
    let fewgroup = FewGroup::from_msgpack(&bytes).with_context(|| shown(&fewgroup_path))?;

    let split = |axis: usize| vec![input.size[axis] / input.nodes[axis] as f64; input.nodes[axis]];
    let problem = Problem {
        mesh: Mesh {
            dx: split(0),
            dy: split(1),
            dz: split(2),
        },
        materials: vec![fewgroup],
        node_material: vec![0; input.nodes.iter().product()],
        boundaries: input.boundaries,
    };
    let solution = solve(&problem, &Options::default())?;
    let leakage: f64 = solution.leakage.iter().sum();
    println!("k-eff {:.6}", solution.k);
    println!("leakage per fission neutron {:.4}", leakage * solution.k);
    tracing::info!(
        iterations = solution.iterations,
        seconds = start.elapsed().as_secs_f64(),
        "core"
    );
    Ok(())
}
