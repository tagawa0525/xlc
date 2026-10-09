//! 炉心の固有値計算を行うコマンド。
//!
//! 使い方: core <ケースのファイル>（ケース、サイクル、炉心の形、台帳の4つの入力。ADR-0105）
//!
//! k_eff、漏れ、集合体の位置ごとの核分裂率は標準出力に、経過（反復の回数、時間）はログとして標準エラーに
//! 出す（`RUST_LOG=debug` で反復ごとの k と変化）。
//!
//! 入力のパスは、それを書いたファイルのあるディレクトリからの相対パス。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

use std::env::VarError;
use std::io::IsTerminal;
use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use anyhow::{Context, anyhow};
use tracing_subscriber::EnvFilter;
use tracing_subscriber::filter::LevelFilter;
use tracing_subscriber::fmt::format::FmtSpan;
use xlc_core::input::{Case, Cycle, InputError, Inventory, Layout, Plant};
use xlc_core::nem::{Options, Problem, Solution, solve};
use xlc_format::fewgroup::FewGroup;

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
    let case = parse(path, Case::parse)?;
    let cycle = parse(&case.cycle, Cycle::parse)?;
    let plant = parse(&cycle.plant, Plant::parse)?;
    let inventory = parse(&cycle.inventory, Inventory::parse)?;
    let layout = Layout::new(&case, &cycle, &plant, &inventory).with_context(|| shown(path))?;
    let constants = layout
        .files
        .iter()
        .map(|file| {
            let bytes = std::fs::read(file).with_context(|| shown(file))?;
            FewGroup::from_msgpack(&bytes).with_context(|| shown(file))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    let problem = layout.problem(&constants).with_context(|| shown(path))?;
    let solution = solve(&problem, &Options::default())?;
    let leakage: f64 = solution.leakage.iter().sum();
    println!("k-eff {:.6}", solution.k);
    println!("leakage per fission neutron {:.4}", leakage * solution.k);
    println!("assembly fission rates (relative to the assembly mean, top row is +y)");
    for row in assembly_fission_rates(&layout, &problem, &solution) {
        let cells: Vec<String> = row
            .iter()
            .map(|r| r.map_or_else(|| format!("{:>6}", "."), |r| format!("{r:6.4}")))
            .collect();
        println!("{}", cells.join(" "));
    }
    tracing::info!(
        iterations = solution.iterations,
        seconds = start.elapsed().as_secs_f64(),
        "core"
    );
    Ok(())
}

fn shown(p: &Path) -> String {
    p.display().to_string()
}

/// ファイルを読み、ファイルのあるディレクトリを起点にして解釈する。
fn parse<T>(
    path: &Path,
    parse: impl FnOnce(&str, &Path) -> Result<T, InputError>,
) -> anyhow::Result<T> {
    let text = std::fs::read_to_string(path).with_context(|| shown(path))?;
    parse(&text, path.parent().unwrap_or(Path::new("."))).with_context(|| shown(path))
}

/// 集合体の位置ごとの、軸方向に積分した核分裂率（集合体の位置の平均で割る）。上の行（+y）から。
fn assembly_fission_rates(
    layout: &Layout,
    problem: &Problem,
    solution: &Solution,
) -> Vec<Vec<Option<f64>>> {
    let (nx, ny) = (problem.mesh.dx.len(), problem.mesh.dy.len());
    let mut rates = vec![vec![0.0; nx]; ny];
    for (n, flux) in solution.flux.iter().enumerate() {
        let (x, y, z) = (n % nx, (n / nx) % ny, n / (nx * ny));
        let m = &problem.materials[problem.node_material[n]];
        let volume = problem.mesh.dx[x] * problem.mesh.dy[y] * problem.mesh.dz[z];
        let fission: f64 = m.fission.iter().zip(flux).map(|(f, phi)| f * phi).sum();
        // 図の上の行が +y
        rates[ny - 1 - y][x] += fission * volume;
    }
    let loaded: Vec<f64> = layout
        .ids
        .iter()
        .zip(&rates)
        .flat_map(|(ids, r)| {
            ids.iter()
                .zip(r)
                .filter(|(id, _)| id.is_some())
                .map(|(_, &r)| r)
        })
        .collect();
    let mean = loaded.iter().sum::<f64>() / loaded.len().max(1) as f64;
    layout
        .ids
        .iter()
        .zip(&rates)
        .map(|(ids, r)| {
            ids.iter()
                .zip(r)
                .map(|(id, &r)| id.as_ref().map(|_| r / mean))
                .collect()
        })
        .collect()
}
