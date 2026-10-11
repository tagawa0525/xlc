//! 炉心の固有値計算を行うコマンド。
//!
//! 使い方: core <ケースのファイル>（ケース、サイクル、炉心の形、台帳の4つの入力。ADR-0105）
//!
//! k_eff、漏れ、集合体の位置ごとの核分裂率は標準出力に、経過（反復の回数、時間）はログとして標準エラーに
//! 出す（`RUST_LOG=debug` で反復ごとの k と変化）。ケースに名前つきの状態があれば、状態を並列に解き、名前の順に
//! `state <名前>` の行に続けて状態ごとに出す（探索した状態は求めたホウ素濃度も。ADR-0145、ADR-0146）。
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
use xlc_core::case::{CaseError, check_search, solve_at, solve_states};
use xlc_core::input::{Case, Cycle, InputError, Inventory, Layout, Plant};
use xlc_core::nem::{Options, Problem, Solution};
use xlc_core::pins::{axial_rows, core_rows, reconstruct};
use xlc_format::pins::{AxialPinPowers, PinPowers};
use xlc_format::table::Table;

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
#[allow(
    clippy::disallowed_methods,
    reason = "RUST_LOG は表示する量だけを選び、計算の結果を変えない（ADR-0137）"
)]
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
    // 4つの入力の組み合わせの誤りなので、どのファイルかを並べて示す
    let inputs = || {
        format!(
            "{} (cycle {}, plant {}, inventory {})",
            shown(path),
            shown(&case.cycle),
            shown(&cycle.plant),
            shown(&cycle.inventory)
        )
    };
    let layout = Layout::new(&case, &cycle, &plant, &inventory).with_context(inputs)?;
    let tables = layout
        .files
        .iter()
        .map(|file| {
            let bytes = std::fs::read(file).with_context(|| shown(file))?;
            Table::from_msgpack(&bytes).with_context(|| shown(file))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    let options = Options::default();
    if case.states.is_empty() {
        let s = solve_at(&layout, &tables, None, &options).map_err(|e| in_inputs(e, inputs()))?;
        report(&layout, &s.problem, &s.solution);
        write_pin_powers(&case, &s.problem, &s.solution)?;
    } else {
        // 状態をすべて並列に解いてから、名前の順に出す（ADR-0146）
        // 解く前にできる検査は、長い計算を始める前に済ませる
        let invalid: Vec<String> = case
            .states
            .iter()
            .filter_map(|state| {
                check_search(state, &options)
                    .err()
                    .map(|e| format!("state {}: {e}", state.name))
            })
            .collect();
        if !invalid.is_empty() {
            return Err(anyhow!("{}", invalid.join("\n")));
        }
        let results = solve_states(&layout, &tables, &case.states, &options);
        // 解けた状態は出し、解けなかった状態のエラーは最後にまとめて返す（長い計算のほかの状態の結果を捨てない）
        let mut errors = Vec::new();
        for (state, result) in case.states.iter().zip(results) {
            match result {
                Ok(s) => {
                    let _span = tracing::info_span!("report", state = %state.name).entered();
                    println!("state {}", state.name);
                    if let Some(searched) = s.searched {
                        println!(
                            "boron {:.4} ppm (searched in {} solves)",
                            searched.ppm, searched.solves
                        );
                    }
                    report(&layout, &s.problem, &s.solution);
                    write_pin_powers(&case, &s.problem, &s.solution)?;
                }
                Err(e) => errors.push(format!(
                    "state {}: {:#}",
                    state.name,
                    in_inputs(e, inputs())
                )),
            }
        }
        if !errors.is_empty() {
            return Err(anyhow!("{}", errors.join("\n")));
        }
    }
    tracing::info!(seconds = start.elapsed().as_secs_f64(), "core");
    Ok(())
}

/// 4つの入力の組み合わせの誤り（`CaseError::Input`）にだけ、どのファイルかを添える。内挿や固有値の反復の
/// エラーは入力の組み合わせの誤りではないので、そのまま返す
fn in_inputs(error: CaseError, inputs: String) -> anyhow::Error {
    match error {
        CaseError::Input(e) => anyhow::Error::new(e).context(inputs),
        e => e.into(),
    }
}

/// k_eff、漏れ、集合体の位置ごとの核分裂率を標準出力に出す。
fn report(layout: &Layout, problem: &Problem, solution: &Solution) {
    let leakage: f64 = solution.leakage.iter().sum();
    // 状態の間の k の差（価値、係数）を出力から取るので、丸めの誤差が固有値の反復の判定（1e-8）より小さくなる
    // 桁まで出す。k の誤差そのものは反復の収束で決まる（ADR-0146）
    println!("k-eff {:.8}", solution.k);
    println!("leakage per fission neutron {:.4}", leakage * solution.k);
    // メッシュの違う計算を面の和で比べられるように、z の面の境界と核分裂率の割合（下から、和が 1）。集合体の表の
    // 後ろに置かないのは、比べるツール（examples/vera/compare_p5.py）が表の見出しから後をすべて表の行として読むため
    println!("axial fission fractions (bottom first: z from, z to (cm), fraction)");
    let mut z = 0.0;
    for (fraction, dz) in layout
        .axial_fission_fractions(problem, solution)
        .iter()
        .zip(&problem.mesh.dz)
    {
        println!("{z:.6} {:.6} {fraction:.8e}", z + dz);
        z += dz;
    }
    println!("assembly fission rates (relative to the assembly mean, top row is +y)");
    println!("(`.` is a reflector site, `-` a position without a node)");
    for (row, empty) in layout
        .assembly_fission_rates(problem, solution)
        .iter()
        .zip(&layout.empty)
    {
        let cells: Vec<String> = row
            .iter()
            .zip(empty)
            .map(|(r, &empty)| match r {
                Some(r) => format!("{r:6.4}"),
                None => format!("{:>6}", if empty { "-" } else { "." }),
            })
            .collect();
        println!("{}", cells.join(" "));
    }
    // 格子計算の均質化のノード平均と比べられるように、直方体の位置の番号（x が最も速く、次に y、z）で出す（全体の
    // 生成率 1 に規格化）。ノードのない位置は出さない
    for (node, flux) in solution.flux.iter().enumerate() {
        if !flux.is_empty() {
            tracing::debug!(node, flux = ?flux, "node average flux");
        }
    }
    tracing::info!(iterations = solution.iterations, "eigenvalue");
}

/// 再構成したピン出力（ADR-0096）を、ケースが指定していれば書く。全体で積分したものと層ごとのもの
fn write_pin_powers(case: &Case, problem: &Problem, solution: &Solution) -> anyhow::Result<()> {
    if case.pin_powers.is_none() && case.axial_pin_powers.is_none() {
        return Ok(());
    }
    let reconstructed = Instant::now();
    let maps = reconstruct(problem, solution).map_err(|e| anyhow!("pin powers: {e}"))?;
    if let Some(file) = &case.pin_powers {
        let rows = core_rows(problem, &maps).map_err(|e| anyhow!("pin powers: {e}"))?;
        let bytes = PinPowers { rows }
            .to_msgpack()
            .with_context(|| shown(file))?;
        std::fs::write(file, bytes).with_context(|| shown(file))?;
    }
    if let Some(file) = &case.axial_pin_powers {
        let layers = axial_rows(problem, &maps).map_err(|e| anyhow!("axial pin powers: {e}"))?;
        let bytes = AxialPinPowers {
            bounds: layers.bounds,
            axial: layers.axial,
            layers: layers.layers,
        }
        .to_msgpack()
        .with_context(|| shown(file))?;
        std::fs::write(file, bytes).with_context(|| shown(file))?;
    }
    tracing::info!(
        seconds = reconstructed.elapsed().as_secs_f64(),
        "pin power reconstruction"
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
