//! 格子計算を行い、少数群定数とセルごとの核分裂率を書き出すコマンド。
//!
//! 使い方: lattice <入力ファイル>
//!
//! 入力ファイル中のライブラリと出力のパスは、入力ファイルのあるディレクトリからの相対パス。
//!
//! k は標準出力に、経過（反復の回数、時間など）はログとして標準エラーに出す。ログの量は
//! 環境変数 `RUST_LOG` で選ぶ（既定は info。debug で反復ごとの k と変化、trace で最も変化の
//! 大きい領域まで）。
//!
//! スレッド数は、環境変数 `RAYON_NUM_THREADS` があればその値、なければ物理コアの数とする
//! （`taskset` やコンテナで使える CPU が絞られていれば、その数を超えない）。
//! MOC のスイープはメモリの転送速度で頭打ちになり、1つのコアを2つのスレッドで分ける
//! （SMT）とかえって遅くなるため。

use std::env::VarError;
use std::io::IsTerminal;
use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use anyhow::{Context, anyhow};
use tracing_subscriber::EnvFilter;
use tracing_subscriber::filter::LevelFilter;
use tracing_subscriber::fmt::format::FmtSpan;
use xlc_lattice::fewgroup::homogenize;
use xlc_lattice::input::{CrossSections, LatticeInput};
use xlc_lattice::macroscopic;
use xlc_lattice::material::Material;
use xlc_lattice::moc::{Options, solve};
use xlc_lattice::power::{cell_fission_rates, to_msgpack};
use xlc_lattice::tracking::Tracking;
use xlc_xs::library::Library;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.as_slice() {
        [input] => init_logging()
            .and_then(|()| use_physical_cores())
            .and_then(|()| run(Path::new(input))),
        _ => Err(anyhow!("usage: lattice <input>")),
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

/// `RAYON_NUM_THREADS` がなければ、rayon のスレッド数を物理コアの数にする。
/// 物理コアの数はマシン全体の値なので、このプロセスに許された CPU の数（論理）で抑える。
fn use_physical_cores() -> anyhow::Result<()> {
    if std::env::var_os("RAYON_NUM_THREADS").is_none() {
        let allowed = std::thread::available_parallelism().context("available CPUs")?;
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_cpus::get_physical().min(allowed.get()))
            .build_global()
            .context("thread pool")?;
    }
    Ok(())
}

fn run(path: &Path) -> anyhow::Result<()> {
    let start = Instant::now();
    let shown = |p: &Path| p.display().to_string();
    let text = std::fs::read_to_string(path).with_context(|| shown(path))?;
    let input = LatticeInput::parse(&text).with_context(|| shown(path))?;
    let base = path.parent().unwrap_or(Path::new("."));
    let (lattice, region_materials) = input.geometry().with_context(|| shown(path))?;

    let (named, groups) = match &input.cross_sections {
        CrossSections::Library(file) => {
            let file = base.join(file);
            let bytes = std::fs::read(&file).with_context(|| shown(&file))?;
            let library = Library::from_msgpack(&bytes).with_context(|| shown(&file))?;
            let named = input
                .materials
                .iter()
                .map(|m| {
                    Material::mix(&library, &m.composition, m.temperature)
                        .map(|mixed| (m.name.clone(), mixed))
                        .map_err(|e| anyhow!("material `{}`: {e}", m.name))
                })
                .collect::<anyhow::Result<Vec<_>>>()?;
            (named, Some(library.groups))
        }
        CrossSections::Macroscopic(file) => {
            let file = base.join(file);
            let text = std::fs::read_to_string(&file).with_context(|| shown(&file))?;
            (
                macroscopic::parse(&text).with_context(|| shown(&file))?,
                None,
            )
        }
    };
    let materials = region_materials
        .iter()
        .map(|name| {
            named
                .iter()
                .find(|(n, _)| n == name)
                .map(|(_, m)| m)
                .ok_or_else(|| anyhow!("material `{name}` is not in the cross sections"))
        })
        .collect::<anyhow::Result<Vec<&Material>>>()?;

    let tracking =
        Tracking::with_boundaries(&lattice, input.azimuthal, input.spacing, input.boundaries)?;
    let solution = solve(&tracking, &materials, &Options::default())?;
    if let Some((output, bounds)) = &input.output {
        let groups = groups
            .as_ref()
            .expect("the input checks that few-group output needs a library");
        let few = homogenize(&solution, &materials, groups, bounds)?;
        let output = base.join(output);
        std::fs::write(&output, few.to_msgpack()?).with_context(|| shown(&output))?;
    }
    if let Some(file) = &input.pin_powers {
        let rates = cell_fission_rates(&solution, &lattice, &materials)?;
        let file = base.join(file);
        std::fs::write(&file, to_msgpack(&rates, lattice.cells_per_side())?)
            .with_context(|| shown(&file))?;
    }
    println!("k {:.6}", solution.k);
    tracing::info!(
        iterations = solution.iterations,
        tracks = tracking.tracks.len(),
        seconds = start.elapsed().as_secs_f64(),
        "lattice"
    );
    Ok(())
}
