//! 格子計算を行い、少数群定数のテーブルとセルごとの核分裂率を書き出すコマンド。入力に分岐（`[branches]`）が
//! あれば、状態ごとに格子計算を解き、軸の直積の順の状態をテーブルにする（ADR-0097）。
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
use xlc_format::fewgroup::FewGroup;
use xlc_format::library::Library;
use xlc_format::table::Table;
use xlc_lattice::fewgroup::{Node, homogenize, homogenize_nodes};
use xlc_lattice::input::{CrossSections, FewGroupFiles, FewGroupOutput, LatticeInput, VOID};
use xlc_lattice::macroscopic;
use xlc_lattice::material::Material;
use xlc_lattice::moc::{Options, available_order, solve};
use xlc_lattice::power::{cell_fission_rates, to_msgpack};
use xlc_lattice::shielding::{ShieldingOptions, shield_fuel};
use xlc_lattice::tracking::Tracking;

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

    let tracking =
        Tracking::with_boundaries(&lattice, input.azimuthal, input.spacing, input.boundaries)?;
    // 物質の一覧と、領域ごとのその番号
    let by_name = |mut named: Vec<(String, Material)>,
                   groups: usize|
     -> anyhow::Result<(Vec<Material>, Vec<usize>)> {
        // 真空の帯の物質。名前は予約語で、断面積のファイルにも同じ名前の物質があってはならない
        if named.iter().any(|(n, _)| n == VOID) {
            return Err(anyhow!(
                "`{VOID}` is reserved for void zones and cannot be a material"
            ));
        }
        named.push((VOID.into(), Material::void(groups)));
        let index = region_materials
            .iter()
            .map(|name| {
                named
                    .iter()
                    .position(|(n, _)| n == name)
                    .ok_or_else(|| anyhow!("material `{name}` is not in the cross sections"))
            })
            .collect::<anyhow::Result<Vec<_>>>()?;
        Ok((named.into_iter().map(|(_, m)| m).collect(), index))
    };
    let library = match &input.cross_sections {
        CrossSections::Library(file) => {
            let file = base.join(file);
            let bytes = std::fs::read(&file).with_context(|| shown(&file))?;
            Some(Library::from_msgpack(&bytes).with_context(|| shown(&file))?)
        }
        CrossSections::Macroscopic(_) => None,
    };
    let groups = library.as_ref().map(|l| l.groups.clone());

    // 分岐の状態ごとに格子計算を解き、出力ごとの少数群定数を状態の順に集める（ADR-0097）
    let states = input.branch_states();
    let mut tables: Vec<Vec<FewGroup>> = Vec::new();
    for (state, state_materials) in states.iter().enumerate() {
        let (unique, region_index) = match (&input.cross_sections, &library) {
            (CrossSections::Library(_), Some(library)) => {
                if input.self_shielding {
                    let shielded = shield_fuel(
                        library,
                        &tracking,
                        state_materials,
                        &region_materials,
                        &ShieldingOptions::default(),
                    )
                    .with_context(|| format!("{}: self-shielding", shown(path)))?;
                    (shielded.materials, shielded.region_material)
                } else {
                    let named = state_materials
                        .iter()
                        .map(|m| {
                            Material::mix(library, &m.composition, m.temperature)
                                .map(|mixed| (m.name.clone(), mixed))
                                .map_err(|e| anyhow!("material `{}`: {e}", m.name))
                        })
                        .collect::<anyhow::Result<Vec<_>>>()?;
                    by_name(named, library.groups.len())?
                }
            }
            (CrossSections::Macroscopic(file), _) => {
                let file = base.join(file);
                let text = std::fs::read_to_string(&file).with_context(|| shown(&file))?;
                let named = macroscopic::parse(&text).with_context(|| shown(&file))?;
                let groups = named.first().map_or(0, |(_, m)| m.groups());
                by_name(named, groups)?
            }
            (CrossSections::Library(_), None) => return Err(anyhow!("library is not read")),
        };
        let materials: Vec<&Material> = region_index.iter().map(|&i| &unique[i]).collect();

        let options = Options {
            scattering_order: input
                .scattering_order
                .unwrap_or_else(|| available_order(&materials)),
            ..Options::default()
        };
        let solution = solve(&tracking, &materials, &options)?;
        if let Some(FewGroupOutput { bounds, files }) = &input.output {
            let groups = groups
                .as_ref()
                .ok_or_else(|| anyhow!("few-group output needs a library"))?;
            let boundaries = tracking.boundaries();
            let fews = match files {
                FewGroupFiles::Lattice(_) => vec![homogenize(
                    &solution,
                    &materials,
                    groups,
                    bounds,
                    &boundaries,
                )?],
                FewGroupFiles::Nodes(entries) => {
                    let nodes: Vec<Node> = entries.iter().map(|(_, node)| node.clone()).collect();
                    homogenize_nodes(
                        &solution,
                        &lattice,
                        &materials,
                        groups,
                        bounds,
                        &boundaries,
                        &nodes,
                    )?
                }
            };
            tables.resize_with(fews.len(), Vec::new);
            for (table, few) in tables.iter_mut().zip(fews) {
                table.push(few);
            }
        }
        // ピン出力は分岐のないときだけ（入力の読み込みで確かめてある）
        if let Some(file) = &input.pin_powers {
            let rates = cell_fission_rates(&solution, &lattice, &materials)?;
            let file = base.join(file);
            std::fs::write(&file, to_msgpack(&rates, &lattice)?).with_context(|| shown(&file))?;
        }
        if input.branches.is_some() {
            println!("state {state} k {:.6}", solution.k);
        } else {
            println!("k {:.6}", solution.k);
        }
        tracing::info!(
            state,
            iterations = solution.iterations,
            tracks = tracking.tracks.len(),
            seconds = start.elapsed().as_secs_f64(),
            "lattice"
        );
    }

    if let Some(FewGroupOutput { files, .. }) = &input.output {
        let outputs: Vec<&Path> = match files {
            FewGroupFiles::Lattice(output) => vec![output],
            FewGroupFiles::Nodes(entries) => entries.iter().map(|(o, _)| o.as_path()).collect(),
        };
        let axes = input
            .branches
            .as_ref()
            .map(|b| b.axes())
            .unwrap_or_default();
        for (output, states) in outputs.into_iter().zip(tables) {
            let output = base.join(output);
            let bytes = Table {
                axes: axes.clone(),
                states,
            }
            .to_msgpack()
            .with_context(|| shown(&output))?;
            std::fs::write(&output, bytes).with_context(|| shown(&output))?;
        }
    }
    Ok(())
}
