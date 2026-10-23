//! 格子計算を行い、少数群定数のテーブル、セルごとの核分裂率、細かい群の媒質（ADR-0128）を書き出すコマンド。
//! 入力に分岐（`[branches]`）があれば、状態ごとに格子計算を解き、軸の直積の順の状態をテーブルにする（ADR-0097）。
//! 入力に燃焼（`[depletion]`）があれば、燃焼点ごとの k∞、数密度、ピン出力、少数群定数を燃焼の履歴のファイルに書く（ADR-0173、ADR-0181）。
//! 入力に燃焼の履歴（`[history]`）があれば、選んだ燃焼点の組成で分岐を解き、燃焼度の軸を持つテーブルにする（ADR-0187）。
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
use xlc_depletion::chain::Chain;
use xlc_format::fewgroup::FewGroup;
use xlc_format::library::Library;
use xlc_format::medium::Medium;
use xlc_lattice::axial::{axial, check_nodes};
use xlc_lattice::branches::{Starts, state_indices, tables};
use xlc_lattice::depletion::{History, LibraryReactions};
use xlc_lattice::fewgroup::{Node, homogenize, homogenize_nodes, pin_shapes, split_constants};
use xlc_lattice::input::SymmetryChoice;
use xlc_lattice::input::{
    CrossSections, FewGroupFiles, FewGroupOutput, LatticeInput, MaterialInput, VOID,
};
use xlc_lattice::macroscopic;
use xlc_lattice::material::Material;
use xlc_lattice::medium::{check_medium, check_orders, homogenize_medium};
use xlc_lattice::moc::{Options, available_order, solve_symmetric_from};
use xlc_lattice::power::{cell_fission_rates, to_msgpack};
use xlc_lattice::shielding::{ShieldingReuse, shield_fuel_reusing};
use xlc_lattice::symmetry::choose;
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

/// `RAYON_NUM_THREADS` がなければ、rayon のスレッド数を物理コアの数にする。
/// 物理コアの数はマシン全体の値なので、このプロセスに許された CPU の数（論理）で抑える。
#[allow(
    clippy::disallowed_methods,
    reason = "RAYON_NUM_THREADS は並列の数だけを選ぶ。出力が数によらないかは未確認（#0060、ADR-0137）"
)]
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
    // 輸送の外の段の時間は debug の区間で測る（FmtSpan::CLOSE で、閉じるときに time.busy を出す）
    let reading = tracing::debug_span!("input").entered();
    let shown = |p: &Path| p.display().to_string();
    let text = std::fs::read_to_string(path).with_context(|| shown(path))?;
    let input = LatticeInput::parse(&text).with_context(|| shown(path))?;
    let base = path.parent().unwrap_or(Path::new("."));
    let (lattice, region_materials) = input.geometry().with_context(|| shown(path))?;
    // 軸方向の定数のノードの並びは、MOC を解く前に確かめる（ADR-0128）
    if let Some(FewGroupOutput {
        files: FewGroupFiles::Nodes(entries),
        axial: true,
        ..
    }) = &input.output
    {
        let nodes: Vec<Node> = entries.iter().map(|(_, node)| node.clone()).collect();
        check_nodes(lattice.rows(), &nodes).map_err(|e| anyhow!("{}: {e}", shown(path)))?;
    }

    drop(reading);
    let tracking = tracing::debug_span!("tracking").in_scope(|| {
        Tracking::with_boundaries(&lattice, input.azimuthal, input.spacing, input.boundaries)
    })?;
    // 成り立つ対称操作（ADR-0112）。形状、トラック、領域ごとの入力の物質の名前は分岐の状態によらないので、解く前に一度だけ
    // 求め、自己遮蔽のピンと塊の組（ADR-0112）と MOC に再利用する（使った対称性と、外した操作の理由をログに出す）。強制した
    // 対称性は、その結果から取り出し、成り立たない操作があれば、自己遮蔽の前にその操作の名前と理由を示して止める。"none" は
    // 検証の基準なので判定せず、自己遮蔽も MOC も対称性なしで解く
    let forced = match input.symmetry {
        SymmetryChoice::Detect => None,
        SymmetryChoice::Force(symmetry) => Some(symmetry),
    };
    let chosen = tracing::debug_span!("symmetry")
        .in_scope(|| choose(&tracking, &region_materials, forced))
        .with_context(|| format!("{}: symmetry", shown(path)))?;
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
    // 断面積の出どころは一度だけ読む。巨視的断面積は分岐を持たない（入力の読み込みで確かめてある）
    enum Source {
        Library(Library),
        Macroscopic(Vec<(String, Material)>),
    }
    let loading = tracing::debug_span!("library").entered();
    let source = match &input.cross_sections {
        CrossSections::Library(file) => {
            let file = base.join(file);
            let bytes = std::fs::read(&file).with_context(|| shown(&file))?;
            Source::Library(Library::from_msgpack(&bytes).with_context(|| shown(&file))?)
        }
        CrossSections::Macroscopic(file) => {
            let file = base.join(file);
            let text = std::fs::read_to_string(&file).with_context(|| shown(&file))?;
            Source::Macroscopic(macroscopic::parse(&text).with_context(|| shown(&file))?)
        }
    };
    drop(loading);
    if let Some(depletion) = &input.depletion {
        let Source::Library(library) = &source else {
            return Err(anyhow!("{}: depletion needs a library", shown(path)));
        };
        let file = base.join(&depletion.chain);
        let bytes = std::fs::read(&file).with_context(|| shown(&file))?;
        let chain = Chain::from_msgpack(&bytes).with_context(|| shown(&file))?;
        let history = xlc_lattice::depletion::run(
            &input,
            depletion,
            library,
            &chain,
            &tracking,
            &region_materials,
            &chosen,
            &LibraryReactions(library),
        )
        .with_context(|| format!("{}: depletion", shown(path)))?;
        for point in &history.points {
            println!("burnup {} k {:.6}", point.burnup, point.k);
        }
        // 時間の内訳（ADR-0173 の「全体像」を測るため）
        let (t, d) = (history.transport, history.depletion);
        tracing::info!(
            solves = history.solves,
            regions = history.regions.len(),
            shielding = t.shielding.as_secs_f64(),
            moc = t.moc.as_secs_f64(),
            rates = t.rates.as_secs_f64(),
            matrix = d.matrix.as_secs_f64(),
            cram = d.cram.as_secs_f64(),
            total_seconds = start.elapsed().as_secs_f64(),
            "depletion"
        );
        let output = base.join(&depletion.output);
        let bytes = history.to_msgpack().with_context(|| shown(&output))?;
        std::fs::write(&output, bytes).with_context(|| shown(&output))?;
        return Ok(());
    }
    let groups = match &source {
        Source::Library(library) => Some(library.groups.clone()),
        Source::Macroscopic(_) => None,
    };
    // 燃焼の履歴の燃焼点からの分岐（ADR-0187）: 解く状態、燃焼の領域の物質の名前、燃焼度の軸を足した軸
    let history_branches = match (&input.history, &source) {
        (Some(history), Source::Library(library)) => {
            let file = base.join(&history.file);
            let bytes = std::fs::read(&file).with_context(|| shown(&file))?;
            let read = History::from_msgpack(&bytes).with_context(|| shown(&file))?;
            Some(
                xlc_lattice::history::history_branches(
                    &input,
                    &read,
                    library,
                    &lattice,
                    &region_materials,
                    &chosen,
                )
                .with_context(|| format!("{}: history", shown(path)))?,
            )
        }
        (Some(_), Source::Macroscopic(_)) => {
            return Err(anyhow!("{}: [history] needs a library", shown(path)));
        }
        (None, _) => None,
    };
    // 状態の軸の値（エラーとログに付ける）
    let axes = match &history_branches {
        Some(h) => &h.axes[..],
        None => input.branches.as_ref().map_or(&[][..], |b| b.axes()),
    };
    // 自己遮蔽と MOC に渡す領域ごとの物質の名前（燃焼の領域は `<物質>@<セル>`）
    let solved_materials = history_branches
        .as_ref()
        .map_or(&region_materials[..], |h| &h.region_materials[..]);
    // 平板の物質の媒質（ADR-0128）。ライブラリの群境界と分岐の軸に合うことを解く前に確かめる
    let media = input
        .media
        .iter()
        .map(|m| {
            let file = base.join(&m.file);
            let bytes = std::fs::read(&file).with_context(|| shown(&file))?;
            let medium = Medium::from_msgpack(&bytes).with_context(|| shown(&file))?;
            let bounds = groups.as_ref().map_or(&[][..], |g| g.bounds());
            check_medium(&medium, bounds, axes)
                .map_err(|e| anyhow!("{}: material `{}`: {e}", shown(&file), m.name))?;
            Ok(medium)
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    let medium_names: Vec<String> = input.media.iter().map(|m| m.name.clone()).collect();
    let counts: Vec<usize> = axes.iter().map(|a| a.points().len()).collect();
    let values = |state: usize| -> String {
        axes.iter()
            .zip(state_indices(&counts, state))
            .filter_map(|(a, i)| a.points().get(i).map(|p| format!("{} {p}", a.name())))
            .collect::<Vec<_>>()
            .join(", ")
    };

    // 分岐の状態ごとに格子計算を解き、出力ごとの少数群定数を状態の順に集める（ADR-0097）
    let plain = input.branch_states();
    let count = history_branches.as_ref().map_or(plain.len(), |h| h.len());
    let mut per_state: Vec<Vec<Vec<FewGroup>>> = Vec::with_capacity(count);
    // 状態ごとの細かい群の媒質の出力（ADR-0128）
    let mut written = Vec::new();
    // 燃料でない物質が直前の状態と同じなら、燃料のピンの Dancoff 係数を持ち越す（ADR-0194）
    let mut reuse = ShieldingReuse::new(&tracking);
    // 状態ごとに、前に解いた近い状態の MOC の収束した状態から反復を始める（ADR-0195）
    let mut starts = Starts::new(&counts);
    for state in 0..count {
        // 履歴からの分岐の物質は状態ごとに作る（ADR-0187）
        let built;
        let state_materials: &[MaterialInput] = match &history_branches {
            Some(h) => {
                built = h.state(state);
                &built
            }
            None => &plain[state],
        };
        let state_start = Instant::now();
        let context = || format!("{}: state {state} ({})", shown(path), values(state));
        let (unique, region_index) = match &source {
            Source::Library(library) => {
                if input.self_shielding {
                    let shielded = shield_fuel_reusing(
                        library,
                        &mut reuse,
                        state_materials,
                        solved_materials,
                        chosen.shielding(),
                        &input.shielding_options(),
                    )
                    .with_context(|| format!("{}: self-shielding", context()))?;
                    (shielded.materials, shielded.region_material)
                } else {
                    let mut named = state_materials
                        .iter()
                        .map(|m| {
                            Material::mix(library, &m.composition, m.temperature)
                                .map(|mixed| (m.name.clone(), mixed))
                                .map_err(|e| anyhow!("material `{}`: {e}", m.name))
                        })
                        .collect::<anyhow::Result<Vec<_>>>()
                        .with_context(context)?;
                    // 混ぜた物質は、この状態の成分から作る（真空の成分は核種を持たない。ADR-0128）
                    for mixture in &input.mixtures {
                        let parts: Vec<(&_, f64)> = mixture
                            .parts
                            .iter()
                            .filter(|(part, _)| part != VOID)
                            .map(|(part, f)| {
                                state_materials
                                    .iter()
                                    .find(|m| &m.name == part)
                                    .map(|m| (m, *f))
                                    .ok_or_else(|| {
                                        anyhow!(
                                            "material `{}`: component `{part}` is not in the state",
                                            mixture.name
                                        )
                                    })
                            })
                            .collect::<anyhow::Result<_>>()
                            .with_context(context)?;
                        let mixed = Material::mixture(library, &parts)
                            .map_err(|e| anyhow!("material `{}`: {e}", mixture.name))
                            .with_context(context)?;
                        named.push((mixture.name.clone(), mixed));
                    }
                    // 媒質は分岐の状態と同じ番号の状態（軸が同じことは確かめてある）
                    for (m, medium) in input.media.iter().zip(&media) {
                        let from = medium
                            .states
                            .get(state)
                            .ok_or_else(|| anyhow!("material `{}`: no state {state}", m.name))
                            .and_then(|s| {
                                Material::from_medium(s)
                                    .map_err(|e| anyhow!("material `{}`: {e}", m.name))
                            })
                            .with_context(context)?;
                        named.push((m.name.clone(), from));
                    }
                    check_orders(&named, &medium_names)
                        .map_err(|e| anyhow!(e))
                        .with_context(context)?;
                    by_name(named, library.groups.len())?
                }
            }
            Source::Macroscopic(named) => {
                let groups = named.first().map_or(0, |(_, m)| m.groups());
                by_name(named.clone(), groups)?
            }
        };
        let materials: Vec<&Material> = region_index.iter().map(|&i| &unique[i]).collect();
        // 媒質（ADR-0128）の領域。領域がすべて1つの核分裂のある媒質のノードの拡散係数は、その媒質の B1 の D（ADR-0141）
        let media_regions: Vec<Option<usize>> = region_materials
            .iter()
            .map(|name| input.media.iter().position(|m| &m.name == name))
            .collect();

        let options = Options {
            scattering_order: input
                .scattering_order
                .unwrap_or_else(|| available_order(&materials)),
            ..Options::default()
        };
        let (solution, converged) = solve_symmetric_from(
            &tracking,
            &materials,
            &options,
            chosen.moc(),
            starts.start(),
        )
        .with_context(context)?;
        starts.solved(converged);
        let emitting = tracing::debug_span!("output").entered();
        if let Some(FewGroupOutput {
            bounds,
            files,
            pins,
            split,
            axial: along_z,
        }) = &input.output
        {
            let groups = groups
                .as_ref()
                .ok_or_else(|| anyhow!("few-group output needs a library"))?;
            let boundaries = tracking.boundaries();
            let fews = match files {
                FewGroupFiles::Lattice(_) => {
                    let mut few = homogenize(&solution, &materials, groups, bounds, &boundaries)?;
                    // ピン出力の再構成の形状関数と角の不連続因子（ADR-0096）
                    if *pins {
                        few.pins = Some(
                            pin_shapes(
                                &solution,
                                &lattice,
                                &materials,
                                groups,
                                bounds,
                                &boundaries,
                            )
                            .with_context(context)?,
                        );
                    }
                    // 集合体を分けた部分ノードの定数（ADR-0108）
                    if *split > 1 {
                        vec![
                            split_constants(
                                &few,
                                &solution,
                                &lattice,
                                &materials,
                                groups,
                                bounds,
                                &boundaries,
                                *split,
                            )
                            .with_context(context)?,
                        ]
                    } else {
                        vec![vec![few]]
                    }
                }
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
                        &media_regions,
                        *pins,
                    )?
                    .into_iter()
                    .map(|few| vec![few])
                    .collect()
                }
            };
            // 格子の x を炉心の z にした定数。z の厚さはノードの x の幅（ADR-0128）
            let fews = if *along_z {
                let widths = lattice.widths();
                let node_widths: Vec<f64> = match files {
                    FewGroupFiles::Lattice(_) => vec![widths.iter().sum()],
                    FewGroupFiles::Nodes(entries) => entries
                        .iter()
                        .map(|(_, node)| widths[node.columns.clone()].iter().sum())
                        .collect(),
                };
                fews.iter()
                    .zip(node_widths)
                    .map(|(subs, width)| {
                        subs.iter()
                            .map(|few| axial(few, width))
                            .collect::<Result<Vec<_>, _>>()
                    })
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|e| anyhow!(e))
                    .with_context(context)?
            } else {
                fews
            };
            per_state.push(fews);
        }
        if input.medium.is_some() {
            written.push(
                homogenize_medium(&solution, &materials, &tracking.boundaries())
                    .with_context(context)?,
            );
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
        drop(emitting);
        tracing::info!(
            state,
            iterations = solution.iterations,
            tracks = tracking.tracks.len(),
            seconds = state_start.elapsed().as_secs_f64(),
            total_seconds = start.elapsed().as_secs_f64(),
            "lattice"
        );
    }

    let _write = tracing::debug_span!("write").entered();
    if let Some(output) = &input.medium {
        let groups = groups
            .as_ref()
            .ok_or_else(|| anyhow!("a medium needs a library"))?;
        let medium = Medium {
            bounds: groups.bounds().to_vec(),
            axes: axes.to_vec(),
            states: written,
        };
        let output = base.join(output);
        let bytes = medium.to_msgpack().with_context(|| shown(&output))?;
        std::fs::write(&output, bytes).with_context(|| shown(&output))?;
    }
    if let Some(FewGroupOutput { files, split, .. }) = &input.output {
        let outputs: Vec<&Path> = match files {
            FewGroupFiles::Lattice(output) => vec![output],
            FewGroupFiles::Nodes(entries) => entries.iter().map(|(o, _)| o.as_path()).collect(),
        };
        let mut tables =
            tables(axes, *split, per_state).map_err(|e| anyhow!("{}: {e}", shown(path)))?;
        // 履歴からの分岐は格子全体の 1 つのテーブル（入力の検査）で、初めの重金属の質量を持つ（ADR-0190）
        if let Some(h) = &history_branches {
            for table in &mut tables {
                table.heavy_metal = Some(h.heavy_metal);
            }
        }
        for (output, table) in outputs.into_iter().zip(tables) {
            let output = base.join(output);
            let bytes = table.to_msgpack().with_context(|| shown(&output))?;
            std::fs::write(&output, bytes).with_context(|| shown(&output))?;
        }
    }
    Ok(())
}
