//! 格子計算を行い、少数群定数とセルごとの核分裂率を書き出すコマンド。
//!
//! 使い方: lattice <入力ファイル>
//!
//! 入力ファイル中のライブラリと出力のパスは、入力ファイルのあるディレクトリからの相対パス。
//!
//! スレッド数は、環境変数 `RAYON_NUM_THREADS` があればその値、なければ物理コアの数とする。
//! MOC のスイープはメモリの転送速度で頭打ちになり、1つのコアを2つのスレッドで分ける
//! （SMT）とかえって遅くなるため。

use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use anyhow::{Context, anyhow};
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
        [input] => use_physical_cores().and_then(|()| run(Path::new(input))),
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

/// `RAYON_NUM_THREADS` がなければ、rayon のスレッド数を物理コアの数にする。
fn use_physical_cores() -> anyhow::Result<()> {
    if std::env::var_os("RAYON_NUM_THREADS").is_none() {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_cpus::get_physical())
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
        std::fs::write(&output, few.to_msgpack()).with_context(|| shown(&output))?;
    }
    if let Some(file) = &input.pin_powers {
        let rates = cell_fission_rates(&solution, &lattice, &materials)?;
        let file = base.join(file);
        std::fs::write(&file, to_msgpack(&rates, lattice.cells_per_side()))
            .with_context(|| shown(&file))?;
    }
    eprintln!(
        "k {:.6} ({} iterations, {} tracks, {:.1} s)",
        solution.k,
        solution.iterations,
        tracking.tracks.len(),
        start.elapsed().as_secs_f64(),
    );
    Ok(())
}
