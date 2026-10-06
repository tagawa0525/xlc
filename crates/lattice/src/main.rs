//! 格子計算を行い、少数群定数とセルごとの核分裂率を書き出すコマンド。
//!
//! 使い方: lattice <入力ファイル>
//!
//! 入力ファイル中のライブラリと出力のパスは、入力ファイルのあるディレクトリからの相対パス。

use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use lattice::fewgroup::homogenize;
use lattice::input::{CrossSections, LatticeInput};
use lattice::macroscopic;
use lattice::material::Material;
use lattice::moc::{Options, solve};
use lattice::power::{cell_fission_rates, to_text};
use lattice::tracking::Tracking;
use xs::library::Library;

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
    let start = Instant::now();
    let read = |p: &Path| std::fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()));
    let input =
        LatticeInput::parse(&read(path)?).map_err(|e| format!("{}: {e}", path.display()))?;
    let base = path.parent().unwrap_or(Path::new("."));
    let (lattice, region_materials) = input.geometry()?;

    let (named, groups) = match &input.cross_sections {
        CrossSections::Library(file) => {
            let file = base.join(file);
            let library =
                Library::parse(&read(&file)?).map_err(|e| format!("{}: {e}", file.display()))?;
            let named = input
                .materials
                .iter()
                .map(|m| {
                    Material::mix(&library, &m.composition, m.temperature)
                        .map(|mixed| (m.name.clone(), mixed))
                        .map_err(|e| format!("material `{}`: {e}", m.name))
                })
                .collect::<Result<Vec<_>, String>>()?;
            (named, Some(library.groups))
        }
        CrossSections::Macroscopic(file) => {
            let file = base.join(file);
            let named = macroscopic::parse(&read(&file)?)
                .map_err(|e| format!("{}: {e}", file.display()))?;
            (named, None)
        }
    };
    let materials = region_materials
        .iter()
        .map(|name| {
            named
                .iter()
                .find(|(n, _)| n == name)
                .map(|(_, m)| m)
                .ok_or_else(|| format!("material `{name}` is not in the cross sections"))
        })
        .collect::<Result<Vec<&Material>, String>>()?;

    let tracking =
        Tracking::with_boundaries(&lattice, input.azimuthal, input.spacing, input.boundaries)?;
    let solution = solve(&tracking, &materials, &Options::default())?;
    if let Some((output, bounds)) = &input.output {
        let groups = groups
            .as_ref()
            .expect("the input checks that output needs a library");
        let few = homogenize(&solution, &materials, groups, bounds)?;
        let output = base.join(output);
        std::fs::write(&output, few.to_text()).map_err(|e| format!("{}: {e}", output.display()))?;
    }
    if let Some(file) = &input.pin_powers {
        let rates = cell_fission_rates(&solution, &lattice, &materials)?;
        let file = base.join(file);
        std::fs::write(&file, to_text(&rates, lattice.cells_per_side()))
            .map_err(|e| format!("{}: {e}", file.display()))?;
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
