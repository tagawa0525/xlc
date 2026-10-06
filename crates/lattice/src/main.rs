//! ピンセルの格子計算を行い、少数群定数を書き出すコマンド。
//!
//! 使い方: lattice <入力ファイル>
//!
//! 入力ファイル中のライブラリと出力のパスは、入力ファイルのあるディレクトリからの相対パス。

use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use lattice::fewgroup::homogenize;
use lattice::geometry::PinCell;
use lattice::input::PinCellInput;
use lattice::material::Material;
use lattice::moc::{Options, solve};
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
        PinCellInput::parse(&read(path)?).map_err(|e| format!("{}: {e}", path.display()))?;
    let base = path.parent().unwrap_or(Path::new("."));
    let library_path = base.join(&input.library);
    let library = Library::parse(&read(&library_path)?)
        .map_err(|e| format!("{}: {e}", library_path.display()))?;

    let fuel = Material::mix(&library, &input.fuel, input.fuel_temperature)
        .map_err(|e| format!("fuel: {e}"))?;
    let moderator = Material::mix(&library, &input.moderator, input.moderator_temperature)
        .map_err(|e| format!("moderator: {e}"))?;
    let cell = PinCell::new(input.pitch, input.radii())?;
    let materials: Vec<&Material> = (0..cell.regions())
        .map(|i| {
            if i < input.fuel_rings {
                &fuel
            } else {
                &moderator
            }
        })
        .collect();
    let tracking = Tracking::new(&cell, input.azimuthal, input.spacing)?;
    let solution = solve(&tracking, &materials, &Options::default())?;
    let few = homogenize(
        &solution,
        &materials,
        &library.groups,
        &input.few_group_bounds,
    )?;

    let output = base.join(&input.output);
    std::fs::write(&output, few.to_text()).map_err(|e| format!("{}: {e}", output.display()))?;
    eprintln!(
        "k-inf {:.6} ({} iterations, {} tracks, {:.1} s) -> {}",
        solution.k,
        solution.iterations,
        tracking.tracks.len(),
        start.elapsed().as_secs_f64(),
        output.display()
    );
    Ok(())
}
