//! 均質な直方体の炉心の固有値計算を行うコマンド。
//!
//! 使い方: core <入力ファイル>
//!
//! 入力ファイル中の少数群定数のパスは、入力ファイルのあるディレクトリからの相対パス。

use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use xlc_core::input::CoreInput;
use xlc_core::nem::{Mesh, Options, Problem, solve};
use xlc_lattice::fewgroup::FewGroup;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [input] = args.as_slice() else {
        eprintln!("usage: core <input>");
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
    let input = CoreInput::parse(&read(path)?).map_err(|e| format!("{}: {e}", path.display()))?;
    let fewgroup_path = path
        .parent()
        .unwrap_or(Path::new("."))
        .join(&input.fewgroup);
    let fewgroup = FewGroup::parse(&read(&fewgroup_path)?)
        .map_err(|e| format!("{}: {e}", fewgroup_path.display()))?;

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
    eprintln!(
        "{} iterations, {:.1} s",
        solution.iterations,
        start.elapsed().as_secs_f64()
    );
    Ok(())
}
