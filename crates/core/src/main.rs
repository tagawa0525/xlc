//! 均質な直方体の炉心の固有値計算を行うコマンド。
//!
//! 使い方: core <入力ファイル>
//!
//! 入力ファイル中の少数群定数のパスは、入力ファイルのあるディレクトリからの相対パス。

use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use anyhow::{Context, anyhow};
use xlc_core::input::CoreInput;
use xlc_core::nem::{Mesh, Options, Problem, solve};
use xlc_lattice::fewgroup::FewGroup;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.as_slice() {
        [input] => run(Path::new(input)),
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
    let solution = solve(&problem, &Options::default()).map_err(|e| anyhow!(e))?;
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
