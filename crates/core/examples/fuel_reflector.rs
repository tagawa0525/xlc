//! 燃料と反射体の2つのノードを x に並べた炉心（ADR-0095 の I 型の格子と同じ体系）を、core の NEM で解く。
//! 炉心の装荷図の入力形式（#0030）が決まるまで、core のライブラリの API で問題を組む。
//!
//! 使い方: cargo run --release -p xlc-core --example fuel_reflector -- <燃料の少数群定数> <反射体の少数群定数> <幅 (cm)>
//!
//! -x、±y、±z は反射、+x は真空（反射体の少数群定数の右の辺の J/φ。ADR-0101）。k と、ノードごとの群の
//! 中性子束と核分裂率を標準出力に出す。

use std::path::Path;

use anyhow::{Context, anyhow};
use xlc_core::nem::{Boundary, Mesh, Options, Problem, solve};
use xlc_format::fewgroup::FewGroup;

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [fuel, reflector, width] = args.as_slice() else {
        return Err(anyhow!(
            "usage: fuel_reflector <fuel few-group> <reflector few-group> <width (cm)>"
        ));
    };
    let read = |path: &str| -> anyhow::Result<FewGroup> {
        let bytes = std::fs::read(Path::new(path)).with_context(|| path.to_string())?;
        FewGroup::from_msgpack(&bytes).with_context(|| path.to_string())
    };
    let h: f64 = width.parse().with_context(|| format!("width {width}"))?;
    let mut boundaries = [Boundary::Reflective; 6];
    boundaries[1] = Boundary::Vacuum;
    let problem = Problem {
        mesh: Mesh {
            dx: vec![h, h],
            dy: vec![h],
            dz: vec![h],
        },
        materials: vec![read(fuel)?, read(reflector)?],
        node_material: vec![0, 1],
        boundaries,
    };
    let solution = solve(&problem, &Options::default())?;
    println!("k-eff {:.6}", solution.k);
    for (n, flux) in solution.flux.iter().enumerate() {
        let m = &problem.materials[problem.node_material[n]];
        let fission: f64 = m.fission.iter().zip(flux).map(|(f, phi)| f * phi).sum();
        let groups: Vec<String> = flux.iter().map(|f| format!("{f:.9e}")).collect();
        println!("node {n} flux {} fission {fission:.9e}", groups.join(" "));
    }
    Ok(())
}
