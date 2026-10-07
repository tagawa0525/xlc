//! 格子の入力の物質の多群定数（無限希釈）を書く（OpenMC の多群モードで同じ定数を解いて、MOC の輸送を検証する用）。
//!
//! 使い方: mg_materials <格子の入力 (TOML)> <出力 (msgpack)>
//!
//! 物質ごとに、全断面積・核分裂・νΣf の行列 [入射の群][出射の群]（核種ごとの νΣf × χ の和）・散乱行列の
//! ルジャンドルのモーメント [l][入射の群][出射の群] (1/cm) を書く。群は高いエネルギーから。自己遮蔽はしない
//! （検証するのは輸送の解き方なので、どちらも同じ定数であればよい）。

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, anyhow};
use serde::Serialize;
use xlc_lattice::input::{CrossSections, LatticeInput};
use xlc_lattice::material::Material;
use xlc_xs::library::Library;
use xlc_xs::msgpack::to_named;

#[derive(Serialize)]
struct Constants {
    total: Vec<f64>,
    fission: Vec<f64>,
    /// [入射の群][出射の群]
    nu_fission: Vec<Vec<f64>>,
    /// [l][入射の群][出射の群]
    scatter: Vec<Vec<Vec<f64>>>,
}

#[derive(Serialize)]
struct Materials {
    format: &'static str,
    /// 群境界 (eV)、降順
    group_bounds: Vec<f64>,
    materials: BTreeMap<String, Constants>,
}

fn constants(m: &Material) -> Constants {
    let ng = m.groups();
    let mut nu_fission = vec![vec![0.0; ng]; ng];
    for term in &m.fission_terms {
        for (row, n) in nu_fission.iter_mut().zip(&term.nu_fission) {
            for (value, chi) in row.iter_mut().zip(&term.chi) {
                *value += n * chi;
            }
        }
    }
    Constants {
        total: m.total.clone(),
        fission: m.fission.clone(),
        nu_fission,
        scatter: std::iter::once(&m.scatter)
            .chain(&m.anisotropic)
            .cloned()
            .collect(),
    }
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let [_, input, output] = args.as_slice() else {
        return Err(anyhow!("usage: mg_materials <lattice input> <output>"));
    };
    let path = Path::new(input);
    let input = LatticeInput::parse(&std::fs::read_to_string(path)?).context("input")?;
    let CrossSections::Library(file) = &input.cross_sections else {
        return Err(anyhow!("the input must use a library"));
    };
    let base = path.parent().unwrap_or(Path::new("."));
    let library = Library::from_msgpack(&std::fs::read(base.join(file))?)?;
    let materials = input
        .materials
        .iter()
        .map(|m| {
            let mixed = Material::mix(&library, &m.composition, m.temperature)
                .with_context(|| format!("material `{}`", m.name))?;
            Ok((m.name.clone(), constants(&mixed)))
        })
        .collect::<anyhow::Result<_>>()?;
    let bytes = to_named(&Materials {
        format: "xlc-mg-materials",
        group_bounds: library.groups.bounds().to_vec(),
        materials,
    })?;
    std::fs::write(output, bytes)?;
    Ok(())
}
