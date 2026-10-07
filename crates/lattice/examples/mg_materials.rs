//! 格子の入力の物質の多群定数（無限希釈）を書く（OpenMC の多群モードで同じ定数を解いて、MOC の輸送を検証する用）。
//!
//! 使い方: mg_materials <格子の入力 (TOML)> <出力 (msgpack)>
//!
//! 物質ごとに、全断面積・核分裂・νΣf の行列 [入射の群][出射の群]（核種ごとの νΣf × χ の和）・散乱行列の
//! ルジャンドルのモーメント [l][入射の群][出射の群] (1/cm) を書く。群は高いエネルギーから。モーメントは入力の
//! `scattering-order` まで（省略するとライブラリにある全て）で、xlc が解く定数と同じにする。入力は
//! `self-shielding = false` であること（自己遮蔽した定数はリングごとに違い、物質ごとには書けない）。

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, anyhow};
use serde::Serialize;
use xlc_lattice::input::{CrossSections, LatticeInput};
use xlc_lattice::material::Material;
use xlc_lattice::moc::available_order;
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
    /// 物質ごとの温度 (K)（解く問題と合っているかを確かめるため）
    temperatures: BTreeMap<String, f64>,
    materials: BTreeMap<String, Constants>,
}

/// 散乱のモーメントを P_`order` まで持つ物質の定数（物質にその次数までなければ `None`）
fn constants(m: &Material, order: usize) -> Option<Constants> {
    Some(Constants {
        total: m.total.clone(),
        fission: m.fission.clone(),
        nu_fission: m.fission_matrix(),
        scatter: std::iter::once(&m.scatter)
            .chain(m.anisotropic.get(..order)?)
            .cloned()
            .collect(),
    })
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let [_, input, output] = args.as_slice() else {
        return Err(anyhow!("usage: mg_materials <lattice input> <output>"));
    };
    let path = Path::new(input);
    let shown = |p: &Path| p.display().to_string();
    let text = std::fs::read_to_string(path).with_context(|| shown(path))?;
    let input = LatticeInput::parse(&text).with_context(|| shown(path))?;
    let CrossSections::Library(file) = &input.cross_sections else {
        return Err(anyhow!("the input must use a library"));
    };
    if input.self_shielding {
        return Err(anyhow!(
            "the input must set `self-shielding = false` (shielded constants differ by ring)"
        ));
    }
    let file = path.parent().unwrap_or(Path::new(".")).join(file);
    let bytes = std::fs::read(&file).with_context(|| shown(&file))?;
    let library = Library::from_msgpack(&bytes).with_context(|| shown(&file))?;
    let mixed = input
        .materials
        .iter()
        .map(|m| {
            let mixed = Material::mix(&library, &m.composition, m.temperature)
                .with_context(|| format!("material `{}`", m.name))?;
            Ok((m.name.clone(), mixed))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    let references: Vec<&Material> = mixed.iter().map(|(_, m)| m).collect();
    let available = available_order(&references);
    let order = input.scattering_order.unwrap_or(available);
    if order > available {
        return Err(anyhow!(
            "scattering order {order} exceeds P{available} of the library"
        ));
    }
    let materials = mixed
        .iter()
        .map(|(name, m)| {
            let c = constants(m, order)
                .ok_or_else(|| anyhow!("material `{name}` has no P{order} scattering"))?;
            Ok((name.clone(), c))
        })
        .collect::<anyhow::Result<_>>()?;
    let temperatures = input
        .materials
        .iter()
        .map(|m| (m.name.clone(), m.temperature))
        .collect();
    let bytes = to_named(&Materials {
        format: "xlc-mg-materials",
        group_bounds: library.groups.bounds().to_vec(),
        temperatures,
        materials,
    })?;
    std::fs::write(output, bytes).with_context(|| output.clone())?;
    Ok(())
}
