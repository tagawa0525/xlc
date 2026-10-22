//! 燃焼の履歴の燃焼点からの分岐（ADR-0187）。
//!
//! 入力の `[history]` が挙げた燃焼点ごとに、燃焼の履歴（ADR-0173）の燃焼の領域の数密度を固定し、`[branches]` の状態
//! （ADR-0097）で輸送する。分岐の状態では燃焼させない。燃焼の領域は燃焼の輸送と同じ `<物質>@<セル>` の物質として解き、
//! 輸送する核種は履歴の `transportable`（数密度が 0 の核種も含め、負の値は 0 にする）。状態は燃焼度がいちばん遅く変わり、
//! 各燃焼度の点の中は分岐の軸の直積の順。
//!
//! ```toml
//! [history]
//! file = "history.msgpack"   # 燃焼の履歴（入力のファイルのあるディレクトリから）
//! burnup = [0, 10, 20, 40]   # 分岐を解く燃焼点 (MWd/kg)。履歴の燃焼度から選ぶ
//! ```

use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{anyhow, bail};
use xlc_format::library::Library;
use xlc_format::table::Axis;

use crate::depletion::{
    History, check_negative, check_region_names, depletion_regions, has_data_at, region_material,
    solved_names,
};
use crate::geometry::Lattice;
use crate::input::{LatticeInput, MaterialInput};
use crate::symmetry::Chosen;

/// 入力の `[history]`（ADR-0187）。
#[derive(Debug, Clone, PartialEq)]
pub struct HistoryInput {
    /// 燃焼の履歴のファイル（入力のファイルのあるディレクトリから）
    pub file: PathBuf,
    /// 分岐を解く燃焼点 (MWd/kg)。0 以上で増え、どれも履歴の燃焼度のどれかと等しい
    pub burnup: Vec<f64>,
}

/// 履歴の燃焼点からの分岐の状態。
#[derive(Debug, Clone, PartialEq)]
pub struct HistoryBranches {
    /// テーブルの軸（燃焼度、続いて分岐の軸）
    pub axes: Vec<Axis>,
    /// 軸の直積の順（燃焼度がいちばん遅く変わる）の、解く物質の組
    pub states: Vec<Vec<MaterialInput>>,
    /// 領域ごとの解く物質の名前（燃焼の領域は `<物質>@<セル>`）
    pub region_materials: Vec<String>,
}

/// 入力 `input` の `[history]` と `[branches]` から、履歴 `history` の燃焼点ごとの分岐の状態を作る。`region_materials` は
/// 入力の領域ごとの物質の名前、`chosen` はその名前で選んだ対称操作。入力の格子が履歴の格子と違えばエラー
pub fn history_branches(
    input: &LatticeInput,
    history: &History,
    library: &Library,
    lattice: &Lattice,
    region_materials: &[String],
    chosen: &Chosen,
) -> anyhow::Result<HistoryBranches> {
    let selection = input
        .history
        .as_ref()
        .ok_or_else(|| anyhow!("the input has no [history]"))?;
    let branches = input
        .branches
        .as_ref()
        .ok_or_else(|| anyhow!("[history] needs [branches] (ADR-0187)"))?;
    // 燃焼させた物質（履歴の燃焼の領域の物質。現れる順）
    let mut depleted: Vec<String> = Vec::new();
    for d in &history.regions {
        if !depleted.contains(&d.material) {
            depleted.push(d.material.clone());
        }
    }
    // 入力から求めた燃焼の領域が、履歴の燃焼の領域（物質、組の代表のセル、体積）と同じこと。MOC の領域の番号は入力のもの
    // を使う（案内管の中身を変えた格子でも、燃料のセルが同じなら同じ領域になる）
    let representatives = chosen.shielding().cell_representatives();
    let regions = depletion_regions(lattice, region_materials, &representatives, &depleted);
    let key =
        |material: &str, cell: usize, volume: f64| (material.to_string(), cell, volume.to_bits());
    let found: Vec<_> = regions
        .iter()
        .map(|d| key(&d.material, d.cell, d.volume))
        .collect();
    let recorded: Vec<_> = history
        .regions
        .iter()
        .map(|d| key(&d.material, d.cell, d.volume))
        .collect();
    if found != recorded {
        bail!(
            "the depletion regions of the lattice (material, cell, volume) differ from those of the history: \
             {} in the lattice, {} in the history; give the lattice of the history (ADR-0187)",
            regions.len(),
            history.regions.len()
        );
    }
    let index: HashMap<&str, usize> = history
        .nuclides
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    // 燃焼させた物質の入力の組成が、履歴の初めの数密度と同じこと
    let first = history
        .points
        .first()
        .filter(|p| p.burnup == 0.0)
        .ok_or_else(|| anyhow!("the history does not start at burnup 0"))?;
    for (d, densities) in history.regions.iter().zip(&first.densities) {
        let material = input
            .materials
            .iter()
            .find(|m| m.name == d.material)
            .ok_or_else(|| anyhow!("depleted material `{}` is not in the input", d.material))?;
        let mut expected = vec![0.0; history.nuclides.len()];
        for (name, density) in &material.composition {
            let i = index.get(name.as_str()).ok_or_else(|| {
                anyhow!(
                    "material `{}`: {name} is not a nuclide of the history",
                    d.material
                )
            })?;
            expected[*i] = *density;
        }
        if &expected != densities {
            bail!(
                "material `{}`: the composition differs from the initial densities of the history; \
                 give the lattice of the history (ADR-0187)",
                d.material
            );
        }
    }
    let points = selection
        .burnup
        .iter()
        .map(|&b| {
            history
                .points
                .iter()
                .position(|p| p.burnup == b)
                .ok_or_else(|| anyhow!("burnup {b} MWd/kg is not a point of the history"))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    let transported = history
        .transportable
        .iter()
        .map(|name| {
            index
                .get(name.as_str())
                .copied()
                .ok_or_else(|| anyhow!("the transportable {name} is not a nuclide of the history"))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    check_region_names(&input.materials, &history.regions)?;
    let renamed = solved_names(lattice, region_materials, &representatives, &depleted);
    for &p in &points {
        check_negative(&history.regions, &history.points[p].densities, |i| {
            history.nuclides[i].clone()
        })?;
    }
    // 燃焼の領域の物質の、分岐の状態ごとの温度。輸送する核種がその温度のデータを持つこと（近い温度で代用しない。ADR-0173 の案 N2）
    let temperatures = branches
        .states()
        .iter()
        .map(|state| {
            history
                .regions
                .iter()
                .map(|d| {
                    let temperature = state
                        .iter()
                        .find(|m| m.name == d.material)
                        .ok_or_else(|| {
                            anyhow!("depleted material `{}` is not in the state", d.material)
                        })?
                        .temperature;
                    if let Some(&i) = transported
                        .iter()
                        .find(|&&i| !has_data_at(library, &history.nuclides[i], temperature))
                    {
                        bail!(
                            "{} (transported in the history) has no library data at {temperature} K of material `{}`",
                            history.nuclides[i],
                            d.material
                        );
                    }
                    Ok(temperature)
                })
                .collect::<anyhow::Result<Vec<_>>>()
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    let mut states = Vec::with_capacity(points.len() * branches.states().len());
    for &p in &points {
        for (state, temperatures) in branches.states().iter().zip(&temperatures) {
            let mut materials: Vec<MaterialInput> = state
                .iter()
                .filter(|m| !depleted.contains(&m.name))
                .cloned()
                .collect();
            for ((d, densities), &temperature) in history
                .regions
                .iter()
                .zip(&history.points[p].densities)
                .zip(temperatures)
            {
                materials.push(region_material(
                    d,
                    temperature,
                    densities,
                    &transported,
                    |i| history.nuclides[i].clone(),
                ));
            }
            states.push(materials);
        }
    }
    let mut axes = vec![Axis::Burnup {
        points: selection.burnup.clone(),
    }];
    axes.extend(branches.axes().iter().cloned());
    Ok(HistoryBranches {
        axes,
        states,
        region_materials: renamed,
    })
}

#[cfg(test)]
mod tests;
