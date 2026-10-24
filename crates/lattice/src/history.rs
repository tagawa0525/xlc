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
//!
//! 状態の自己遮蔽は `[branches]` の `shielding` で選ぶ（[`crate::branches::BranchShielder`]。ADR-0201）。

use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{anyhow, bail};
use xlc_format::library::Library;
use xlc_format::table::Axis;

use crate::depletion::{
    DepletionRegion, HEAVY_METAL_Z, History, check_negative, check_region_names, depletion_regions,
    has_data_at, heavy_metal_mass, region_material, solved_names,
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

/// 履歴の燃焼点からの分岐の状態。状態ごとの物質は `state` で作る。
#[derive(Debug, Clone, PartialEq)]
pub struct HistoryBranches {
    /// テーブルの軸（燃焼度、続いて分岐の軸）
    pub axes: Vec<Axis>,
    /// 領域ごとの解く物質の名前（燃焼の領域は `<物質>@<セル>`）
    pub region_materials: Vec<String>,
    /// 初めの重金属の質量 (g/cm、格子全体の高さ 1 cm あたり)。テーブルの `heavy_metal`（ADR-0190）
    pub heavy_metal: f64,
    /// 履歴の燃焼の領域
    regions: Vec<DepletionRegion>,
    /// 選んだ燃焼度の点ごとの、燃焼の領域の数密度（[点][領域][履歴の核種]）
    densities: Vec<Vec<Vec<f64>>>,
    /// 履歴の核種の名前
    nuclides: Vec<String>,
    /// 輸送する核種（履歴の核種の番号）
    transported: Vec<usize>,
    /// 分岐の状態ごとの、燃焼させない物質と、燃焼の領域ごとの温度 (K)
    branch: Vec<(Vec<MaterialInput>, Vec<f64>)>,
    /// 分岐の状態ごとの参照の番号（`Branches::reference`。ADR-0201）
    reference_of: Vec<Option<usize>>,
    /// 参照の状態（ホウ素が公称。`Branches::references` の順）の、燃焼させない物質と、燃焼の領域ごとの温度 (K)
    references: Vec<(Vec<MaterialInput>, Vec<f64>)>,
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
        let at = found.iter().zip(&recorded).position(|(a, b)| a != b);
        let detail = match at {
            Some(k) => format!(
                "region {k}: {:?} in the lattice, {:?} in the history (material, cell, volume)",
                (&regions[k].material, regions[k].cell, regions[k].volume),
                (
                    &history.regions[k].material,
                    history.regions[k].cell,
                    history.regions[k].volume
                )
            ),
            None => format!(
                "{} regions in the lattice, {} in the history",
                regions.len(),
                history.regions.len()
            ),
        };
        bail!(
            "the depletion regions of the lattice differ from those of the history: {detail}; give the lattice of \
             the history (ADR-0187)"
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
    // 初めの重金属の質量は、燃焼（ADR-0173）と同じ式で、燃焼度 0 の数密度から。履歴はチェーンの原子番号を持たないので、
    // ライブラリの ZA から引く（燃焼度 0 の核種は入力の組成で、上で入力と一致を確かめた。輸送するのでライブラリにある）
    let heavy_metal = heavy_metal_mass(&history.regions, &first.densities, |i| {
        let name = &history.nuclides[i];
        let data = library
            .nuclides
            .iter()
            .find(|n| &n.name == name)
            .ok_or_else(|| anyhow!("{name}: no mass in the library"))?;
        // ZA = 1000 Z + A
        let z = (data.za / 1000.0).floor();
        Ok((z >= f64::from(HEAVY_METAL_Z)).then_some(data.awr))
    })?;
    // 輸送する核種は、燃焼のライブラリが燃焼の温度で持つチェーンの核種（ADR-0173 の案 N1）。入力のライブラリと燃焼させた物質の温度
    // （公称）で求め直した集合が履歴と同じこと。違えば、燃焼と違うライブラリか温度で、公称の分岐が履歴の点の輸送でなくなる
    let mut nominal: Vec<f64> = Vec::new();
    for name in &depleted {
        let t = input
            .materials
            .iter()
            .find(|m| &m.name == name)
            .map(|m| m.temperature)
            .ok_or_else(|| anyhow!("depleted material `{name}` is not in the input"))?;
        if !nominal.contains(&t) {
            nominal.push(t);
        }
    }
    let expected: Vec<&String> = history
        .nuclides
        .iter()
        .filter(|n| nominal.iter().all(|&t| has_data_at(library, n, t)))
        .collect();
    if expected.iter().copied().ne(history.transportable.iter()) {
        bail!(
            "the library has data at {nominal:?} K for {} chain nuclides, the history transported {}; give the library \
             and the depletion temperature of the history (ADR-0187)",
            expected.len(),
            history.transportable.len()
        );
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
    // 輸送しない核種（崩壊と生成だけを追った核種）。燃焼の輸送と同じく輸送に入れず、燃焼点ごとに示す（ADR-0173）
    let in_library = |i: usize| {
        library
            .nuclides
            .iter()
            .any(|n| n.name == history.nuclides[i])
    };
    let total_volume: f64 = history.regions.iter().map(|d| d.volume).sum();
    for &p in &points {
        let densities = &history.points[p].densities;
        let present = |i: usize| densities.iter().any(|n| n[i] > 0.0);
        let total = |i: usize| -> f64 {
            history
                .regions
                .iter()
                .zip(densities)
                .map(|(d, n)| d.volume * n[i].max(0.0))
                .sum::<f64>()
                / total_volume
        };
        let left: Vec<usize> = (0..history.nuclides.len())
            .filter(|i| !transported.contains(i) && present(*i))
            .collect();
        let (lacking, absent): (Vec<usize>, Vec<usize>) =
            left.into_iter().partition(|&i| in_library(i));
        if !lacking.is_empty() {
            tracing::warn!(
                burnup = history.points[p].burnup,
                nuclides = lacking
                    .iter()
                    .map(|&i| history.nuclides[i].as_str())
                    .collect::<Vec<_>>()
                    .join(" "),
                density = lacking.iter().map(|&i| total(i)).sum::<f64>(),
                "not transported: no library data at the temperature of the depletion"
            );
        }
        tracing::info!(
            burnup = history.points[p].burnup,
            nuclides = absent.len(),
            density = absent.iter().map(|&i| total(i)).sum::<f64>(),
            "not transported: not in the library"
        );
    }
    // 分岐の状態（と参照の状態）ごとの、燃焼させない物質と、燃焼の領域の温度。輸送する核種がその温度のデータを持つこと（近い
    // 温度で代用しない。ADR-0173 の案 N2）。温度は物質ごとなので、物質ごとに確かめる
    let split = |state: &[MaterialInput]| -> anyhow::Result<(Vec<MaterialInput>, Vec<f64>)> {
        let mut by_material: Vec<(&str, f64)> = Vec::with_capacity(depleted.len());
        for name in &depleted {
            let temperature = state
                .iter()
                .find(|m| &m.name == name)
                .ok_or_else(|| anyhow!("depleted material `{name}` is not in the state"))?
                .temperature;
            if let Some(&i) = transported
                .iter()
                .find(|&&i| !has_data_at(library, &history.nuclides[i], temperature))
            {
                bail!(
                    "{} (transported in the history) has no library data at {temperature} K of material `{name}`",
                    history.nuclides[i]
                );
            }
            by_material.push((name, temperature));
        }
        let temperatures = history
            .regions
            .iter()
            .map(|d| {
                by_material
                    .iter()
                    .find(|(name, _)| *name == d.material)
                    .map(|&(_, t)| t)
                    .ok_or_else(|| anyhow!("depleted material `{}` has no temperature", d.material))
            })
            .collect::<anyhow::Result<Vec<_>>>()?;
        let kept: Vec<MaterialInput> = state
            .iter()
            .filter(|m| !depleted.contains(&m.name))
            .cloned()
            .collect();
        Ok((kept, temperatures))
    };
    let branch = branches
        .states()
        .iter()
        .map(|state| split(state))
        .collect::<anyhow::Result<Vec<_>>>()?;
    let references = branches
        .references()
        .iter()
        .map(|state| split(state))
        .collect::<anyhow::Result<Vec<_>>>()?;
    let mut axes = vec![Axis::Burnup {
        points: selection.burnup.clone(),
    }];
    axes.extend(branches.axes().iter().cloned());
    Ok(HistoryBranches {
        axes,
        region_materials: renamed,
        heavy_metal,
        regions: history.regions.clone(),
        densities: points
            .iter()
            .map(|&p| history.points[p].densities.clone())
            .collect(),
        nuclides: history.nuclides.clone(),
        transported,
        branch,
        reference_of: (0..branches.states().len())
            .map(|b| branches.reference(b))
            .collect(),
        references,
    })
}

impl HistoryBranches {
    /// 状態の数（燃焼度の点の数 × 分岐の状態の数）
    pub fn len(&self) -> usize {
        self.densities.len() * self.branch.len()
    }

    /// 状態がないか（`history_branches` は作らない）
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// 状態 `index`（軸の直積の順。燃焼度がいちばん遅く変わる）の解く物質。状態ごとに作り、すべての状態の物質を同時に
    /// 持たない（本番の大きさで、核種の名前の複製が状態の数 × 領域 × 核種に増えるため）
    pub fn state(&self, index: usize) -> Vec<MaterialInput> {
        let (point, b) = (index / self.branch.len(), index % self.branch.len());
        self.materials(point, &self.branch[b])
    }

    /// 状態 `index` の自己遮蔽を参照の状態から作るなら、その燃焼度の点と参照の番号（ADR-0201）
    pub fn reference(&self, index: usize) -> Option<(usize, usize)> {
        let (point, b) = (index / self.branch.len(), index % self.branch.len());
        self.reference_of[b].map(|r| (point, r))
    }

    /// 燃焼度の点 `point` の、参照 `reference`（ホウ素が公称）の解く物質
    pub fn reference_state(&self, point: usize, reference: usize) -> Vec<MaterialInput> {
        self.materials(point, &self.references[reference])
    }

    /// 燃焼度の点 `point` の、燃焼させない物質と燃焼の領域の温度が `branch` の状態の解く物質
    fn materials(
        &self,
        point: usize,
        branch: &(Vec<MaterialInput>, Vec<f64>),
    ) -> Vec<MaterialInput> {
        let (kept, temperatures) = branch;
        let mut materials = kept.clone();
        for ((d, densities), &temperature) in self
            .regions
            .iter()
            .zip(&self.densities[point])
            .zip(temperatures)
        {
            materials.push(region_material(
                d,
                temperature,
                densities,
                &self.transported,
                |i| self.nuclides[i].clone(),
            ));
        }
        materials
    }
}

#[cfg(test)]
mod tests;
