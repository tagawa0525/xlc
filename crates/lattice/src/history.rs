//! 燃焼の履歴の燃焼点からの分岐（ADR-0187）。
//!
//! 入力の `[history]` が挙げた燃焼点ごとに、燃焼の履歴（ADR-0173）の燃焼の領域の数密度を固定し、`[branches]` の状態
//! （ADR-0097）で輸送する。分岐の状態では燃焼させない。燃焼の領域は燃焼の輸送と同じ `<物質>@<セル>` の物質として解き、
//! 輸送する核種は履歴の `transportable`（数密度が 0 の核種も含め、負の値は 0 にする）。状態は燃焼度がいちばん遅く変わり、
//! 各燃焼度の点の中は分岐の軸の直積の順。解くのは、公称の状態（`[branches]` の公称）をすべての燃焼点で、ほかの部分格子の状態を
//! `off-nominal-burnup` の燃焼点で（ADR-0200）。
//!
//! ```toml
//! [history]
//! file = "history.msgpack"   # 燃焼の履歴（入力のファイルのあるディレクトリから）
//! burnup = [0, 10, 20, 40]   # 分岐を解く燃焼点 (MWd/kg)。履歴の燃焼度から選ぶ。公称の状態はすべての点で解く
//! off-nominal-burnup = [0, 20, 40]   # 公称でない状態を解く燃焼点（burnup から選ぶ。省略すると burnup のすべて。ADR-0200）
//! ```
//!
//! 状態の自己遮蔽は `[branches]` の `shielding` で選ぶ（[`crate::branches::BranchShielder`]。ADR-0201）。

use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{anyhow, bail};
use xlc_format::library::{Library, Shielded};
use xlc_format::table::{Axis, MicroNuclide, MicroReaction};

use crate::depletion::{
    DepletionRegion, HEAVY_METAL_Z, History, channel_cross_sections, check_negative,
    check_region_names, depletion_regions, has_data_at, heavy_metal_mass, region_material,
    solved_names, split_at, yield_cutoff,
};
use crate::fewgroup::coarse_groups;
use crate::geometry::Lattice;
use crate::input::{LatticeInput, MaterialInput};
use crate::moc::Solution;
use crate::shielding::{FuelShielding, shielded_nuclide_with};
use crate::symmetry::Chosen;
use crate::ufg::NuclideTransfer;

/// 微視的断面積を足す MOC の領域の組: 核種ごとの実効断面積、UFG の範囲の移行、重みを付けた細かい群の中性子束の和
type Member<'a> = (&'a [Shielded], Option<&'a [NuclideTransfer]>, Vec<f64>);

/// 入力の `[history]`（ADR-0187）。
#[derive(Debug, Clone, PartialEq)]
pub struct HistoryInput {
    /// 燃焼の履歴のファイル（入力のファイルのあるディレクトリから）
    pub file: PathBuf,
    /// 分岐を解く燃焼点 (MWd/kg)。0 以上で増え、どれも履歴の燃焼度のどれかと等しい
    pub burnup: Vec<f64>,
    /// 公称でない状態を解く燃焼点 (MWd/kg)。`burnup` から選び、増える（ADR-0200）
    pub off_nominal: Vec<f64>,
    /// ミクロ燃焼で追う核種（ADR-0216）。空なら微視的断面積を作らない
    pub micro: Vec<String>,
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
    /// 分岐の状態ごとに、解くか（`Branches::solved`）
    solved: Vec<bool>,
    /// 公称の分岐の状態（`Branches::nominal_state`）
    nominal: Option<usize>,
    /// 燃焼度の点ごとに、公称でない状態を解くか
    off_nominal: Vec<bool>,
    /// ミクロ燃焼で追う核種の名前（ADR-0216）
    micro_nuclides: Vec<String>,
    /// 追う核種の、輸送する核種（`transported`。燃焼の領域の物質の核種の並び）の中の番号
    micro_positions: Vec<usize>,
    /// [燃焼度の点][追う核種] の均質化した数密度 (1/(barn cm))
    micro_densities: Vec<Vec<f64>>,
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
    // 公称でない状態の燃焼点を間引くと、その差を公称の状態から内挿するので、公称の状態が要る（ADR-0200）
    let off_nominal: Vec<bool> = selection
        .burnup
        .iter()
        .map(|b| selection.off_nominal.contains(b))
        .collect();
    if off_nominal.contains(&false) && branches.nominal_state().is_none() {
        bail!(
            "off-nominal-burnup: the off-nominal states are interpolated from the nominal state, the input materials, \
             which are not at a point of every branch axis (ADR-0200)"
        );
    }
    // ミクロ燃焼で追う核種は、履歴が輸送した核種（ADR-0216。輸送しない核種は自己遮蔽した断面積を持たない）
    let micro_positions = selection
        .micro
        .iter()
        .map(|name| {
            let i = index.get(name.as_str()).ok_or_else(|| {
                anyhow!("micro nuclide {name} is not a nuclide of the history (ADR-0216)")
            })?;
            transported.iter().position(|t| t == i).ok_or_else(|| {
                anyhow!(
                    "micro nuclide {name} is not transported in the history (no library data at the depletion \
                     temperature. ADR-0216)"
                )
            })
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    // 均質化した数密度は、燃焼の領域の数密度（負は輸送と同じく 0）を体積で重みを付けて格子の体積で割ったもの
    let lattice_volume: f64 = lattice.volumes().iter().sum();
    let micro_densities = points
        .iter()
        .map(|&p| {
            micro_positions
                .iter()
                .map(|&k| {
                    let i = transported[k];
                    history
                        .regions
                        .iter()
                        .zip(&history.points[p].densities)
                        .map(|(d, n)| d.volume * n[i].max(0.0))
                        .sum::<f64>()
                        / lattice_volume
                })
                .collect()
        })
        .collect();
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
        solved: (0..branches.states().len())
            .map(|b| branches.solved(b))
            .collect(),
        nominal: branches.nominal_state(),
        off_nominal,
        micro_nuclides: selection.micro.clone(),
        micro_positions,
        micro_densities,
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

    /// 状態 `index`（軸の直積の順）を解くか: 分岐の部分格子の状態で、公称の状態か、公称でない状態を解く燃焼点のもの（ADR-0200）
    pub fn solved(&self, index: usize) -> bool {
        let (point, b) = (index / self.branch.len(), index % self.branch.len());
        self.solved.get(b).copied().unwrap_or(false)
            && (self.nominal == Some(b) || self.off_nominal.get(point).copied().unwrap_or(false))
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

    /// ミクロ燃焼で追う核種の名前（入力の `micro` の順）
    pub fn micro_nuclides(&self) -> &[String] {
        &self.micro_nuclides
    }

    /// [燃焼度の点][追う核種] の、格子で均質化した数密度 (1/(barn cm))
    pub fn micro_densities(&self) -> &[Vec<f64>] {
        &self.micro_densities
    }

    /// 状態 `index` の追う核種の、少数群（境界 `bounds`）の微視的断面積 (barn)（ADR-0216）。`shielding` は状態の燃料の自己遮蔽、
    /// `solution` はその MOC の解。少数群定数（`fewgroup::homogenize`）と同じく MOC の中性子束と体積で重みを付け、物質と同じく SPH
    /// 因子を掛けた核種の群定数を、燃焼の領域の数密度で重みを付けて足し、均質化した数密度と少数群の中性子束で割る。均質化した
    /// 数密度が 0 の核種は、燃焼の領域で 1 として重みを付ける
    pub fn micro(
        &self,
        index: usize,
        library: &Library,
        shielding: &FuelShielding,
        lattice: &Lattice,
        solution: &Solution,
        bounds: &[f64],
    ) -> anyhow::Result<Vec<MicroNuclide>> {
        let (point, b) = (index / self.branch.len(), index % self.branch.len());
        let temperatures = &self.branch[b].1;
        let fine = library.groups.len();
        let group_of = coarse_groups(&library.groups, bounds)?;
        let coarse = bounds.len() - 1;
        // 領域ごとの SPH 因子 μ（物理的な中性子束は MOC の中性子束に μ を掛けたもの）
        let ones = vec![vec![1.0; fine]; lattice.regions()];
        let sph = shielding.physical_flux(lattice, &ones)?;
        if solution.flux.len() != lattice.regions() || solution.volumes.len() != lattice.regions() {
            bail!(
                "the solution has {} regions for {} of the lattice",
                solution.flux.len(),
                lattice.regions()
            );
        }
        let mut flux = vec![0.0; coarse];
        for (phi, volume) in solution.flux.iter().zip(&solution.volumes) {
            for (g, f) in phi.iter().enumerate() {
                flux[group_of[g]] += volume * f;
            }
        }
        if flux.iter().any(|&f| !(f > 0.0 && f.is_finite())) {
            bail!("the few-group fluxes {flux:?} must be positive");
        }
        let lattice_volume: f64 = lattice.volumes().iter().sum();
        let depleted_volume: f64 = self.regions.iter().map(|d| d.volume).sum();
        let mut out = Vec::with_capacity(self.micro_positions.len());
        for (j, &k) in self.micro_positions.iter().enumerate() {
            let name = &self.nuclides[self.transported[k]];
            let mean = self.micro_densities[point][j];
            let mut total = vec![0.0; coarse];
            let mut nu_fission = vec![0.0; coarse];
            let mut fission = vec![0.0; coarse];
            let mut scatter = vec![vec![0.0; coarse]; coarse];
            // 細かい群の核分裂率（収率の重み。ADR-0183）
            let mut fine_fission = vec![0.0; fine];
            let mut reactions: Option<Vec<MicroReaction>> = None;
            for ((d, n), &temperature) in self
                .regions
                .iter()
                .zip(&self.densities[point])
                .zip(temperatures)
            {
                let weight = if mean > 0.0 {
                    n[self.transported[k]].max(0.0)
                } else {
                    1.0
                };
                if weight == 0.0 {
                    continue;
                }
                let (_, data) = library.lookup(name, temperature)?;
                // 実効断面積と移行が同じ MOC の領域（同じ配列を指す。自己遮蔽の結果）は、重みを付けた細かい群の中性子束を
                // 足してから 1 度だけ核種の群定数にする
                let mut members: Vec<Member> = Vec::new();
                for &r in &d.regions {
                    let effective = shielding.region_cross_sections(r).ok_or_else(|| {
                        anyhow!(
                            "region {r} of material `{}` has no shielded cross sections",
                            d.material
                        )
                    })?;
                    if effective.len() != self.transported.len() {
                        bail!(
                            "region {r} of material `{}`: {} shielded nuclides for {} transported",
                            d.material,
                            effective.len(),
                            self.transported.len()
                        );
                    }
                    let transfers = shielding.region_transfers(r);
                    let weighted: Vec<f64> = (0..fine)
                        .map(|g| weight * solution.volumes[r] * solution.flux[r][g] * sph[r][g])
                        .collect();
                    let same = |(e, t, _): &&mut Member| {
                        std::ptr::eq(*e, effective)
                            && match (t, transfers) {
                                (Some(a), Some(b)) => std::ptr::eq(*a, b),
                                (a, b) => a.is_none() && b.is_none(),
                            }
                    };
                    match members.iter_mut().find(|m| same(m)) {
                        Some((_, _, sum)) => {
                            for (s, w) in sum.iter_mut().zip(&weighted) {
                                *s += w;
                            }
                        }
                        None => members.push((effective, transfers, weighted)),
                    }
                }
                for (effective, transfers, weighted) in members {
                    let shielded = &effective[k];
                    let nuclide = shielded_nuclide_with(
                        library,
                        name,
                        temperature,
                        shielded,
                        transfers.map(|t| &t[k]),
                    )?;
                    for (g, &w) in weighted.iter().enumerate() {
                        let big = group_of[g];
                        total[big] += w * nuclide.total[g];
                        nu_fission[big] += w * nuclide.nu_fission[g];
                        fission[big] += w * nuclide.fission[g];
                        for (to, x) in nuclide.scatter[g].iter().enumerate() {
                            scatter[big][group_of[to]] += w * x;
                        }
                        fine_fission[g] += w * shielded.fission[g];
                    }
                    let channels = data
                        .reactions
                        .iter()
                        .map(|reaction| {
                            channel_cross_sections(name, reaction, shielded)
                                .map(|c| (reaction.mt, c))
                                .map_err(|e| anyhow!(e))
                        })
                        .collect::<anyhow::Result<Vec<_>>>()?;
                    let sums = reactions.get_or_insert_with(|| {
                        channels
                            .iter()
                            .map(|(mt, c)| MicroReaction {
                                mt: *mt,
                                channels: vec![vec![0.0; coarse]; c.len()],
                            })
                            .collect()
                    });
                    let shape = |r: &[MicroReaction]| -> Vec<(u32, usize)> {
                        r.iter().map(|r| (r.mt, r.channels.len())).collect()
                    };
                    if shape(sums)
                        != channels
                            .iter()
                            .map(|(mt, c)| (*mt, c.len()))
                            .collect::<Vec<_>>()
                    {
                        bail!(
                            "{name}: the chain reactions of the library differ between temperatures"
                        );
                    }
                    for (sum, (_, c)) in sums.iter_mut().zip(&channels) {
                        for (out, sigma) in sum.channels.iter_mut().zip(c) {
                            for (g, &w) in weighted.iter().enumerate() {
                                out[group_of[g]] += w * sigma[g];
                            }
                        }
                    }
                }
            }
            // 均質化した数密度 × 少数群の中性子束で割る。数密度が 0 なら、燃焼の領域の体積の割合（重み 1 の数密度）
            let density = if mean > 0.0 {
                mean
            } else {
                depleted_volume / lattice_volume
            };
            let per = |values: Vec<f64>| -> Vec<f64> {
                values
                    .iter()
                    .zip(&flux)
                    .map(|(v, f)| v / (density * f))
                    .collect()
            };
            let scatter = scatter
                .into_iter()
                .zip(&flux)
                .map(|(row, f)| row.iter().map(|v| v / (density * f)).collect())
                .collect();
            let reactions = reactions
                .unwrap_or_default()
                .into_iter()
                .map(|r| MicroReaction {
                    mt: r.mt,
                    channels: r.channels.into_iter().map(&per).collect(),
                })
                .collect();
            let (below, above) = split_at(library.groups.bounds(), &fine_fission, yield_cutoff());
            out.push(MicroNuclide {
                total: per(total),
                nu_fission: per(nu_fission),
                fission: per(fission),
                scatter,
                reactions,
                thermal_fission: (below + above > 0.0).then(|| below / (below + above)),
            });
        }
        Ok(out)
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
