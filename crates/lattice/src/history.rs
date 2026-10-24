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
//! 状態の自己遮蔽は `[branches]` の `shielding` で選ぶ（[`BranchShielder`]。ADR-0201）。

use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{Context, anyhow, bail};
use xlc_format::library::{Library, LookupError, Shielded};
use xlc_format::table::Axis;

use crate::branches::BranchShielding;
use crate::depletion::{
    DepletionRegion, HEAVY_METAL_Z, History, check_negative, check_region_names, depletion_regions,
    has_data_at, heavy_metal_mass, region_material, same_cross_sections, solved_names,
};
use crate::geometry::Lattice;
use crate::input::{LatticeInput, MaterialInput};
use crate::material::Material;
use crate::shielding::{
    FuelShielding, ShieldingOptions, ShieldingReuse, shield_fuel_reusing, shielded_nuclide,
};
use crate::symmetry::{Chosen, Detected};
use crate::tracking::Tracking;

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
    /// 状態の自己遮蔽の解き方（ADR-0201）
    shielding: BranchShielding,
    /// [燃料の温度の格子点] の参照の状態（ホウ素と減速材が公称）の、燃焼させない物質と、燃焼の領域ごとの温度 (K)
    references: Vec<(Vec<MaterialInput>, Vec<f64>)>,
    /// ホウ素と減速材の軸が変える物質
    varied: Vec<String>,
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
    // 公称の減速材の自己遮蔽は、燃焼させた物質を参照の自己遮蔽から取り、ホウ素と減速材の軸が変える物質だけを作り直すので、燃焼させた
    // 物質がその軸に挙がれば状態の温度にならない
    if branches.shielding() == BranchShielding::NominalModerator
        && let Some(name) = branches.varied().iter().find(|n| depleted.contains(n))
    {
        bail!("depleted material `{name}` is changed by the boron or moderator axis");
    }
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
        shielding: branches.shielding(),
        references,
        varied: branches.varied().to_vec(),
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

    /// 状態 `index` の自己遮蔽を参照の状態から作るなら、その燃焼度の点と燃料の温度の格子点の番号（ADR-0201）
    pub fn reference(&self, index: usize) -> Option<(usize, usize)> {
        match self.shielding {
            BranchShielding::EveryState => None,
            // 燃料の温度の軸は最後の軸（いちばん速く変わる）
            BranchShielding::NominalModerator => {
                Some((index / self.branch.len(), index % self.references.len()))
            }
        }
    }

    /// 燃焼度の点 `point` の、燃料の温度の格子点 `fuel` の参照の状態（ホウ素と減速材が公称）の解く物質
    pub fn reference_state(&self, point: usize, fuel: usize) -> Vec<MaterialInput> {
        self.materials(point, &self.references[fuel])
    }

    /// ホウ素と減速材の軸が変える物質
    pub fn varied(&self) -> &[String] {
        &self.varied
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

/// 履歴からの分岐の状態を順に自己遮蔽する（ADR-0201）。`EveryState` は状態ごとに解き、`NominalModerator` は燃焼度の点と燃料の
/// 温度の格子点ごとに参照の状態（ホウ素と減速材が公称）で 1 度解いて持ち、状態の物質をそこから作る。燃料の組成によらない前処理
/// （Dancoff 係数）は `ShieldingReuse` で持ち越す
pub struct BranchShielder<'a> {
    tracking: &'a Tracking,
    reuse: ShieldingReuse<'a>,
    /// 今の燃焼度の点の、燃料の温度の格子点ごとの参照の自己遮蔽と、それに渡した物質
    held: Vec<(usize, usize, FuelShielding, Vec<MaterialInput>)>,
    /// 自己遮蔽を解いた回数
    shieldings: usize,
}

impl<'a> BranchShielder<'a> {
    /// トラッキング `tracking` の格子の状態を自己遮蔽する（まだ何も持たない）
    pub fn new(tracking: &'a Tracking) -> Self {
        Self {
            tracking,
            reuse: ShieldingReuse::new(tracking),
            held: Vec::new(),
            shieldings: 0,
        }
    }

    /// 自己遮蔽を解いた回数
    pub fn shieldings(&self) -> usize {
        self.shieldings
    }

    /// 状態 `index` の物質の一覧と、[領域] のその番号。状態は燃焼度の点の順に渡す（参照は今の燃焼度の点のものだけを持つ）
    pub fn shield(
        &mut self,
        library: &Library,
        branches: &HistoryBranches,
        index: usize,
        symmetry: &Detected,
        options: &ShieldingOptions,
    ) -> anyhow::Result<(Vec<Material>, Vec<usize>)> {
        let materials = branches.state(index);
        let names = &branches.region_materials;
        let Some((point, fuel)) = branches.reference(index) else {
            self.shieldings += 1;
            let shielded = shield_fuel_reusing(
                library,
                &mut self.reuse,
                &materials,
                names,
                symmetry,
                options,
            )?;
            return Ok((shielded.materials, shielded.region_material));
        };
        let lattice = &self.tracking.lattice;
        self.held.retain(|(p, ..)| *p == point);
        if !self.held.iter().any(|(_, f, ..)| *f == fuel) {
            let nominal = branches.reference_state(point, fuel);
            self.shieldings += 1;
            let shielded =
                shield_fuel_reusing(library, &mut self.reuse, &nominal, names, symmetry, options)?;
            // 作り直しの前提（自己遮蔽の物質が、領域の実効断面積と組成から作った群定数に SPH 因子を掛けたもの）を、公称の組成で
            // 作り直した物質が自己遮蔽の物質とビットで同じことで確かめる。崩れれば、状態の物質が自己遮蔽の作り方と食い違う
            let rebuilt = remix_varied(
                library,
                &shielded,
                lattice,
                names,
                &branches.varied,
                &nominal,
                &nominal,
            )
            .context("remixing the nominal materials")?;
            if let Some((m, _)) = rebuilt
                .iter()
                .find(|(m, material)| *material != shielded.materials[*m])
            {
                bail!(
                    "material {m} remixed from its shielded cross sections differs from the self-shielding at the \
                     nominal boron and moderator"
                );
            }
            self.held.push((point, fuel, shielded, nominal));
        }
        let (_, _, reference, nominal) = self
            .held
            .iter()
            .find(|(_, f, ..)| *f == fuel)
            .ok_or_else(|| anyhow!("no reference shielding for fuel temperature {fuel}"))?;
        let mut unique = reference.materials.clone();
        if materials != *nominal {
            let rebuilt = remix_varied(
                library,
                reference,
                lattice,
                names,
                &branches.varied,
                nominal,
                &materials,
            )
            .context("remixing the materials of the boron and moderator axes")?;
            for (m, material) in rebuilt {
                unique[m] = material;
            }
        }
        Ok((unique, reference.region_material.clone()))
    }
}

/// 公称の状態の物質 `nominal` で解いた自己遮蔽 `reference` の実効断面積と SPH 因子のまま、物質の名前が `varied` の領域の物質を、状態の
/// 物質 `state` の組成と温度で作り直す（ADR-0201）。`region_materials` は [領域] の物質の名前で、`nominal` と `state` は同じ名前を
/// 同じ順に持つ。返すのは、作り直した `reference.materials` の番号と物質。自己遮蔽した物質は、自己遮蔽と同じく、核種の実効断面積から
/// 作った群定数を数密度で混ぜて SPH 因子を掛け、自己遮蔽しなかった物質は `Material::mix`（無限希釈）で作る。温度が公称と違えば、
/// 核種の実効断面積は状態の温度の無限希釈の値に公称の実効断面積と無限希釈の値の比を掛けたもの（[`carry_shielding`]）、公称の組成に
/// ない核種は状態の温度の無限希釈の値
pub(crate) fn remix_varied(
    library: &Library,
    reference: &FuelShielding,
    lattice: &Lattice,
    region_materials: &[String],
    varied: &[String],
    nominal: &[MaterialInput],
    state: &[MaterialInput],
) -> anyhow::Result<Vec<(usize, Material)>> {
    if nominal.len() != state.len() || nominal.iter().zip(state).any(|(a, b)| a.name != b.name) {
        bail!("the nominal and the state materials differ in their names");
    }
    if region_materials.len() != reference.region_material.len() {
        bail!(
            "{} region materials for {} regions",
            region_materials.len(),
            reference.region_material.len()
        );
    }
    let ng = library.groups.len();
    // 領域ごとの SPH 因子 μ（[領域][群]）。物理的な中性子束は MOC の中性子束に μ を掛けたものなので、1 の中性子束から引く
    let ones = vec![vec![1.0; ng]; lattice.regions()];
    let sph = reference.physical_flux(lattice, &ones)?;
    // 物質の番号ごとの、作り直す物質の名前と最初の領域。共有する領域の実効断面積と SPH 因子は最初の領域のものを使うので、同じこと。
    // 作り直さない物質の領域と共有していれば、作り直せない
    let mut owner: Vec<Option<(&str, usize)>> = vec![None; reference.materials.len()];
    let mut kept = vec![false; reference.materials.len()];
    for (r, name) in region_materials.iter().enumerate() {
        let m = reference.region_material[r];
        if !varied.contains(name) {
            kept[m] = true;
            continue;
        }
        match owner[m] {
            None => owner[m] = Some((name, r)),
            Some((other, _)) if other != name => {
                bail!("material {m} is shared by `{other}` and `{name}`")
            }
            Some((_, f)) => {
                if sph[r] != sph[f]
                    || !same_cross_sections(
                        reference.region_cross_sections(r),
                        reference.region_cross_sections(f),
                    )
                {
                    bail!(
                        "regions {f} and {r} share material {m} with different cross sections or SPH factors"
                    );
                }
            }
        }
    }
    if let Some(m) = (0..owner.len()).find(|&m| owner[m].is_some() && kept[m]) {
        bail!("material {m} is shared with a material the branches keep");
    }
    let infinite = *library.sigma0.first().ok_or(LookupError::NoSigma0)?;
    let dilute =
        |name: &str, temperature: f64| library.shielded(name, temperature, &vec![infinite; ng]);
    owner
        .iter()
        .enumerate()
        .filter_map(|(m, o)| o.map(|(name, r)| (m, name, r)))
        .map(|(m, name, r)| -> anyhow::Result<(usize, Material)> {
            let k = nominal
                .iter()
                .position(|x| x.name == name)
                .ok_or_else(|| anyhow!("material `{name}` is not in the state"))?;
            let (before, after) = (&nominal[k], &state[k]);
            // 自己遮蔽しなかった物質の領域は、自己遮蔽に渡した物質の一覧と同じ番号の、`Material::mix` の物質（自己遮蔽した物質は
            // 一覧と真空の後に足される）。自己遮蔽した群定数の作り方（`shielded_nuclide`）で作ると、散乱行列の補正の丸めで食い違う
            if m == k {
                return Ok((
                    m,
                    Material::mix(library, &after.composition, after.temperature)?,
                ));
            }
            let cross_sections = reference
                .region_cross_sections(r)
                .ok_or_else(|| anyhow!("region {r} has no shielded cross sections"))?;
            if cross_sections.len() != before.composition.len() {
                bail!(
                    "region {r}: {} shielded nuclides for the {} nuclides of `{name}`",
                    cross_sections.len(),
                    before.composition.len()
                );
            }
            let nuclides = after
                .composition
                .iter()
                .map(|(nuclide, density)| {
                    let at = before.composition.iter().position(|(n, _)| n == nuclide);
                    let shielded = match at {
                        Some(i) if after.temperature == before.temperature => {
                            cross_sections[i].clone()
                        }
                        Some(i) => carry_shielding(
                            &dilute(nuclide, before.temperature)?,
                            &dilute(nuclide, after.temperature)?,
                            &cross_sections[i],
                        )
                        .with_context(|| format!("{nuclide} of `{name}`"))?,
                        None => dilute(nuclide, after.temperature)?,
                    };
                    Ok((
                        shielded_nuclide(library, nuclide, after.temperature, &shielded)?,
                        *density,
                    ))
                })
                .collect::<anyhow::Result<Vec<_>>>()?;
            let mixed = Material::from_nuclides(&nuclides, ng)?.scaled(&sph[r])?;
            Ok((m, mixed))
        })
        .collect()
}

/// 公称の温度の無限希釈の値 `nominal` に対する実効断面積 `shielded` の比（自己遮蔽因子）を、状態の温度の無限希釈の値 `state` に
/// 掛けた実効断面積（ADR-0201）。遮蔽していない群（比が 1）は `state` の値そのもの。無限希釈の値が 0 の群は、実効断面積も 0 で
/// なければ比が決まらないのでエラー
pub(crate) fn carry_shielding(
    nominal: &Shielded,
    state: &Shielded,
    shielded: &Shielded,
) -> anyhow::Result<Shielded> {
    let carry = |reaction: &str, n: &[f64], s: &[f64], e: &[f64]| -> anyhow::Result<Vec<f64>> {
        if n.len() != s.len() || n.len() != e.len() {
            bail!(
                "{reaction}: {} groups, {} at the state, {} shielded",
                n.len(),
                s.len(),
                e.len()
            );
        }
        n.iter()
            .zip(s)
            .zip(e)
            .enumerate()
            .map(|(g, ((&n, &s), &e))| {
                if n != 0.0 {
                    Ok(s * (e / n))
                } else if e == 0.0 {
                    Ok(s)
                } else {
                    bail!("{reaction}: group {g} is {e} shielded but 0 at infinite dilution")
                }
            })
            .collect()
    };
    Ok(Shielded {
        total: carry("total", &nominal.total, &state.total, &shielded.total)?,
        elastic: carry(
            "elastic",
            &nominal.elastic,
            &state.elastic,
            &shielded.elastic,
        )?,
        fission: carry(
            "fission",
            &nominal.fission,
            &state.fission,
            &shielded.fission,
        )?,
        nu_fission: carry(
            "nu-fission",
            &nominal.nu_fission,
            &state.nu_fission,
            &shielded.nu_fission,
        )?,
        capture: carry(
            "capture",
            &nominal.capture,
            &state.capture,
            &shielded.capture,
        )?,
    })
}

#[cfg(test)]
mod tests;
