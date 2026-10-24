//! 格子計算の分岐（ADR-0097）。入力の `[branches]` の軸（ホウ素濃度、減速材の温度と密度、燃料の温度）から、
//! 軸の直積の順（最初の軸がいちばん遅く変わる）の物質の組を作る。
//!
//! ```toml
//! [branches]
//! boron = { ppm = [900, 1300, 1700], materials = ["water"] }   # ホウ酸水のホウ素濃度
//! # 温度と、ホウ酸水の密度（格子点ごと）。挙げた物質の温度を変え、ホウ酸水なら密度も変える
//! moderator = { temperature = [550, 565, 580], density = [0.76972, 0.74271, 0.71190], materials = ["water", "clad"] }
//! fuel = { temperature = [550, 565, 580], materials = ["fuel"] }
//! shielding = "nominal-boron" # 自己遮蔽（"nominal-boron" か "every-state"。省略すると nominal-boron。ADR-0201）
//! ```
//!
//! 軸は boron、moderator、fuel の順に並べ、3つとも書く（core は状態を与えるとき3つの軸をすべて求める）。温度の軸が変える物質は名前で挙げ、物質の種類から
//! 推測しない。同じ物質を2つの温度の軸に挙げられない。
//!
//! 体積の割合で混ぜた物質（ADR-0128）は、軸が変える成分を持てばその軸に挙げ、持たなければ挙げない（混ぜた物質は
//! 状態ごとの成分から作る）。細かい群の媒質は燃料の温度の軸にだけ挙げ、分岐の状態ごとに媒質の同じ状態を使う。

use std::collections::BTreeMap;
use std::rc::Rc;

use anyhow::{Context, anyhow, bail};
use serde::Deserialize;
use xlc_format::fewgroup::FewGroup;
use xlc_format::library::{Library, LookupError};
use xlc_format::table::{Axis, Table};

use crate::depletion::same_cross_sections;
use crate::geometry::Lattice;
use crate::input::{MaterialInput, MediumInput, Mixture};
use crate::material::Material;
use crate::shielding::{
    FuelShielding, ShieldingOptions, ShieldingReuse, shield_fuel_reusing, shielded_nuclide,
};
use crate::symmetry::Detected;
use crate::tracking::Tracking;
use crate::water::borated_water;

/// ホウ酸水で書いた物質の量（ADR-0106）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Water {
    pub density: f64,
    pub boron: f64,
    pub b10: f64,
}

/// 分岐の状態の自己遮蔽の解き方（ADR-0201）
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BranchShielding {
    /// 状態ごとに、その状態の物質で解き直す
    EveryState,
    /// 減速材と燃料の温度の格子点ごとに（履歴からの分岐は燃焼度の点ごとにも）、ホウ素を公称（入力の物質）にした参照の状態で
    /// 1 度だけ解く。状態の物質は、参照の自己遮蔽の物質のうち、ホウ素の軸が変える物質だけを状態の組成で作り直す（[`BranchShielder`]）
    NominalBoron,
}

/// 分岐の軸と、軸の直積の順の物質の組。
#[derive(Debug, Clone, PartialEq)]
pub struct Branches {
    axes: Vec<Axis>,
    states: Vec<Vec<MaterialInput>>,
    shielding: BranchShielding,
    /// [減速材と燃料の温度の格子点の直積] の参照の物質の組（ホウ素の軸だけを入力のままにしたもの）
    references: Vec<Vec<MaterialInput>>,
    /// ホウ素の軸が変える物質
    varied: Vec<String>,
}

impl Branches {
    /// テーブルの軸（boron、moderator、fuel の順）
    pub fn axes(&self) -> &[Axis] {
        &self.axes
    }

    /// 軸の直積の順の物質の組
    pub fn states(&self) -> &[Vec<MaterialInput>] {
        &self.states
    }

    /// 状態の自己遮蔽の解き方
    pub fn shielding(&self) -> BranchShielding {
        self.shielding
    }

    /// [減速材と燃料の温度の格子点の直積（減速材が遅く変わる）] の参照の物質の組: 減速材と燃料の温度の軸はその格子点、ホウ素の
    /// 軸が変える物質は入力のまま
    pub fn references(&self) -> &[Vec<MaterialInput>] {
        &self.references
    }

    /// 状態 `state` の自己遮蔽を参照から作るなら、その参照の `references` の番号（`every-state` は `None`）
    pub fn reference(&self, state: usize) -> Option<usize> {
        match self.shielding {
            BranchShielding::EveryState => None,
            // ホウ素は最初の軸（いちばん遅く変わる）なので、ほかの軸の添字は状態の番号をそれらの点の数の積で割った余り
            BranchShielding::NominalBoron => Some(state % self.references.len()),
        }
    }

    /// ホウ素の軸が変える物質
    pub fn varied(&self) -> &[String] {
        &self.varied
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct BranchesEntry {
    boron: Option<BoronEntry>,
    moderator: Option<ModeratorEntry>,
    fuel: Option<FuelEntry>,
    shielding: Option<ShieldingEntry>,
}

/// 分岐の状態の自己遮蔽の解き方（ADR-0201）
#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
enum ShieldingEntry {
    EveryState,
    NominalBoron,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BoronEntry {
    ppm: Vec<f64>,
    materials: Vec<String>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ModeratorEntry {
    temperature: Vec<f64>,
    density: Vec<f64>,
    materials: Vec<String>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct FuelEntry {
    temperature: Vec<f64>,
    materials: Vec<String>,
}

/// 軸の点の誤り（空、有限でない、範囲の外、増えない）。
fn check_points(axis: &str, points: &[f64], valid: impl Fn(f64) -> bool) -> Result<(), String> {
    if points.is_empty()
        || !points.iter().all(|&p| p.is_finite() && valid(p))
        || !points.windows(2).all(|w| w[0] < w[1])
    {
        return Err(format!(
            "branches.{axis}: points must be valid and increasing, got {points:?}"
        ));
    }
    Ok(())
}

/// 挙げた物質が空でなく、重ならず、定義されていること（核種かホウ酸水で書いた物質、混ぜた物質、媒質）。
/// 混ぜた物質は、軸が変える成分を持つものだけを、すべて挙げること（ADR-0128）
fn check_materials(
    axis: &str,
    listed: &[String],
    materials: &[MaterialInput],
    mixtures: &[Mixture],
    media: &[MediumInput],
) -> Result<(), String> {
    if listed.is_empty() {
        return Err(format!("branches.{axis}: list the materials it changes"));
    }
    for (i, name) in listed.iter().enumerate() {
        if listed[..i].contains(name) {
            return Err(format!(
                "branches.{axis}: material `{name}` is listed twice"
            ));
        }
        if !materials.iter().any(|m| &m.name == name)
            && !mixtures.iter().any(|m| &m.name == name)
            && !media.iter().any(|m| &m.name == name)
        {
            return Err(format!("branches.{axis}: unknown material `{name}`"));
        }
    }
    for mixture in mixtures {
        let changed = mixture.parts.iter().any(|(part, _)| listed.contains(part));
        match (changed, listed.contains(&mixture.name)) {
            (true, false) => {
                return Err(format!(
                    "branches.{axis}: the mixture `{}` has a component it changes; list it too",
                    mixture.name
                ));
            }
            (false, true) => {
                return Err(format!(
                    "branches.{axis}: the mixture `{}` has no component it changes",
                    mixture.name
                ));
            }
            _ => {}
        }
    }
    Ok(())
}

/// 媒質は燃料の温度の軸にだけ挙げる（媒質の状態を選ぶ。ADR-0128）
fn check_media(axis: &str, listed: &[String], media: &[MediumInput]) -> Result<(), String> {
    for medium in media {
        match (axis == "fuel", listed.contains(&medium.name)) {
            (true, false) => {
                return Err(format!(
                    "branches.fuel: list the medium `{}` (the branches choose its states)",
                    medium.name
                ));
            }
            (false, true) => {
                return Err(format!(
                    "branches.{axis}: the medium `{}` belongs to the fuel axis only",
                    medium.name
                ));
            }
            _ => {}
        }
    }
    Ok(())
}

/// ホウ酸水の物質がすべて軸に挙がっていること（テーブルの状態は炉心の状態なので、一部だけを変えると
/// テーブルと計算が食い違う）。
fn check_all_waters(
    axis: &str,
    listed: &[String],
    waters: &BTreeMap<String, Water>,
) -> Result<(), String> {
    match waters.keys().find(|name| !listed.contains(name)) {
        Some(name) => Err(format!(
            "branches.{axis}: the water material `{name}` must be listed too"
        )),
        None => Ok(()),
    }
}

/// 状態ごとの出力の部分ノードの少数群定数（[状態][出力][部分ノード]。分けていなければ部分ノードは1つ）を、出力
/// ごとのテーブル（状態は軸の直積の順）にまとめる。
pub fn tables(
    axes: &[Axis],
    split: usize,
    per_state: Vec<Vec<Vec<FewGroup>>>,
) -> Result<Vec<Table>, String> {
    let expected: usize = axes.iter().map(|a| a.points().len()).product();
    if per_state.len() != expected {
        return Err(format!(
            "{} states for axes of {expected} points in all",
            per_state.len()
        ));
    }
    let outputs = per_state.first().map_or(0, Vec::len);
    let mut tables: Vec<Table> = (0..outputs)
        .map(|_| Table {
            axes: axes.to_vec(),
            split,
            states: Vec::with_capacity(expected * split * split),
            heavy_metal: None,
        })
        .collect();
    for (state, fews) in per_state.into_iter().enumerate() {
        if fews.len() != outputs {
            return Err(format!(
                "state {state} has {} outputs, state 0 has {outputs}",
                fews.len()
            ));
        }
        for (table, subs) in tables.iter_mut().zip(fews) {
            if subs.len() != split * split {
                return Err(format!(
                    "state {state} has {} sub-nodes, split {split} needs {}",
                    subs.len(),
                    split * split
                ));
            }
            table.states.extend(subs);
        }
    }
    Ok(tables)
}

/// 軸ごとの変更
enum Change {
    Boron(Vec<f64>, Vec<String>),
    Moderator(Vec<f64>, Vec<f64>, Vec<String>),
    Fuel(Vec<f64>, Vec<String>),
}

impl BranchesEntry {
    /// 軸を確かめ、軸の直積の順の物質の組（核種かホウ酸水で書いた物質）を作る。`waters` はホウ酸水で書いた物質の量、
    /// `mixtures` と `media` は軸に挙げる規則を確かめるだけ（ADR-0128）。
    pub(crate) fn build(
        self,
        materials: &[MaterialInput],
        waters: &BTreeMap<String, Water>,
        mixtures: &[Mixture],
        media: &[MediumInput],
    ) -> Result<Branches, String> {
        let mut axes = Vec::new();
        let mut changes = Vec::new();
        if let Some(b) = self.boron {
            check_points("boron", &b.ppm, |p| (0.0..1.0e6).contains(&p))?;
            check_materials("boron", &b.materials, materials, mixtures, media)?;
            check_media("boron", &b.materials, media)?;
            let mut b10 = None;
            for name in b
                .materials
                .iter()
                .filter(|n| !mixtures.iter().any(|m| &m.name == *n))
            {
                let water = waters.get(name).ok_or_else(|| {
                    format!("branches.boron: material `{name}` is not written as water")
                })?;
                if b10.is_some_and(|f| f != water.b10) {
                    return Err("branches.boron: the materials have different b10".into());
                }
                b10 = Some(water.b10);
            }
            let b10 = b10.ok_or("branches.boron: list the materials it changes")?;
            check_all_waters("boron", &b.materials, waters)?;
            axes.push(Axis::Boron {
                ppm: b.ppm.clone(),
                b10,
            });
            changes.push(Change::Boron(b.ppm, b.materials));
        }
        if let Some(m) = self.moderator {
            check_points("moderator", &m.temperature, |t| t > 0.0)?;
            if m.density.len() != m.temperature.len()
                || !m.density.iter().all(|&d| d > 0.0 && d.is_finite())
            {
                return Err(format!(
                    "branches.moderator: {} positive densities are needed, got {:?}",
                    m.temperature.len(),
                    m.density
                ));
            }
            check_materials("moderator", &m.materials, materials, mixtures, media)?;
            check_media("moderator", &m.materials, media)?;
            check_all_waters("moderator", &m.materials, waters)?;
            if waters.is_empty() {
                return Err("branches.moderator: no material is written as water".into());
            }
            axes.push(Axis::Moderator {
                temperature: m.temperature.clone(),
                density: m.density.clone(),
            });
            changes.push(Change::Moderator(m.temperature, m.density, m.materials));
        }
        if let Some(f) = self.fuel {
            check_points("fuel", &f.temperature, |t| t > 0.0)?;
            check_materials("fuel", &f.materials, materials, mixtures, media)?;
            check_media("fuel", &f.materials, media)?;
            if let Some(Change::Moderator(_, _, listed)) =
                changes.iter().find(|c| matches!(c, Change::Moderator(..)))
                // 混ぜた物質は成分から作るので、成分の温度が両方の軸で変わってよい
                && let Some(name) = f
                    .materials
                    .iter()
                    .find(|n| listed.contains(n) && !mixtures.iter().any(|m| &m.name == *n))
            {
                return Err(format!(
                    "branches: material `{name}` is in both the moderator and the fuel axes"
                ));
            }
            axes.push(Axis::Fuel {
                temperature: f.temperature.clone(),
            });
            changes.push(Change::Fuel(f.temperature, f.materials));
        }
        if axes.is_empty() {
            return Err("branches: no axes".into());
        }
        // core は状態を与えるとき3つの軸をすべて求める（ADR-0097）
        if axes.len() != 3 {
            return Err("branches: give all three axes (boron, moderator and fuel)".into());
        }

        // 軸ごとの格子点の添字（`None` の軸は入力のまま）の物質の組
        let state_of = |at: &[Option<usize>]| -> Result<Vec<MaterialInput>, String> {
            let mut state = materials.to_vec();
            let mut state_waters = waters.clone();
            for (change, i) in changes.iter().zip(at) {
                let Some(i) = *i else { continue };
                for m in &mut state {
                    match change {
                        Change::Boron(ppm, listed) if listed.contains(&m.name) => {
                            if let Some(w) = state_waters.get_mut(&m.name) {
                                w.boron = ppm[i];
                            }
                        }
                        Change::Moderator(temperature, density, listed)
                            if listed.contains(&m.name) =>
                        {
                            m.temperature = temperature[i];
                            if let Some(w) = state_waters.get_mut(&m.name) {
                                w.density = density[i];
                            }
                        }
                        Change::Fuel(temperature, listed) if listed.contains(&m.name) => {
                            m.temperature = temperature[i];
                        }
                        _ => {}
                    }
                }
            }
            for m in &mut state {
                if let Some(w) = state_waters.get(&m.name) {
                    m.composition = borated_water(w.density, w.boron, w.b10)
                        .map_err(|e| format!("material `{}`: {e}", m.name))?
                        .into_iter()
                        .collect();
                }
            }
            Ok(state)
        };
        let counts: Vec<usize> = axes.iter().map(|a| a.points().len()).collect();
        let total: usize = counts.iter().product();
        let states = (0..total)
            .map(|index| {
                let at: Vec<Option<usize>> = state_indices(&counts, index)
                    .into_iter()
                    .map(Some)
                    .collect();
                state_of(&at)
            })
            .collect::<Result<Vec<_>, _>>()?;
        // 参照はホウ素の軸（最初の軸）を入力のままにし、ほかの軸の格子点の直積を同じ順に並べる
        let rest: usize = counts[1..].iter().product();
        let references = (0..rest)
            .map(|index| {
                let at: Vec<Option<usize>> = std::iter::once(None)
                    .chain(state_indices(&counts[1..], index).into_iter().map(Some))
                    .collect();
                state_of(&at)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let varied = changes
            .iter()
            .find_map(|change| match change {
                Change::Boron(_, listed) => Some(listed.clone()),
                _ => None,
            })
            .ok_or("branches: no boron axis")?;
        Ok(Branches {
            axes,
            states,
            shielding: match self.shielding {
                None | Some(ShieldingEntry::NominalBoron) => BranchShielding::NominalBoron,
                Some(ShieldingEntry::EveryState) => BranchShielding::EveryState,
            },
            references,
            varied,
        })
    }
}

/// 分岐の状態を順に自己遮蔽する（ADR-0201）。参照の指定のない状態は、その状態の物質で解く。参照のある状態は、参照の状態（ホウ素が
/// 公称）で 1 度解いた自己遮蔽を持ち、ホウ素の軸が変える物質だけを状態の組成で作り直す。参照は今の燃焼度の点のものだけを持つ
/// （状態は燃焼度の点の順に渡す）。燃料の組成によらない前処理（Dancoff 係数）は `ShieldingReuse` で持ち越す
pub struct BranchShielder<'a> {
    tracking: &'a Tracking,
    reuse: ShieldingReuse<'a>,
    /// [領域] の自己遮蔽と MOC に渡す物質の名前
    region_materials: &'a [String],
    /// ホウ素の軸が変える物質
    varied: &'a [String],
    symmetry: &'a Detected,
    options: ShieldingOptions,
    /// 今の燃焼度の点の、参照の番号ごとの自己遮蔽と、それに渡した物質
    held: Vec<(usize, usize, FuelShielding, Vec<MaterialInput>)>,
    /// 自己遮蔽を解いた回数
    shieldings: usize,
}

impl<'a> BranchShielder<'a> {
    /// トラッキング `tracking` の格子の、領域ごとの物質の名前が `region_materials` の状態を、対称操作 `symmetry` と条件
    /// `options` で自己遮蔽する。`varied` はホウ素の軸が変える物質
    pub fn new(
        tracking: &'a Tracking,
        region_materials: &'a [String],
        varied: &'a [String],
        symmetry: &'a Detected,
        options: ShieldingOptions,
    ) -> Self {
        Self {
            tracking,
            reuse: ShieldingReuse::new(tracking),
            region_materials,
            varied,
            symmetry,
            options,
            held: Vec::new(),
            shieldings: 0,
        }
    }

    /// 自己遮蔽を解いた回数
    pub fn shieldings(&self) -> usize {
        self.shieldings
    }

    /// 物質 `state` の状態の、物質の一覧と [領域] のその番号。`reference` は参照から作るときの燃焼度の点と参照の番号で、その参照の
    /// 物質はまだ解いていなければ `nominal` で作る
    pub fn shield(
        &mut self,
        library: &Library,
        state: &[MaterialInput],
        reference: Option<(usize, usize)>,
        nominal: impl FnOnce() -> Vec<MaterialInput>,
    ) -> anyhow::Result<(Vec<Material>, Vec<usize>)> {
        let names = self.region_materials;
        let Some((point, index)) = reference else {
            self.shieldings += 1;
            let shielded = shield_fuel_reusing(
                library,
                &mut self.reuse,
                state,
                names,
                self.symmetry,
                &self.options,
            )?;
            return Ok((shielded.materials, shielded.region_material));
        };
        let lattice = &self.tracking.lattice;
        self.held.retain(|(p, ..)| *p == point);
        if !self.held.iter().any(|(_, i, ..)| *i == index) {
            let nominal = nominal();
            self.shieldings += 1;
            let shielded = shield_fuel_reusing(
                library,
                &mut self.reuse,
                &nominal,
                names,
                self.symmetry,
                &self.options,
            )?;
            // 作り直しの前提（自己遮蔽の物質が、領域の実効断面積と組成から作った群定数に SPH 因子を掛けたもの）を、公称の組成で
            // 作り直した物質が自己遮蔽の物質とビットで同じことで確かめる。崩れれば、状態の物質が自己遮蔽の作り方と食い違う
            let rebuilt = remix_varied(
                library,
                &shielded,
                lattice,
                names,
                self.varied,
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
                     nominal boron"
                );
            }
            self.held.push((point, index, shielded, nominal));
        }
        let (_, _, shielded, nominal) = self
            .held
            .iter()
            .find(|(_, i, ..)| *i == index)
            .ok_or_else(|| anyhow!("no reference shielding {index}"))?;
        // 参照から作り直すのはホウ素の軸の物質だけなので、ほかの物質は参照と同じこと（違えば、別の減速材か燃料の温度の参照）
        if nominal.len() != state.len() {
            bail!(
                "{} materials in the state, {} in the reference",
                state.len(),
                nominal.len()
            );
        }
        if let Some(m) = state
            .iter()
            .zip(nominal)
            .find(|(m, n)| m != n && !self.varied.contains(&m.name))
            .map(|(m, _)| m)
        {
            bail!(
                "material `{}` of the state differs from the reference at the nominal boron",
                m.name
            );
        }
        let mut unique = shielded.materials.clone();
        if state != nominal.as_slice() {
            let rebuilt = remix_varied(
                library,
                shielded,
                lattice,
                names,
                self.varied,
                nominal,
                state,
            )
            .context("remixing the materials of the boron axis")?;
            for (m, material) in rebuilt {
                unique[m] = material;
            }
        }
        Ok((unique, shielded.region_material.clone()))
    }
}

/// 公称の状態の物質 `nominal` で解いた自己遮蔽 `reference` の実効断面積と SPH 因子のまま、物質の名前が `varied` の領域の物質を、状態の
/// 物質 `state` の組成で作り直す（ADR-0201）。`region_materials` は [領域] の物質の名前で、`nominal` と `state` は同じ名前を同じ順に
/// 持ち、作り直す物質の温度は同じであること（実効断面積はその温度の値）。返すのは、作り直した `reference.materials` の番号と物質。
/// 自己遮蔽した物質は、自己遮蔽と同じく、核種の実効断面積から作った群定数を数密度で混ぜて SPH 因子を掛け、自己遮蔽しなかった物質は
/// `Material::mix`（無限希釈）で作る。公称の組成にない核種（公称のホウ素が 0 の B-10 など）は無限希釈の値
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
            if after.temperature != before.temperature {
                bail!(
                    "material `{name}` is at {} K in the state and {} K at the nominal boron",
                    after.temperature,
                    before.temperature
                );
            }
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
                    let shielded = match before.composition.iter().position(|(n, _)| n == nuclide) {
                        Some(i) => cross_sections[i].clone(),
                        None => {
                            library.shielded(nuclide, after.temperature, &vec![infinite; ng])?
                        }
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

/// 状態を軸の直積の順（最初の軸がいちばん遅く変わる）に解くときの、各状態の反復を始める解（ADR-0195）。状態 0 は始める解が
/// なく、ほかの状態は、添字が 0 でない最後の軸を 1 つ戻し、それより後の軸を 0 にした状態（その 1 つの軸だけが違う、前に解いた
/// 状態）の解から始める。直前に解いた状態は、その軸に加えて後の軸もすべて違う。
#[derive(Debug)]
pub struct Starts<T> {
    /// [軸] の格子点の数
    counts: Vec<usize>,
    /// 次に解く状態
    state: usize,
    /// [軸] の、その軸より後の軸の添字がすべて 0 の状態のうち、最後に解いたものの解
    saved: Vec<Option<Rc<T>>>,
}

impl<T> Starts<T> {
    /// `counts` は [軸] の格子点の数。状態 0 から順に解く
    pub fn new(counts: &[usize]) -> Self {
        Self {
            counts: counts.to_vec(),
            state: 0,
            saved: vec![None; counts.len()],
        }
    }

    /// 次に解く状態の反復を始める解
    pub fn start(&self) -> Option<&T> {
        self.moved(self.state)
            .and_then(|axis| self.saved[axis].as_deref())
    }

    /// 次に解く状態の解 `value` を残し、その次の状態に進む
    pub fn solved(&mut self, value: T) {
        let value = Rc::new(value);
        // 添字が 0 でない最後の軸より後の軸はすべて 0 なので、その軸と後の軸の解になる
        let first = self.moved(self.state).unwrap_or(0);
        for slot in &mut self.saved[first..] {
            *slot = Some(Rc::clone(&value));
        }
        self.state += 1;
    }

    /// 次に解く状態と、残した解（同じ解は 1 つ）と、[軸] の残した解の番号（チェックポイントに書く。ADR-0196）
    pub fn saved(&self) -> (usize, Vec<&T>, Vec<Option<usize>>) {
        let mut solutions: Vec<&Rc<T>> = Vec::new();
        let slots = self
            .saved
            .iter()
            .map(|slot| {
                slot.as_ref().map(|value| {
                    solutions
                        .iter()
                        .position(|s| Rc::ptr_eq(s, value))
                        .unwrap_or_else(|| {
                            solutions.push(value);
                            solutions.len() - 1
                        })
                })
            })
            .collect();
        (
            self.state,
            solutions.into_iter().map(|s| &**s).collect(),
            slots,
        )
    }

    /// `saved` で書いた、次に解く状態 `state`、解 `solutions`、[軸] の解の番号 `slots` から作り直す。`counts` は [軸] の格子点の
    /// 数。状態 0 より後で解の欠けた軸があれば（壊れたチェックポイント）、平らな初期値から黙って始めずにエラー
    pub fn restored(
        counts: &[usize],
        state: usize,
        solutions: Vec<T>,
        slots: &[Option<usize>],
    ) -> Result<Self, String> {
        if slots.len() != counts.len() {
            return Err(format!(
                "the saved starts have {} axes, the run has {}",
                slots.len(),
                counts.len()
            ));
        }
        let solutions: Vec<Rc<T>> = solutions.into_iter().map(Rc::new).collect();
        let saved = slots
            .iter()
            .map(|slot| {
                slot.map(|i| {
                    solutions
                        .get(i)
                        .cloned()
                        .ok_or_else(|| format!("start {i} of {} solutions", solutions.len()))
                })
                .transpose()
            })
            .collect::<Result<Vec<_>, _>>()?;
        let starts = Self {
            counts: counts.to_vec(),
            state,
            saved,
        };
        // 状態 0 を解くと、すべての軸が解を持つ（`solved`）。欠けた軸は、後の状態を平らな初期値から黙って始めさせる
        let present = starts.saved.iter().filter(|s| s.is_some()).count();
        match (state, present) {
            (0, 0) => Ok(starts),
            (0, _) => Err("the saved starts of state 0 must be empty".into()),
            (_, n) if n == counts.len() => Ok(starts),
            _ => Err(format!(
                "no start for state {state} on {} of {} axes",
                counts.len() - present,
                counts.len()
            )),
        }
    }

    /// 状態 `state` の、添字が 0 でない最後の軸
    fn moved(&self, state: usize) -> Option<usize> {
        state_indices(&self.counts, state)
            .iter()
            .rposition(|&i| i > 0)
    }
}

/// 軸の直積の順（最初の軸がいちばん遅く変わる）の状態 `state` の、[軸] の格子点の添字。`counts` は [軸] の格子点の数
pub fn state_indices(counts: &[usize], state: usize) -> Vec<usize> {
    let mut at = vec![0; counts.len()];
    let mut rest = state;
    for (k, &n) in counts.iter().enumerate().rev() {
        at[k] = rest % n;
        rest /= n;
    }
    at
}

#[cfg(test)]
mod tests;
