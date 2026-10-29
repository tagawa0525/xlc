//! 格子計算の分岐（ADR-0097、ADR-0200）。入力の `[branches]` の軸（ホウ素濃度、減速材の温度と密度、燃料の温度）から、
//! 軸の直積の順（最初の軸がいちばん遅く変わる）の物質の組を作り、そのうち解く状態を部分格子（`subgrids`）で選ぶ。
//!
//! ```toml
//! [branches]
//! boron = { ppm = [900, 1300, 1700], materials = ["water"] }   # ホウ酸水のホウ素濃度
//! # 温度と、ホウ酸水の密度（格子点ごと）。挙げた物質の温度を変え、ホウ酸水なら密度も変える
//! moderator = { temperature = [550, 565, 580], density = [0.76972, 0.74271, 0.71190], materials = ["water", "clad"] }
//! fuel = { temperature = [550, 565, 580], materials = ["fuel"] }
//! shielding = "nominal-boron" # 自己遮蔽（"nominal-boron" か "every-state"。省略すると nominal-boron。ADR-0201）
//! # 解く部分格子の軸の組（cut-HDMR。ADR-0200）。省略すると、[history] があれば 2 軸の組すべて、なければ 3 つの軸（直積。ADR-0097）
//! subgrids = [["boron", "moderator"], ["boron", "fuel"], ["moderator", "fuel"]]
//! ```
//!
//! 軸は boron、moderator、fuel の順に並べ、3つとも書く（core は状態を与えるとき3つの軸をすべて求める）。温度の軸が変える物質は名前で挙げ、物質の種類から
//! 推測しない。同じ物質を2つの温度の軸に挙げられない。
//!
//! 公称の状態は入力の物質の状態（燃焼の履歴の条件）で、軸ごとに、入力の物質をそのまま残す格子点（点が 1 つの軸はその点）。
//! 解く状態は、公称と違う軸の組が `subgrids` の組のどれかの部分集合になる状態（公称の状態、1 軸ずつの分岐、組の部分格子）。
//! 組のどれかが、点が 2 つ以上の軸をすべて持てば直積で、公称の状態は要らない。そうでなければ、入力の物質がどの軸でも格子点に
//! あること。
//!
//! 体積の割合で混ぜた物質（ADR-0128）は、軸が変える成分を持てばその軸に挙げ、持たなければ挙げない（混ぜた物質は
//! 状態ごとの成分から作る）。細かい群の媒質は燃料の温度の軸にだけ挙げ、分岐の状態ごとに媒質の同じ状態を使う。

use std::collections::BTreeMap;
use std::rc::Rc;

use anyhow::{Context, anyhow, bail};
use serde::Deserialize;
use xlc_format::fewgroup::FewGroup;
use xlc_format::library::{Library, LookupError};
use xlc_format::table::{Axis, Table, component};

use crate::depletion::{same_cross_sections, same_transfers};
use crate::geometry::Lattice;
use crate::input::{MaterialInput, MediumInput, Mixture};
use crate::material::Material;
use crate::shielding::{
    FuelShielding, ShieldingOptions, ShieldingReuse, shield_fuel_reusing, shielded_nuclide_with,
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
    /// [軸] の公称の格子点（入力の物質が格子点にない軸があれば `None`）
    nominal: Option<Vec<usize>>,
    /// [状態（軸の直積の順）] を解くか
    solved: Vec<bool>,
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

    /// 状態 `state`（軸の直積の順の番号。部分格子の状態だけを解くときも同じ番号）の自己遮蔽を参照から作るなら、その参照の
    /// `references` の番号（`every-state` と、ホウ素が 1 点で状態ごとに参照が違い再利用するものがないときは `None`）
    pub fn reference(&self, state: usize) -> Option<usize> {
        match self.shielding {
            BranchShielding::EveryState => None,
            BranchShielding::NominalBoron if self.references.len() == self.states.len() => None,
            // ホウ素は最初の軸（いちばん遅く変わる）なので、ほかの軸の添字は状態の番号をそれらの点の数の積で割った余り
            BranchShielding::NominalBoron => Some(state % self.references.len()),
        }
    }

    /// ホウ素の軸が変える物質
    pub fn varied(&self) -> &[String] {
        &self.varied
    }

    /// [軸] の公称の格子点（入力の物質の状態）
    pub fn nominal(&self) -> Option<&[usize]> {
        self.nominal.as_deref()
    }

    /// 公称の状態の、軸の直積の順の番号
    pub fn nominal_state(&self) -> Option<usize> {
        let counts: Vec<usize> = self.axes.iter().map(Axis::len).collect();
        self.nominal.as_ref().map(|n| {
            n.iter()
                .zip(&counts)
                .fold(0, |index, (&i, &c)| index * c + i)
        })
    }

    /// 状態 `state`（軸の直積の順）を解くか
    pub fn solved(&self, state: usize) -> bool {
        self.solved.get(state).copied().unwrap_or(false)
    }

    /// 軸の直積のすべての状態を解くか
    pub fn is_product(&self) -> bool {
        self.solved.iter().all(|&s| s)
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct BranchesEntry {
    boron: Option<BoronEntry>,
    moderator: Option<ModeratorEntry>,
    fuel: Option<FuelEntry>,
    shielding: Option<ShieldingEntry>,
    subgrids: Option<Vec<Vec<String>>>,
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

/// 解いた状態ごとの出力の部分ノードの少数群定数（[状態][出力][部分ノード]。分けていなければ部分ノードは1つ）を、出力
/// ごとのテーブルにまとめる。`points` は解いた状態の、軸の順の格子点の番号（テーブルの点）。
pub fn tables(
    axes: &[Axis],
    split: usize,
    points: &[Vec<usize>],
    per_state: Vec<Vec<Vec<FewGroup>>>,
) -> Result<Vec<Table>, String> {
    if per_state.len() != points.len() {
        return Err(format!(
            "{} states for {} points",
            per_state.len(),
            points.len()
        ));
    }
    let outputs = per_state.first().map_or(0, Vec::len);
    let mut tables: Vec<Table> = (0..outputs)
        .map(|_| Table {
            axes: axes.to_vec(),
            split,
            points: points.to_vec(),
            states: Vec::with_capacity(points.len() * split * split),
            heavy_metal: None,
            micro: None,
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

/// 軸の名前の組の bit（boron、moderator、fuel の順）
fn subgrid_mask(axis_set: &[String]) -> Result<u32, String> {
    let mut mask = 0;
    for name in axis_set {
        let bit = match name.as_str() {
            "boron" => 0,
            "moderator" => 1,
            "fuel" => 2,
            other => {
                return Err(format!(
                    "branches.subgrids: unknown axis `{other}` (boron, moderator or fuel)"
                ));
            }
        };
        if mask & (1 << bit) != 0 {
            return Err(format!(
                "branches.subgrids: `{name}` appears twice in a subgrid"
            ));
        }
        mask |= 1 << bit;
    }
    if mask == 0 {
        return Err("branches.subgrids: a subgrid needs axes".into());
    }
    Ok(mask)
}

/// 軸ごとの変更
enum Change {
    Boron(Vec<f64>, Vec<String>),
    Moderator(Vec<f64>, Vec<f64>, Vec<String>),
    Fuel(Vec<f64>, Vec<String>),
}

impl BranchesEntry {
    /// 軸を確かめ、軸の直積の順の物質の組（核種かホウ酸水で書いた物質）と、公称の状態と、解く状態を作る。`waters` はホウ酸水で
    /// 書いた物質の量、`mixtures` と `media` は軸に挙げる規則を確かめるだけ（ADR-0128）。
    /// `history` は燃焼の履歴からの分岐か（`subgrids` の既定を決める。ADR-0200）
    pub(crate) fn build(
        self,
        materials: &[MaterialInput],
        waters: &BTreeMap<String, Water>,
        mixtures: &[Mixture],
        media: &[MediumInput],
        history: bool,
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
            axes.push(Axis::boron(b.ppm.clone(), b10));
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
            axes.push(Axis::moderator(m.temperature.clone(), m.density.clone()));
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
            axes.push(Axis::fuel(f.temperature.clone()));
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
        // 公称の状態: 軸ごとに、その軸だけを変えても入力の物質が変わらない格子点（点が 1 つの軸はその点）
        let input = state_of(&[None, None, None])?;
        let mut per_axis = Vec::with_capacity(counts.len());
        for (k, &count) in counts.iter().enumerate() {
            if count == 1 {
                per_axis.push(Some(0));
                continue;
            }
            let mut found = Vec::new();
            for i in 0..count {
                let mut at = vec![None; counts.len()];
                at[k] = Some(i);
                if state_of(&at)? == input {
                    found.push(i);
                }
            }
            per_axis.push(if let [only] = found[..] {
                Some(only)
            } else {
                None
            });
        }
        // 解く部分格子の軸の組と、その部分集合（成分）
        let subgrids = match self.subgrids {
            Some(sets) if sets.is_empty() => {
                return Err("branches.subgrids: give the axis sets of the subgrids".into());
            }
            Some(sets) => sets
                .iter()
                .map(|set| subgrid_mask(set))
                .collect::<Result<Vec<_>, _>>()?,
            // 燃焼度の軸を持つテーブルは 2 軸の組すべて（ADR-0200）、M3 のテーブルは直積（ADR-0097）
            None if history => vec![0b011, 0b101, 0b110],
            None => vec![0b111],
        };
        let components = |mask: u32| subgrids.iter().any(|&s| s & mask == mask);
        let varying = counts
            .iter()
            .enumerate()
            .filter(|&(_, &c)| c > 1)
            .fold(0, |mask, (k, _)| mask | (1 << k));
        // 点が 2 つ以上の軸は、どれかの組が持つ（持たなければ、その軸の分岐を解かず、core が公称の値を黙って使う）
        let covered = subgrids.iter().fold(0, |mask, &s| mask | s);
        if let Some(axis) = (0..counts.len()).find(|&k| varying & !covered & (1 << k) != 0) {
            return Err(format!(
                "branches.subgrids: no subgrid holds the {} axis of {} points",
                axes[axis].name(),
                counts[axis]
            ));
        }
        let nominal: Option<Vec<usize>> = per_axis.iter().copied().collect();
        let solved = if components(varying) {
            vec![true; total]
        } else {
            let at = nominal.as_ref().ok_or_else(|| {
                let names: Vec<&str> = axes
                    .iter()
                    .zip(&per_axis)
                    .filter(|(_, n)| n.is_none())
                    .map(|(a, _)| a.name())
                    .collect();
                format!(
                    "branches.subgrids: the subgrids hang from the nominal state, the input materials, which are not \
                     at a point of the {names:?} axes (ADR-0200)"
                )
            })?;
            (0..total)
                .map(|index| components(component(&state_indices(&counts, index), at)))
                .collect()
        };
        let axes = match &nominal {
            Some(at) => axes
                .into_iter()
                .zip(at)
                .map(|(axis, &i)| {
                    let value = axis.points()[i];
                    axis.with_nominal(value)
                })
                .collect(),
            None => axes,
        };
        Ok(Branches {
            axes,
            states,
            shielding: match self.shielding {
                None | Some(ShieldingEntry::NominalBoron) => BranchShielding::NominalBoron,
                Some(ShieldingEntry::EveryState) => BranchShielding::EveryState,
            },
            references,
            varied,
            nominal,
            solved,
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

    /// 物質 `state` の状態の自己遮蔽。`reference` は参照から作るときの燃焼度の点と参照の番号で、その参照の物質はまだ解いていなければ
    /// `nominal` で作る。参照から作った状態は参照の自己遮蔽を持ち（作り直すのはホウ素の軸の物質だけで、燃料の実効断面積は参照と同じ）、
    /// 公称と違う状態は作り直した物質の一覧を添える
    pub fn shield_state(
        &mut self,
        library: &Library,
        state: &[MaterialInput],
        reference: Option<(usize, usize)>,
        nominal: impl FnOnce() -> Vec<MaterialInput>,
    ) -> anyhow::Result<StateShielding<'_>> {
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
            return Ok(StateShielding {
                fuel: Fuel::Solved(Box::new(shielded)),
                remixed: None,
            });
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
        let remixed = if state == nominal.as_slice() {
            None
        } else {
            let mut unique = shielded.materials.clone();
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
            Some(unique)
        };
        Ok(StateShielding {
            fuel: Fuel::Reference(shielded),
            remixed,
        })
    }
}

/// 分岐の 1 つの状態の自己遮蔽（`BranchShielder::shield_state`）: 燃料の自己遮蔽と、ホウ素の軸の物質を作り直したときの物質の一覧
#[derive(Debug)]
pub struct StateShielding<'a> {
    fuel: Fuel<'a>,
    remixed: Option<Vec<Material>>,
}

impl StateShielding<'_> {
    /// 物質の一覧
    pub fn materials(&self) -> &[Material] {
        self.remixed.as_deref().unwrap_or(&self.fuel.materials)
    }

    /// [領域] の `materials` の番号
    pub fn region_material(&self) -> &[usize] {
        &self.fuel.region_material
    }

    /// 状態の燃料の核種の実効断面積と SPH 因子を持つ自己遮蔽
    pub fn fuel(&self) -> &FuelShielding {
        &self.fuel
    }
}

/// 状態の燃料の自己遮蔽（`BranchShielder::shield_state`）: その状態で解いたものか、参照のもの
#[derive(Debug)]
pub enum Fuel<'a> {
    Solved(Box<FuelShielding>),
    Reference(&'a FuelShielding),
}

impl std::ops::Deref for Fuel<'_> {
    type Target = FuelShielding;

    fn deref(&self) -> &FuelShielding {
        match self {
            Self::Solved(s) => s,
            Self::Reference(s) => s,
        }
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
                    || !same_transfers(reference.region_transfers(r), reference.region_transfers(f))
                {
                    bail!(
                        "regions {f} and {r} share material {m} with different cross sections, transfers or SPH factors"
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
            let transfers = reference.region_transfers(r);
            if transfers.is_some_and(|t| t.len() != cross_sections.len()) {
                bail!("region {r}: the transfers do not match the shielded nuclides");
            }
            let nuclides = after
                .composition
                .iter()
                .map(|(nuclide, density)| {
                    let index = before.composition.iter().position(|(n, _)| n == nuclide);
                    let shielded = match index {
                        Some(i) => cross_sections[i].clone(),
                        None => {
                            library.shielded(nuclide, after.temperature, &vec![infinite; ng])?
                        }
                    };
                    // 公称の組成にない核種は無限希釈で、UFG の範囲の移行も持たない
                    let transfer = index.and_then(|i| transfers.map(|t| &t[i]));
                    Ok((
                        shielded_nuclide_with(
                            library,
                            nuclide,
                            after.temperature,
                            &shielded,
                            transfer,
                        )?,
                        *density,
                    ))
                })
                .collect::<anyhow::Result<Vec<_>>>()?;
            let mixed = Material::from_nuclides(&nuclides, ng)?.scaled(&sph[r])?;
            Ok((m, mixed))
        })
        .collect()
}

/// 状態を軸の直積の順（最初の軸がいちばん遅く変わる）に解くときの、各状態の反復を始める解（ADR-0195）。最初に解く状態は始める解
/// がなく、ほかの状態は、添字が 0 でない最後の軸を 1 つ戻し、それより後の軸を 0 にした状態（その 1 つの軸だけが違う、前に解いた
/// 状態）の解から始める。直前に解いた状態は、その軸に加えて後の軸もすべて違う。部分格子の状態だけを解き（ADR-0200）、その状態を
/// 解いていなければ、直前に解いた状態の解から始める（平らな初期値より近い）。
#[derive(Debug)]
pub struct Starts<T> {
    /// [軸] の格子点の数
    counts: Vec<usize>,
    /// 解く状態（軸の直積の順の番号。増える順）
    present: Vec<usize>,
    /// 次に解く状態の `present` の中の番号（解いた状態の数）
    next: usize,
    /// [軸] の、その軸より後の軸の添字がすべて 0 の解いた状態のうち、最後に解いたものの番号と解。最後の軸は直前に解いた状態
    saved: Vec<Option<(usize, Rc<T>)>>,
}

impl<T> Starts<T> {
    /// `counts` は [軸] の格子点の数。軸の直積のすべての状態を、状態 0 から順に解く
    pub fn new(counts: &[usize]) -> Self {
        Self::over(counts, (0..counts.iter().product()).collect())
    }

    /// `counts` は [軸] の格子点の数、`present` は解く状態（軸の直積の順の番号。増える順）
    pub fn over(counts: &[usize], present: Vec<usize>) -> Self {
        Self {
            counts: counts.to_vec(),
            present,
            next: 0,
            saved: vec![None; counts.len()],
        }
    }

    /// 次に解く状態の反復を始める解
    pub fn start(&self) -> Option<&T> {
        let state = *self.present.get(self.next)?;
        let at = state_indices(&self.counts, state);
        let axis = at.iter().rposition(|&i| i > 0)?;
        let mut back = at;
        back[axis] -= 1;
        back[axis + 1..].iter_mut().for_each(|i| *i = 0);
        let back = state_index(&self.counts, &back);
        match &self.saved[axis] {
            Some((solved, value)) if *solved == back => Some(value),
            _ => self.saved.last()?.as_ref().map(|(_, value)| &**value),
        }
    }

    /// 次に解く状態の解 `value` を残し、その次の状態に進む。解く状態をすべて解いた後に呼べばエラー（呼ぶ側の誤り）
    pub fn solved(&mut self, value: T) -> Result<(), String> {
        let value = Rc::new(value);
        let state = *self.present.get(self.next).ok_or_else(|| {
            format!(
                "solved more states than the {} the run solves",
                self.present.len()
            )
        })?;
        // 添字が 0 でない最後の軸より後の軸はすべて 0 なので、その軸と後の軸の解になる
        let first = self.moved(state).unwrap_or(0);
        for slot in &mut self.saved[first..] {
            *slot = Some((state, Rc::clone(&value)));
        }
        self.next += 1;
        Ok(())
    }

    /// 解いた状態の数と、残した解（同じ解は 1 つ）と、[軸] の残した解の番号（チェックポイントに書く。ADR-0196）
    pub fn saved(&self) -> (usize, Vec<&T>, Vec<Option<usize>>) {
        let mut solutions: Vec<&Rc<T>> = Vec::new();
        let slots = self
            .saved
            .iter()
            .map(|slot| {
                slot.as_ref().map(|(_, value)| {
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
            self.next,
            solutions.into_iter().map(|s| &**s).collect(),
            slots,
        )
    }

    /// `saved` で書いた、解いた状態の数 `next`、解 `solutions`、[軸] の解の番号 `slots` から作り直す。`counts` は [軸] の格子点の
    /// 数、`present` は解く状態。軸ごとの解の状態は解いた状態の並びから決まるので、その有無が書いた番号と合わなければ（壊れた
    /// チェックポイント）、近い解の代わりに平らな初期値や別の解から黙って始めずにエラー
    pub fn restored(
        counts: &[usize],
        present: Vec<usize>,
        next: usize,
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
        if next > present.len() {
            return Err(format!(
                "{next} states saved for a run of {}",
                present.len()
            ));
        }
        let mut starts = Self::over(counts, present);
        // 解いた状態の並びから決まる、軸ごとの解の状態
        let mut expected: Vec<Option<usize>> = vec![None; counts.len()];
        for &state in &starts.present[..next] {
            let first = starts.moved(state).unwrap_or(0);
            for slot in &mut expected[first..] {
                *slot = Some(state);
            }
        }
        if next == 0 && slots.iter().any(Option::is_some) {
            return Err("the saved starts of state 0 must be empty".into());
        }
        let missing = expected
            .iter()
            .zip(slots)
            .filter(|(e, s)| e.is_some() != s.is_some())
            .count();
        if missing > 0 {
            return Err(format!(
                "no start for state {next} on {missing} of {} axes",
                counts.len()
            ));
        }
        let solutions: Vec<Rc<T>> = solutions.into_iter().map(Rc::new).collect();
        starts.saved = expected
            .iter()
            .zip(slots)
            .map(|(state, slot)| match (state, slot) {
                (Some(state), Some(i)) => solutions
                    .get(*i)
                    .map(|value| Some((*state, Rc::clone(value))))
                    .ok_or_else(|| format!("start {i} of {} solutions", solutions.len())),
                _ => Ok(None),
            })
            .collect::<Result<Vec<_>, _>>()?;
        starts.next = next;
        Ok(starts)
    }

    /// 状態 `state` の、添字が 0 でない最後の軸
    fn moved(&self, state: usize) -> Option<usize> {
        state_indices(&self.counts, state)
            .iter()
            .rposition(|&i| i > 0)
    }
}

/// [軸] の格子点の添字 `at` の、軸の直積の順（最初の軸がいちばん遅く変わる）の状態の番号（`state_indices` の逆）
pub fn state_index(counts: &[usize], at: &[usize]) -> usize {
    at.iter()
        .zip(counts)
        .fold(0, |index, (&i, &c)| index * c + i)
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
