//! 格子計算の分岐（ADR-0097）。入力の `[branches]` の軸（ホウ素濃度、減速材の温度と密度、燃料の温度）から、
//! 軸の直積の順（最初の軸がいちばん遅く変わる）の物質の組を作る。
//!
//! ```toml
//! [branches]
//! boron = { ppm = [900, 1300, 1700], materials = ["water"] }   # ホウ酸水のホウ素濃度
//! # 温度と、ホウ酸水の密度（格子点ごと）。挙げた物質の温度を変え、ホウ酸水なら密度も変える
//! moderator = { temperature = [550, 565, 580], density = [0.76972, 0.74271, 0.71190], materials = ["water", "clad"] }
//! fuel = { temperature = [550, 565, 580], materials = ["fuel"] }
//! ```
//!
//! 軸は boron、moderator、fuel の順に並べ、3つとも書く（core は状態を与えるとき3つの軸をすべて求める）。温度の軸が変える物質は名前で挙げ、物質の種類から
//! 推測しない。同じ物質を2つの温度の軸に挙げられない。
//!
//! 体積の割合で混ぜた物質（ADR-0128）は、軸が変える成分を持てばその軸に挙げ、持たなければ挙げない（混ぜた物質は
//! 状態ごとの成分から作る）。細かい群の媒質は燃料の温度の軸にだけ挙げ、分岐の状態ごとに媒質の同じ状態を使う。

use std::collections::BTreeMap;
use std::rc::Rc;

use serde::Deserialize;
use xlc_format::fewgroup::FewGroup;
use xlc_format::table::{Axis, Table};

use crate::input::{MaterialInput, MediumInput, Mixture};
use crate::water::borated_water;

/// ホウ酸水で書いた物質の量（ADR-0106）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Water {
    pub density: f64,
    pub boron: f64,
    pub b10: f64,
}

/// 分岐の軸と、軸の直積の順の物質の組。
#[derive(Debug, Clone, PartialEq)]
pub struct Branches {
    axes: Vec<Axis>,
    states: Vec<Vec<MaterialInput>>,
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
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct BranchesEntry {
    boron: Option<BoronEntry>,
    moderator: Option<ModeratorEntry>,
    fuel: Option<FuelEntry>,
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

        let counts: Vec<usize> = axes.iter().map(|a| a.points().len()).collect();
        let total: usize = counts.iter().product();
        let mut states = Vec::with_capacity(total);
        for index in 0..total {
            let at = state_indices(&counts, index);
            let mut state = materials.to_vec();
            let mut state_waters = waters.clone();
            for (change, &i) in changes.iter().zip(&at) {
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
            states.push(state);
        }
        Ok(Branches { axes, states })
    }
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
    /// 数。次に解く状態の始める解がなければ（壊れたチェックポイント）、平らな初期値から黙って始めずにエラー
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
        if let Some(axis) = starts.moved(state)
            && starts.saved[axis].is_none()
        {
            return Err(format!("no start for state {state} on axis {axis}"));
        }
        Ok(starts)
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
