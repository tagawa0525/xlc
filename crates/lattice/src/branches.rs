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
//! 軸は boron、moderator、fuel の順に並べ、書かない軸は持たない。温度の軸が変える物質は名前で挙げ、物質の種類から
//! 推測しない。同じ物質を2つの温度の軸に挙げられない。

use std::collections::BTreeMap;

use serde::Deserialize;
use xlc_format::table::Axis;

use crate::input::MaterialInput;
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
    pub fn axes(&self) -> Vec<Axis> {
        self.axes.clone()
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

/// 挙げた物質が空でなく、重ならず、定義されていること。
fn check_materials(
    axis: &str,
    listed: &[String],
    materials: &[MaterialInput],
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
        if !materials.iter().any(|m| &m.name == name) {
            return Err(format!("branches.{axis}: unknown material `{name}`"));
        }
    }
    Ok(())
}

/// 軸ごとの変更
enum Change {
    Boron(Vec<f64>, Vec<String>),
    Moderator(Vec<f64>, Vec<f64>, Vec<String>),
    Fuel(Vec<f64>, Vec<String>),
}

impl BranchesEntry {
    /// 軸を確かめ、軸の直積の順の物質の組を作る。`waters` はホウ酸水で書いた物質の量。
    pub(crate) fn build(
        self,
        materials: &[MaterialInput],
        waters: &BTreeMap<String, Water>,
    ) -> Result<Branches, String> {
        let mut axes = Vec::new();
        let mut changes = Vec::new();
        if let Some(b) = self.boron {
            check_points("boron", &b.ppm, |p| (0.0..1.0e6).contains(&p))?;
            check_materials("boron", &b.materials, materials)?;
            let mut b10 = None;
            for name in &b.materials {
                let water = waters.get(name).ok_or_else(|| {
                    format!("branches.boron: material `{name}` is not written as water")
                })?;
                if b10.is_some_and(|f| f != water.b10) {
                    return Err("branches.boron: the materials have different b10".into());
                }
                b10 = Some(water.b10);
            }
            let b10 = b10.ok_or("branches.boron: list the materials it changes")?;
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
            check_materials("moderator", &m.materials, materials)?;
            if !m.materials.iter().any(|n| waters.contains_key(n)) {
                return Err("branches.moderator: no material written as water".into());
            }
            axes.push(Axis::Moderator {
                temperature: m.temperature.clone(),
                density: m.density.clone(),
            });
            changes.push(Change::Moderator(m.temperature, m.density, m.materials));
        }
        if let Some(f) = self.fuel {
            check_points("fuel", &f.temperature, |t| t > 0.0)?;
            check_materials("fuel", &f.materials, materials)?;
            if let Some(Change::Moderator(_, _, listed)) =
                changes.iter().find(|c| matches!(c, Change::Moderator(..)))
                && let Some(name) = f.materials.iter().find(|n| listed.contains(n))
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

        let counts: Vec<usize> = axes.iter().map(|a| a.points().len()).collect();
        let total: usize = counts.iter().product();
        let mut states = Vec::with_capacity(total);
        for index in 0..total {
            // 最初の軸がいちばん遅く変わる
            let mut at = vec![0; counts.len()];
            let mut rest = index;
            for (k, &n) in counts.iter().enumerate().rev() {
                at[k] = rest % n;
                rest /= n;
            }
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

#[cfg(test)]
mod tests;
