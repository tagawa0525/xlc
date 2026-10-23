//! 炉心の燃焼: ケースのステップに沿ってノードの燃焼度を進める（ADR-0190）。
//!
//! ノード n の初めの重金属の質量 Mₙ は、テーブルの `heavy_metal` (g/cm) × ノードの高さ。定格の出力密度 p (W/g) で EFPD を
//! ΔE 進めると、炉心の平均の燃焼度の増分は Δ = p ΔE / 1000 (MWd/kg) で、ノードの増分は Δ fₙ M / Mₙ（fₙ はノードの核分裂率
//! Σf φ V の割合、M = Σ Mₙ）。ステップは中点の形（ADR-0173 の CE/CM と同じ）: 初めの燃焼度で解いた割合で中点の燃焼度を、
//! 中点で解いた割合で終わりの燃焼度を求める。燃焼点（0 と各ステップの終わり）の値は、その点の燃焼度で解いたもの。

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xlc_format::msgpack::{self, MsgpackError};
use xlc_format::table::Table;

use crate::case::{CaseError, Solved, solve_burnups};
use crate::input::{Boron, CaseState, DepletionInput, Layout};
use crate::nem::{Options, fission_rates};

const FORMAT: &str = "xlc-core-history";

/// 炉心の燃焼の履歴（ADR-0190）。ノードの値は [集合体][軸方向のノード（下から）]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CoreHistory {
    /// 定格の出力密度 (W/g、炉心の平均、初めの重金属)
    pub specific_power: f64,
    /// 軸方向のノードの境界 (cm、下から)
    pub axial_mesh: Vec<f64>,
    /// 集合体の ID（名前の順）
    pub assemblies: Vec<String>,
    /// ノードの初めの重金属の質量 (g)
    pub heavy_metal: Vec<Vec<f64>>,
    pub points: Vec<CorePoint>,
}

/// 燃焼点の値（その点の燃焼度で解いたもの）
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CorePoint {
    /// 累積の EFPD
    pub efpd: f64,
    /// 炉心の平均の燃焼度 (MWd/kg。質量の重み)
    pub burnup: f64,
    pub k: f64,
    /// 集合体ごとの燃焼度 (MWd/kg。質量の重み)
    pub assembly_burnups: Vec<f64>,
    /// ノードの燃焼度 (MWd/kg)
    pub burnups: Vec<Vec<f64>>,
    /// ノードの出力の割合（核分裂率の割合。炉心の和が 1）
    pub powers: Vec<Vec<f64>>,
}

/// 書くファイルの形（`CoreHistory` に `format` を足したもの）
#[derive(Serialize)]
struct FileOut<'a> {
    format: &'static str,
    #[serde(flatten)]
    history: &'a CoreHistory,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct FileIn {
    format: String,
    specific_power: f64,
    axial_mesh: Vec<f64>,
    assemblies: Vec<String>,
    heavy_metal: Vec<Vec<f64>>,
    points: Vec<CorePoint>,
}

impl CoreHistory {
    /// フィールド名付きの map 形式の MessagePack: `{format: "xlc-core-history", specific_power, axial_mesh, assemblies,
    /// heavy_metal, points: [{efpd, burnup, k, assembly_burnups, burnups, powers}]}`
    pub fn to_msgpack(&self) -> Result<Vec<u8>, MsgpackError> {
        msgpack::to_named(&FileOut {
            format: FORMAT,
            history: self,
        })
    }

    /// `to_msgpack` の形を読む（形式の名前を確かめる）
    pub fn from_msgpack(bytes: &[u8]) -> Result<Self, CaseError> {
        let file: FileIn = msgpack::from_value(
            msgpack::read_value(bytes).map_err(|e| CaseError::Depletion(e.to_string()))?,
        )
        .map_err(|e| CaseError::Depletion(e.to_string()))?;
        if file.format != FORMAT {
            return Err(CaseError::Depletion(format!(
                "expected format {FORMAT}, found {}",
                file.format
            )));
        }
        Ok(Self {
            specific_power: file.specific_power,
            axial_mesh: file.axial_mesh,
            assemblies: file.assemblies,
            heavy_metal: file.heavy_metal,
            points: file.points,
        })
    }
}

/// 燃焼度 `burnups` (MWd/kg) を、炉心の平均の増分 `increment` (MWd/kg) だけ、核分裂率 `rates` の割合と質量 `masses` (g) で
/// 進める（ノードの増分は increment × fₙ × M / Mₙ）。核分裂するのに質量のないノードと、核分裂のない炉心はエラー
pub fn advance(
    burnups: &[f64],
    masses: &[f64],
    increment: f64,
    rates: &[f64],
) -> Result<Vec<f64>, CaseError> {
    let total: f64 = rates.iter().sum();
    if !(total > 0.0 && total.is_finite()) {
        return Err(CaseError::Depletion(format!(
            "the core fission rate is {total}"
        )));
    }
    let mass: f64 = masses.iter().sum();
    burnups
        .iter()
        .zip(masses)
        .zip(rates)
        .enumerate()
        .map(|(n, ((&b, &m), &r))| {
            if m > 0.0 {
                Ok(b + increment * (r / total) * (mass / m))
            } else if r > 0.0 {
                Err(CaseError::Depletion(format!(
                    "node {n} has fission but no heavy-metal mass; its table needs heavy_metal (ADR-0190)"
                )))
            } else {
                Ok(b)
            }
        })
        .collect()
}

/// 中点の形の 1 ステップ（ADR-0190）: 初めの燃焼度 `burnups` の核分裂率 `rates` で中点の燃焼度を求め、`solve`（燃焼度 →
/// 核分裂率）で中点の核分裂率を求めて、終わりの燃焼度を返す
pub fn midpoint(
    burnups: &[f64],
    masses: &[f64],
    increment: f64,
    rates: &[f64],
    mut solve: impl FnMut(&[f64]) -> Result<Vec<f64>, CaseError>,
) -> Result<Vec<f64>, CaseError> {
    let half = advance(burnups, masses, 0.5 * increment, rates)?;
    let middle = solve(&half)?;
    advance(burnups, masses, increment, &middle)
}

/// 質量 `m` の重みの平均の燃焼度
fn weighted(b: &[f64], m: &[f64]) -> f64 {
    b.iter().zip(m).map(|(b, m)| b * m).sum::<f64>() / m.iter().sum::<f64>()
}

/// ケースの燃焼（ADR-0190）: 状態 `state` で、ノードの燃焼度を 0 から `depletion` のステップに沿って進め、燃焼点ごとの値を
/// 返す。`layout` は燃焼のケースから作った並び（`radial-nodes` = 1）。点は解くたびに `on_point` にも渡す（後のステップで止まっても、
/// それまでの点を出せるように）
pub fn deplete(
    layout: &Layout,
    tables: &[Table],
    state: &CaseState,
    depletion: &DepletionInput,
    options: &Options,
    mut on_point: impl FnMut(&CorePoint),
) -> Result<CoreHistory, CaseError> {
    let fail = |m: String| Err(CaseError::Depletion(m));
    if layout.split != 1 {
        return fail(format!("radial-nodes = {} cannot deplete", layout.split));
    }
    if tables.len() != layout.files.len() {
        return fail(format!(
            "{} tables for {} few-group constant files",
            tables.len(),
            layout.files.len()
        ));
    }
    let ppm = match &state.boron {
        Boron::Given(ppm) => *ppm,
        Boron::Search(_) => return fail("boron searches are not taken in depletion steps".into()),
    };
    if state.burnup.is_some() {
        return fail(format!("state {} gives a burnup", state.name));
    }
    let point_state = state.at(ppm);
    // 状態のバンクの位置で制御棒を入れた並び（ADR-0127）
    let layout = &layout.with_banks(&state.banks)?;
    let mesh = &layout.mesh;
    let [nx, ny, nz] = mesh.counts();
    let sites = nx * ny * nz;
    // 位置ごとの初めの重金属の質量 (g)
    let masses: Vec<f64> = (0..sites)
        .map(|site| {
            let [_, _, z] = mesh.position(site);
            layout
                .site_file(site)
                .and_then(|f| tables[f].heavy_metal)
                .map_or(0.0, |m| m * mesh.dz[z])
        })
        .collect();
    // 集合体の ID（名前の順）と、その位置（x, y）。図の上の行が +y
    let rows = layout.ids.len();
    let mut placed: BTreeMap<&str, [usize; 2]> = BTreeMap::new();
    for (r, row) in layout.ids.iter().enumerate() {
        for (x, id) in row.iter().enumerate() {
            if let Some(id) = id {
                placed.insert(id, [x, rows - 1 - r]);
            }
        }
    }
    let by_assembly = |values: &[f64]| -> Vec<Vec<f64>> {
        placed
            .values()
            .map(|&[x, y]| (0..nz).map(|z| values[mesh.site([x, y, z])]).collect())
            .collect()
    };
    let heavy_metal = by_assembly(&masses);
    if let Some((id, _)) = placed
        .keys()
        .zip(&heavy_metal)
        .find(|(_, m)| m.iter().sum::<f64>() <= 0.0)
    {
        return fail(format!(
            "assembly {id} has no heavy-metal mass; its tables need heavy_metal (ADR-0190)"
        ));
    }
    let solve = |burnups: &[f64]| -> Result<(Solved, Vec<f64>), CaseError> {
        let solved = solve_burnups(layout, tables, &point_state, burnups, options)?;
        let rates = fission_rates(&solved.problem, &solved.solution);
        Ok((solved, rates))
    };
    let record = |efpd: f64, burnups: &[f64], solved: &Solved, rates: &[f64]| {
        let burnup_rows = by_assembly(burnups);
        let total: f64 = rates.iter().sum();
        let shares: Vec<f64> = rates.iter().map(|r| r / total).collect();
        CorePoint {
            efpd,
            burnup: weighted(burnups, &masses),
            k: solved.solution.k,
            assembly_burnups: burnup_rows
                .iter()
                .zip(&heavy_metal)
                .map(|(b, m)| weighted(b, m))
                .collect(),
            burnups: burnup_rows,
            powers: by_assembly(&shares),
        }
    };

    let mut burnups = vec![0.0; sites];
    let mut efpd = 0.0;
    let (mut solved, mut rates) = solve(&burnups)?;
    let mut points = vec![record(efpd, &burnups, &solved, &rates)];
    on_point(&points[0]);
    tracing::info!(efpd, k = solved.solution.k, "depletion point");
    for step in &depletion.steps {
        let _span = tracing::info_span!("depletion step", efpd = step.efpd).entered();
        // 炉心の平均の燃焼度の増分 (MWd/kg。W/g × d = 1e-3 MWd/kg)。ステップの EFPD の差を足し合わせず、終わりの累積の EFPD の
        // 燃焼度と今の炉心の平均の差にする（ステップの丸めの誤差を積み上げない。積み上げると、テーブルの最後の点で終わるケースが
        // その点を越えて止まる）
        let target = depletion.specific_power * step.efpd * 1.0e-3;
        let increment = target - weighted(&burnups, &masses);
        burnups = midpoint(&burnups, &masses, increment, &rates, |b| {
            solve(b).map(|(_, r)| r)
        })?;
        efpd = step.efpd;
        (solved, rates) = solve(&burnups)?;
        let point = record(efpd, &burnups, &solved, &rates);
        on_point(&point);
        points.push(point);
        tracing::info!(efpd, k = solved.solution.k, "depletion point");
    }
    Ok(CoreHistory {
        specific_power: depletion.specific_power,
        axial_mesh: layout.axial_mesh().to_vec(),
        assemblies: placed.keys().map(|id| id.to_string()).collect(),
        heavy_metal,
        points,
    })
}

#[cfg(test)]
mod tests;
