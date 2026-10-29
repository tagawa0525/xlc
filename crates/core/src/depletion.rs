//! 炉心の燃焼: ケースのステップに沿ってノードの燃焼度を進める（ADR-0190）。
//!
//! ノード n の初めの重金属の質量 Mₙ は、テーブルの `heavy_metal` (g/cm) × ノードの高さ。定格の出力密度 p (W/g) で EFPD を
//! ΔE 進めると、炉心の平均の燃焼度の増分は Δ = p ΔE / 1000 (MWd/kg) で、ノードの増分は Δ fₙ M / Mₙ（fₙ はノードの核分裂率
//! Σf φ V の割合、M = Σ Mₙ）。ステップは中点の形（ADR-0173 の CE/CM と同じ）: 初めの燃焼度で解いた割合で中点の燃焼度を、
//! 中点で解いた割合で終わりの燃焼度を求める。燃焼点（0 と各ステップの終わり）の値は、その点の燃焼度で解いたもの。状態が
//! ホウ素の探索を持てば、点と中点の解のすべてで k が目標になるホウ素濃度を探す（ADR-0191）。

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xlc_depletion::chain::Chain;
use xlc_format::msgpack::{self, MsgpackError};
use xlc_format::table::Table;

use crate::case::{CaseError, Solved, check_search, search_boron, solve_burnups};
use crate::input::{Boron, BoronSearch, CaseState, DepletionInput, Layout};
use crate::nem::{Options, fission_rates};

const FORMAT: &str = "xlc-core-history";

/// 炉心の燃焼の履歴（ADR-0190）。ノードの値は [集合体][軸方向のノード（下から）]
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct CoreHistory {
    /// 定格の出力密度 (W/g、炉心の平均、初めの重金属)
    pub specific_power: f64,
    /// すべてのステップで解いた状態
    pub state: CoreState,
    /// 軸方向のノードの境界 (cm、下から)
    pub axial_mesh: Vec<f64>,
    /// 集合体の ID（名前の順）
    pub assemblies: Vec<String>,
    /// ノードの初めの重金属の質量 (g)
    pub heavy_metal: Vec<Vec<f64>>,
    /// ミクロ燃焼（ADR-0216）で追う核種の名前（テーブルの順）。ミクロ燃焼しなければ `None`
    #[serde(skip_serializing_if = "Option::is_none")]
    pub nuclides: Option<Vec<String>>,
    pub points: Vec<CorePoint>,
}

/// 解いた炉心の状態（ケースの状態。燃焼度はノードごと）。`boron` と `search` は、入力と同じくどちらか一方だけを持つ
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CoreState {
    /// 与えたホウ素濃度 (ppm)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub boron: Option<f64>,
    /// 臨界ホウ素の探索（ADR-0191）
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub search: Option<BoronSearch>,
    /// ホウ素の中の B-10 の原子分率
    pub b10: f64,
    /// 減速材の温度 (K)
    pub moderator: f64,
    /// 燃料の温度 (K)
    pub fuel: f64,
    /// バンク → 位置（ステップ）
    pub banks: BTreeMap<String, u32>,
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
    /// その点を解いたホウ素濃度 (ppm。与えた値か、探索で求めた値)
    pub boron: f64,
    /// その点で解いた固有値計算の回数（与えた濃度なら 1、探索なら ADR-0145 の回数。範囲の両端を含み、下端で許容の中なら 1）
    pub solves: usize,
    /// 集合体ごとの燃焼度 (MWd/kg。質量の重み)
    pub assembly_burnups: Vec<f64>,
    /// ノードの燃焼度 (MWd/kg)
    pub burnups: Vec<Vec<f64>>,
    /// ノードの出力の割合（核分裂率の割合。炉心の和が 1）
    pub powers: Vec<Vec<f64>>,
    /// ミクロ燃焼のノードの追う核種の数密度 (1/(barn cm)。[集合体][軸方向のノード][核種]。核分裂のないノードは空)。
    /// ミクロ燃焼しなければ `None`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub densities: Option<Vec<Vec<Vec<f64>>>>,
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
    state: CoreState,
    axial_mesh: Vec<f64>,
    assemblies: Vec<String>,
    heavy_metal: Vec<Vec<f64>>,
    #[serde(default)]
    nuclides: Option<Vec<String>>,
    points: Vec<CorePoint>,
}

impl CoreHistory {
    /// フィールド名付きの map 形式の MessagePack: `{format: "xlc-core-history", specific_power, state: {boron | search: {boron:
    /// [下端, 上端], k, tolerance, max_solves}, b10, moderator, fuel, banks}, axial_mesh, assemblies,
    /// heavy_metal, points: [{efpd, burnup, k, boron, solves, assembly_burnups, burnups, powers}]}`
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
        if file.state.boron.is_some() == file.state.search.is_some() {
            return Err(CaseError::Depletion(
                "the history state must hold either boron or search".into(),
            ));
        }
        if let Some(search) = &file.state.search {
            search
                .check()
                .map_err(|m| CaseError::Depletion(format!("the history state: {m}")))?;
        }
        Ok(Self {
            specific_power: file.specific_power,
            state: file.state,
            axial_mesh: file.axial_mesh,
            assemblies: file.assemblies,
            heavy_metal: file.heavy_metal,
            nuclides: file.nuclides,
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

/// エラーに、燃焼のどの解かを添える。入力の組み合わせの誤り（どのファイルかを main が添える）と、場所をすでに持つエラーはそのまま
fn locate(error: CaseError, at: &dyn Fn() -> String) -> CaseError {
    match error {
        CaseError::Input(_) | CaseError::At { .. } => error,
        e => CaseError::At {
            at: at(),
            cause: Box::new(e),
        },
    }
}

/// 燃焼度の分布を 1 回解いた結果: 解、そのホウ素濃度 (ppm) と固有値計算の回数、位置ごとの核分裂率
struct Solve {
    solved: Solved,
    ppm: f64,
    solves: usize,
    rates: Vec<f64>,
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
    chain: Option<&Chain>,
    options: &Options,
    mut on_point: impl FnMut(&CorePoint),
) -> Result<CoreHistory, CaseError> {
    let _ = chain;
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
    if state.burnup.is_some() {
        return fail(format!("state {} gives a burnup", state.name));
    }
    // 長い計算を始める前に
    check_search(state, options)?;
    // 状態のバンクの位置で制御棒を入れた並び（ADR-0127）
    let layout = &layout.with_banks(&state.banks)?;
    let mesh = &layout.mesh;
    let [nx, ny, nz] = mesh.counts();
    let sites = nx * ny * nz;
    // 核分裂するテーブルは質量を持つこと（燃焼度の定義の質量がないので推測しない）。解く前に、ファイルを示して止める
    for (file, table) in layout.files.iter().zip(tables) {
        let fissile = table
            .states
            .iter()
            .any(|s| s.nu_fission.iter().any(|&v| v > 0.0));
        if fissile && table.heavy_metal.is_none() {
            return fail(format!(
                "{}: the table has fission but no heavy_metal; make it from a depletion history (ADR-0190)",
                file.display()
            ));
        }
    }
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
    // 燃焼度 `burnups` で解く（探索の状態は、k が目標になるホウ素濃度で）。`at` はエラーに添える解の場所
    let solve = |burnups: &[f64], at: &dyn Fn() -> String| -> Result<Solve, CaseError> {
        let solve_ppm = |ppm: f64| solve_burnups(layout, tables, &state.at(ppm), burnups, options);
        let solved = match &state.boron {
            Boron::Given(ppm) => solve_ppm(*ppm).map(|solved| (solved, *ppm, 1)),
            Boron::Search(search) => search_boron(search, solve_ppm)
                .map(|(solved, found)| (solved, found.ppm, found.solves)),
        };
        let (solved, ppm, solves) = solved.map_err(|e| locate(e, at))?;
        let rates = fission_rates(&solved.problem, &solved.solution);
        Ok(Solve {
            solved,
            ppm,
            solves,
            rates,
        })
    };
    let record = |efpd: f64, burnups: &[f64], solve: &Solve| {
        let Solve {
            solved,
            ppm,
            solves,
            rates,
        } = solve;
        let burnup_rows = by_assembly(burnups);
        let total: f64 = rates.iter().sum();
        let shares: Vec<f64> = rates.iter().map(|r| r / total).collect();
        CorePoint {
            efpd,
            burnup: weighted(burnups, &masses),
            k: solved.solution.k,
            boron: *ppm,
            solves: *solves,
            assembly_burnups: burnup_rows
                .iter()
                .zip(&heavy_metal)
                .map(|(b, m)| weighted(b, m))
                .collect(),
            burnups: burnup_rows,
            powers: by_assembly(&shares),
            densities: None,
        }
    };

    let mut burnups = vec![0.0; sites];
    let mut efpd = 0.0;
    let mut solve_at_point = solve(&burnups, &|| "the point at 0 EFPD".into())?;
    let mut points = vec![record(efpd, &burnups, &solve_at_point)];
    on_point(&points[0]);
    tracing::info!(
        efpd,
        k = points[0].k,
        boron = points[0].boron,
        solves = points[0].solves,
        "depletion point"
    );
    for step in &depletion.steps {
        let _span = tracing::info_span!("depletion step", efpd = step.efpd).entered();
        // 炉心の平均の燃焼度の増分 (MWd/kg。W/g × d = 1e-3 MWd/kg)。ステップの EFPD の差を足し合わせず、終わりの累積の EFPD の
        // 燃焼度と今の炉心の平均の差にする（ステップの丸めの誤差を積み上げない。積み上げると、テーブルの最後の点で終わるケースが
        // その点を越えて止まる）
        let target = depletion.specific_power * step.efpd * 1.0e-3;
        let increment = target - weighted(&burnups, &masses);
        burnups = midpoint(&burnups, &masses, increment, &solve_at_point.rates, |b| {
            let middle = solve(b, &|| {
                format!("the midpoint of the step to {} EFPD", step.efpd)
            })?;
            tracing::info!(
                k = middle.solved.solution.k,
                boron = middle.ppm,
                solves = middle.solves,
                "depletion midpoint"
            );
            Ok(middle.rates)
        })
        // 中点の解のエラーは場所を持つ。燃焼度の増分（`advance`）のエラーに、ステップを添える
        .map_err(|e| locate(e, &|| format!("the step to {} EFPD", step.efpd)))?;
        efpd = step.efpd;
        solve_at_point = solve(&burnups, &|| format!("the point at {efpd} EFPD"))?;
        let point = record(efpd, &burnups, &solve_at_point);
        on_point(&point);
        tracing::info!(
            efpd,
            k = point.k,
            boron = point.boron,
            solves = point.solves,
            "depletion point"
        );
        points.push(point);
    }
    let (boron, search) = match state.boron {
        Boron::Given(ppm) => (Some(ppm), None),
        Boron::Search(search) => (None, Some(search)),
    };
    Ok(CoreHistory {
        specific_power: depletion.specific_power,
        state: CoreState {
            boron,
            search,
            b10: state.b10,
            moderator: state.moderator,
            fuel: state.fuel,
            banks: state.banks.clone(),
        },
        axial_mesh: layout.axial_mesh().to_vec(),
        assemblies: placed.keys().map(|id| id.to_string()).collect(),
        heavy_metal,
        nuclides: None,
        points,
    })
}

#[cfg(test)]
mod tests;
