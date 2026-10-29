//! 炉心の燃焼: ケースのステップに沿ってノードの燃焼度を進める（ADR-0190）。
//!
//! ノード n の初めの重金属の質量 Mₙ は、テーブルの `heavy_metal` (g/cm) × ノードの高さ。定格の出力密度 p (W/g) で EFPD を
//! ΔE 進めると、炉心の平均の燃焼度の増分は Δ = p ΔE / 1000 (MWd/kg) で、ノードの増分は Δ fₙ M / Mₙ（fₙ はノードの核分裂率
//! Σf φ V の割合、M = Σ Mₙ）。ステップは中点の形（ADR-0173 の CE/CM と同じ）: 初めの燃焼度で解いた割合で中点の燃焼度を、
//! 中点で解いた割合で終わりの燃焼度を求める。燃焼点（0 と各ステップの終わり）の値は、その点の燃焼度で解いたもの。状態が
//! ホウ素の探索を持てば、点と中点の解のすべてで k が目標になるホウ素濃度を探す（ADR-0191）。
//!
//! ステップは状態と出力の割合を持つ（ADR-0218）。燃焼点は、そこから始まるステップの状態で解く（最後の点は最後のステップの
//! 状態）ので、ステップの初めの点、中点、修正子はどれもそのステップの状態の解になる。出力の割合はステップの時間（EFPD の増分 /
//! 割合）とノードの出力密度を決め、燃焼度の増分は変えない。崩壊だけのステップは、ミクロ燃焼の組成を崩壊の項だけで進める。

use std::collections::BTreeMap;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use xlc_depletion::CramOrder;
use xlc_depletion::chain::Chain;
use xlc_format::msgpack::{self, MsgpackError};
use xlc_format::table::{Axis, Micro, Table};

use crate::case::{
    CaseError, Compositions, Solved, check_search, search_boron, solve_compositions,
};
use crate::input::{Boron, BoronSearch, CaseState, DepletionInput, Layout, StepKind};
use crate::micro::MicroChain;
use crate::nem::{Options, fission_rates};
use crate::rods::mix_micro;
use crate::state::interpolate_micro;

const FORMAT: &str = "xlc-core-history";

/// 炉心の燃焼の履歴（ADR-0190）。ノードの値は [集合体][軸方向のノード（下から）]
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct CoreHistory {
    /// 定格の出力密度 (W/g、炉心の平均、初めの重金属)
    pub specific_power: f64,
    /// 解いた状態（名前 → 状態。ステップの使った状態だけ）
    pub states: BTreeMap<String, CoreState>,
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
    /// 累積の EFPD（崩壊だけのステップの前後で同じ）
    pub efpd: f64,
    /// 燃焼点 0 からの経過の日数（崩壊だけのステップの日数を含む）
    pub days: f64,
    /// その点で終わるステップの出力の割合（定格に対する比。崩壊だけのステップは 0）。燃焼点 0 は `None`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub power: Option<f64>,
    /// その点を解いた状態の名前（その点から始まるステップの状態。最後の点は最後のステップの状態。ADR-0218）
    pub state: String,
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
    states: BTreeMap<String, CoreState>,
    axial_mesh: Vec<f64>,
    assemblies: Vec<String>,
    heavy_metal: Vec<Vec<f64>>,
    #[serde(default)]
    nuclides: Option<Vec<String>>,
    points: Vec<CorePoint>,
}

impl CoreHistory {
    /// フィールド名付きの map 形式の MessagePack: `{format: "xlc-core-history", specific_power, states: {<名前>: {boron | search:
    /// {boron: [下端, 上端], k, tolerance, max_solves}, b10, moderator, fuel, banks}}, axial_mesh, assemblies, heavy_metal,
    /// nuclides?, points: [{efpd, days, power?, state, burnup, k, boron, solves, assembly_burnups, burnups, powers, densities?}]}`
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
        for (name, state) in &file.states {
            if state.boron.is_some() == state.search.is_some() {
                return Err(CaseError::Depletion(format!(
                    "the history state {name} must hold either boron or search"
                )));
            }
            if let Some(search) = &state.search {
                search
                    .check()
                    .map_err(|m| CaseError::Depletion(format!("the history state {name}: {m}")))?;
            }
        }
        // ミクロ燃焼の履歴は、どの点も [集合体][軸方向のノード] ごとに追う核種の数の有限な数密度（核分裂のないノードは空）を持つ
        let layers = file.axial_mesh.len().saturating_sub(1);
        for (k, point) in file.points.iter().enumerate() {
            let fits = match (&file.nuclides, &point.densities) {
                (None, None) => true,
                (Some(names), Some(rows)) => {
                    rows.len() == file.assemblies.len()
                        && rows.iter().all(|r| r.len() == layers)
                        && rows.iter().flatten().all(|n| {
                            (n.is_empty() || n.len() == names.len())
                                && n.iter().all(|v| v.is_finite())
                        })
                }
                _ => false,
            };
            if !fits {
                return Err(CaseError::Depletion(format!(
                    "point {k}: the densities must match the micro nuclides of the history"
                )));
            }
            if !file.states.contains_key(&point.state) {
                return Err(CaseError::Depletion(format!(
                    "point {k}: the history has no state {}",
                    point.state
                )));
            }
            // 燃焼点 0 だけが、そこで終わるステップを持たない
            if point.power.is_some() == (k == 0) {
                return Err(CaseError::Depletion(format!(
                    "point {k}: only the first point has no step power"
                )));
            }
            // 出力の割合は 0 以上の有限の値で、0 は EFPD の進まない崩壊だけのステップだけ。EFPD と日数は減らない（ADR-0218）
            if let (Some(power), Some(before)) =
                (point.power, k.checked_sub(1).map(|j| &file.points[j]))
                && !(power.is_finite()
                    && power >= 0.0
                    && (power == 0.0) == (point.efpd == before.efpd)
                    && point.efpd >= before.efpd
                    && point.days >= before.days)
            {
                return Err(CaseError::Depletion(format!(
                    "point {k}: the step power {power}, EFPD {} and days {} do not follow the point before \
                     ({} EFPD, {} days)",
                    point.efpd, point.days, before.efpd, before.days
                )));
            }
        }
        Ok(Self {
            specific_power: file.specific_power,
            states: file.states,
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

/// 燃焼度の分布を 1 回解いた結果: 解いた状態（`states` の番号）、解、そのホウ素濃度 (ppm) と固有値計算の回数、位置ごとの核分裂率
struct Solve {
    state: usize,
    solved: Solved,
    ppm: f64,
    solves: usize,
    rates: Vec<f64>,
}

/// 質量 `m` の重みの平均の燃焼度
fn weighted(b: &[f64], m: &[f64]) -> f64 {
    b.iter().zip(m).map(|(b, m)| b * m).sum::<f64>() / m.iter().sum::<f64>()
}

/// 状態 `state` の、燃焼の履歴の形
fn core_state(state: &CaseState) -> CoreState {
    let (boron, search) = match state.boron {
        Boron::Given(ppm) => (Some(ppm), None),
        Boron::Search(search) => (None, Some(search)),
    };
    CoreState {
        boron,
        search,
        b10: state.b10,
        moderator: state.moderator,
        fuel: state.fuel,
        banks: state.banks.clone(),
    }
}

/// ケースの燃焼（ADR-0190、ADR-0218）: 状態 `states`（ケースの名前つきの状態。ステップが番号で指す）で、ノードの燃焼度を 0 から
/// `depletion` のステップに沿って進め、燃焼点ごとの値を返す。`layout` は燃焼のケースから作った並び（`radial-nodes` = 1）。点は
/// 解くたびに `on_point` にも渡す（後のステップで止まっても、それまでの点を出せるように）
pub fn deplete(
    layout: &Layout,
    tables: &[Table],
    states: &[CaseState],
    depletion: &DepletionInput,
    chain: Option<&Chain>,
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
    if states.is_empty() || depletion.steps.is_empty() {
        return fail(format!(
            "the depletion needs a state and a step, got {} states and {} steps",
            states.len(),
            depletion.steps.len()
        ));
    }
    if let Some(step) = depletion.steps.iter().find(|s| s.state >= states.len()) {
        return fail(format!(
            "a step uses state {} of {} states",
            step.state,
            states.len()
        ));
    }
    // 燃焼点 i を解く状態: そこから始まるステップの状態。最後の点は最後のステップの状態（ADR-0218 の案 P2）
    let point_state = |i: usize| {
        depletion
            .steps
            .get(i)
            .or(depletion.steps.last())
            .map_or(0, |s| s.state)
    };
    // ステップの使う状態（番号の順）
    let used: Vec<usize> = (0..states.len())
        .filter(|&k| depletion.steps.iter().any(|s| s.state == k) || point_state(0) == k)
        .collect();
    for &k in &used {
        let state = &states[k];
        if state.burnup.is_some() {
            return fail(format!("state {} gives a burnup", state.name));
        }
        // 長い計算を始める前に
        check_search(state, options)?;
    }
    // 状態のバンクの位置で制御棒を入れた並び（ADR-0127）。使わない状態は並びを作らない（`None`）
    let layouts: Vec<Option<Layout>> = (0..states.len())
        .map(|k| {
            used.contains(&k)
                .then(|| layout.with_banks(&states[k].banks))
                .transpose()
        })
        .collect::<Result<_, _>>()?;
    let first = layouts[point_state(0)]
        .as_ref()
        .ok_or_else(|| CaseError::Depletion("the first point has no state".into()))?;
    // 組成と燃焼度は位置ごとに持つので、どの状態の並びも同じファイルと位置であること
    if let Some((k, _)) = layouts.iter().enumerate().find(|(_, l)| {
        l.as_ref()
            .is_some_and(|l| l.files != first.files || l.mesh.counts() != first.mesh.counts())
    }) {
        return fail(format!(
            "state {} lays out other files or sites than state {}",
            states[k].name,
            states[point_state(0)].name
        ));
    }
    let mesh = &first.mesh;
    let [nx, ny, nz] = mesh.counts();
    let sites = nx * ny * nz;
    // 核分裂するテーブルは質量を持つこと（燃焼度の定義の質量がないので推測しない）。解く前に、ファイルを示して止める
    for (file, table) in first.files.iter().zip(tables) {
        if fissile(table) && table.heavy_metal.is_none() {
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
            first
                .site_file(site)
                .and_then(|f| tables[f].heavy_metal)
                .map_or(0.0, |m| m * mesh.dz[z])
        })
        .collect();
    // 集合体の ID（名前の順）と、その位置（x, y）。図の上の行が +y
    let rows = first.ids.len();
    let mut placed: BTreeMap<&str, [usize; 2]> = BTreeMap::new();
    for (r, row) in first.ids.iter().enumerate() {
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
    // ミクロ燃焼（ADR-0216）のノードと、初めの組成（テーブルの燃焼度 0 の数密度。補助の核種は 0）
    let micro = match (&depletion.micro, chain) {
        (Some(_), Some(chain)) => Some(MicroNodes::new(&layouts, tables, chain)?),
        (None, None) => None,
        (Some(_), None) => {
            return fail("the case asks for micro depletion, but no chain is given".into());
        }
        (None, Some(_)) => {
            return fail("a chain is given, but the case asks for no micro depletion".into());
        }
    };
    if micro.is_none()
        && depletion
            .steps
            .iter()
            .any(|s| matches!(s.kind, StepKind::Decay { .. }))
    {
        return fail("a decay step needs micro depletion (ADR-0218)".into());
    }
    let mut densities = micro
        .as_ref()
        .map(|m| m.initial(point_state(0), &layouts, tables))
        .transpose()?;
    if let Some(m) = &micro {
        tracing::info!(
            tracked = m.chain.nuclides.len(),
            auxiliary = m.chain.auxiliary(),
            nodes = m.uses[point_state(0)].iter().flatten().count(),
            "micro depletion"
        );
    }
    // 状態 `k` の燃焼度 `burnups`（ミクロ燃焼なら組成 `n` も）で解く（探索の状態は、k が目標になるホウ素濃度で）。`at` は
    // エラーに添える解の場所
    let solve = |k: usize,
                 burnups: &[f64],
                 n: Option<&[Option<Vec<f64>>]>,
                 at: &dyn Fn() -> String|
     -> Result<Solve, CaseError> {
        let state = &states[k];
        let layout = layouts[k]
            .as_ref()
            .ok_or_else(|| CaseError::Depletion(format!("state {} is not laid out", state.name)))?;
        let compositions = micro.as_ref().zip(n).map(|(m, densities)| Compositions {
            chain: &m.chain,
            densities,
        });
        let solve_ppm = |ppm: f64| {
            solve_compositions(
                layout,
                tables,
                &state.at(ppm),
                burnups,
                compositions.as_ref(),
                options,
            )
        };
        let solved = match &state.boron {
            Boron::Given(ppm) => solve_ppm(*ppm).map(|solved| (solved, *ppm, 1)),
            Boron::Search(search) => search_boron(search, solve_ppm)
                .map(|(solved, found)| (solved, found.ppm, found.solves)),
        };
        let (solved, ppm, solves) = solved.map_err(|e| locate(e, at))?;
        let rates = fission_rates(&solved.problem, &solved.solution);
        Ok(Solve {
            state: k,
            solved,
            ppm,
            solves,
            rates,
        })
    };
    // 解 `solve`（燃焼度 `burnups`、組成 `n`）の反応率で、組成 `from` を時間 `t` (s) だけ進める。ノードの出力密度は出力の割合
    // `power` × 定格の出力密度 × 解の核分裂率の割合 × 炉心の質量 / ノードの体積（燃焼度を進める `advance` と同じエネルギーを、
    // 割合で伸びた時間で出す。ADR-0218）
    let mass: f64 = masses.iter().sum();
    let deplete_nodes = |m: &MicroNodes,
                         solve: &Solve,
                         power: f64,
                         burnups: &[f64],
                         n: &[Option<Vec<f64>>],
                         from: &[Option<Vec<f64>>],
                         t: f64,
                         at: &dyn Fn() -> String|
     -> Result<Vec<Option<Vec<f64>>>, CaseError> {
        let total: f64 = solve.rates.iter().sum();
        let state = states[solve.state].at(solve.ppm);
        let (Some(layout), uses) = (&layouts[solve.state], &m.uses[solve.state]) else {
            return Err(CaseError::Depletion(format!(
                "state {} is not laid out",
                states[solve.state].name
            )));
        };
        let start = std::time::Instant::now();
        let next = (0..sites)
            .into_par_iter()
            .map(|site| {
                let (Some(node_use), Some(n), Some(from)) = (uses[site], &n[site], &from[site])
                else {
                    return Ok(None);
                };
                let [x, y, z] = mesh.position(site);
                let volume = mesh.dx[x] * mesh.dy[y] * mesh.dz[z];
                let power_density =
                    power * depletion.specific_power * solve.rates[site] / total * mass / volume;
                let s = crate::state::State {
                    burnup: Some(burnups[site]),
                    ..state
                };
                let micro_of = |table_use: usize| {
                    let t_use = &layout.tables[table_use];
                    interpolate_micro(&tables[t_use.file], &s, t_use.rod.as_deref())
                        .map_err(|e| format!("{}: {e}", layout.files[t_use.file].display()))
                };
                let sigma = match node_use {
                    NodeUse::Uniform(table_use) => micro_of(table_use)?,
                    // 区間の微視的断面積を、ノードの少数群定数を混ぜた重みで混ぜる（ADR-0217）
                    NodeUse::Partial(index) => {
                        let micro = layout.partial[index]
                            .segments
                            .iter()
                            .map(|&(table_use, _)| micro_of(table_use))
                            .collect::<Result<Vec<_>, String>>()?;
                        let weights = solve.solved.segment_weights.get(index).ok_or_else(|| {
                            format!(
                                "site {site}: the solve has no weights of the partially rodded node"
                            )
                        })?;
                        mix_micro(&micro, weights, &solve.solved.solution.flux[site])
                            .map_err(|e| format!("site {site}: {e}"))?
                    }
                };
                let rates =
                    m.chain
                        .rates(&sigma, &solve.solved.solution.flux[site], n, power_density)?;
                let a = m.chain.matrix(&rates)?;
                let next =
                    xlc_depletion::deplete(&a, from, t, CRAM_ORDER).map_err(|e| e.to_string())?;
                check_negative(&m.chain, &next)?;
                Ok(Some(next))
            })
            .collect::<Result<Vec<_>, String>>()
            .map_err(|e| locate(CaseError::Depletion(e), at))?;
        tracing::debug!(
            seconds = start.elapsed().as_secs_f64(),
            "micro depletion of the nodes"
        );
        Ok(next)
    };
    // 崩壊だけのステップ: 組成 `n` を、崩壊の項だけの行列で時間 `t` (s) だけ進める（ADR-0218）
    let decay_nodes = |m: &MicroNodes,
                       n: &[Option<Vec<f64>>],
                       t: f64,
                       at: &dyn Fn() -> String|
     -> Result<Vec<Option<Vec<f64>>>, CaseError> {
        let a = m
            .chain
            .decay_matrix()
            .map_err(|e| locate(CaseError::Depletion(e), at))?;
        n.par_iter()
            .map(|n| {
                n.as_ref()
                    .map(|n| {
                        let next = xlc_depletion::deplete(&a, n, t, CRAM_ORDER)
                            .map_err(|e| e.to_string())?;
                        check_negative(&m.chain, &next)?;
                        Ok(next)
                    })
                    .transpose()
            })
            .collect::<Result<Vec<_>, String>>()
            .map_err(|e| locate(CaseError::Depletion(e), at))
    };
    let record = |efpd: f64,
                  days: f64,
                  power: Option<f64>,
                  burnups: &[f64],
                  solve: &Solve,
                  n: Option<&[Option<Vec<f64>>]>| {
        let Solve {
            state,
            solved,
            ppm,
            solves,
            rates,
        } = solve;
        let burnup_rows = by_assembly(burnups);
        let total: f64 = rates.iter().sum();
        let shares: Vec<f64> = rates.iter().map(|r| r / total).collect();
        let densities = micro.as_ref().zip(n).map(|(m, n)| {
            placed
                .values()
                .map(|&[x, y]| {
                    (0..nz)
                        .map(|z| {
                            n[mesh.site([x, y, z])]
                                .as_ref()
                                .map_or_else(Vec::new, |v| m.chain.tracked(v))
                        })
                        .collect()
                })
                .collect()
        });
        CorePoint {
            efpd,
            days,
            power,
            state: states[*state].name.clone(),
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
            densities,
        }
    };
    let report = |point: &CorePoint| {
        tracing::info!(
            efpd = point.efpd,
            days = point.days,
            state = point.state,
            k = point.k,
            boron = point.boron,
            solves = point.solves,
            "depletion point"
        );
    };

    let mut burnups = vec![0.0; sites];
    let mut efpd = 0.0;
    let mut days = 0.0;
    let mut solve_at_point = solve(point_state(0), &burnups, densities.as_deref(), &|| {
        "the point at 0 EFPD".into()
    })?;
    let mut points = vec![record(
        efpd,
        days,
        None,
        &burnups,
        &solve_at_point,
        densities.as_deref(),
    )];
    on_point(&points[0]);
    report(&points[0]);
    for (i, step) in depletion.steps.iter().enumerate() {
        let (power, at_point) = match step.kind {
            StepKind::Decay { days: d } => {
                let _span = tracing::info_span!("decay step", efpd, days = d).entered();
                let at_step = || format!("the decay step of {d} days after {efpd} EFPD");
                if let (Some(m), Some(n)) = (&micro, &densities) {
                    densities = Some(decay_nodes(m, n, d * DAY, &at_step)?);
                }
                days += d;
                (
                    0.0,
                    format!("the point after the decay to {days} days at {efpd} EFPD"),
                )
            }
            StepKind::Power {
                efpd: end,
                power: fraction,
            } => {
                let _span = tracing::info_span!("depletion step", efpd = end).entered();
                // 炉心の平均の燃焼度の増分 (MWd/kg。W/g × d = 1e-3 MWd/kg)。ステップの EFPD の差を足し合わせず、終わりの累積の
                // EFPD の燃焼度と今の炉心の平均の差にする（ステップの丸めの誤差を積み上げない。積み上げると、テーブルの最後の点で
                // 終わるケースがその点を越えて止まる）。出力の割合によらない（EFPD は割合を含む。ADR-0218）
                let target = depletion.specific_power * end * 1.0e-3;
                let increment = target - weighted(&burnups, &masses);
                // ステップの時間 (s)。EFPD の差を出力の割合で割った日数
                let step_days = (end - efpd) / fraction;
                let time = step_days * DAY;
                let at_step = || format!("the step to {end} EFPD");
                // 中点の解と、その燃焼度と組成（ミクロ燃焼の修正子の反応率に使う。ADR-0173 の CE/CM）
                let mut middle_of_step = None;
                let start = burnups.clone();
                burnups = midpoint(&burnups, &masses, increment, &solve_at_point.rates, |b| {
                    // 予測子: 初めの解の反応率で、中点の組成へ
                    let half = match (&micro, &densities) {
                        (Some(m), Some(n)) => Some(deplete_nodes(
                            m,
                            &solve_at_point,
                            fraction,
                            &start,
                            n,
                            n,
                            0.5 * time,
                            &at_step,
                        )?),
                        _ => None,
                    };
                    let middle = solve(step.state, b, half.as_deref(), &|| {
                        format!("the midpoint of the step to {end} EFPD")
                    })?;
                    tracing::info!(
                        k = middle.solved.solution.k,
                        boron = middle.ppm,
                        solves = middle.solves,
                        "depletion midpoint"
                    );
                    let rates = middle.rates.clone();
                    middle_of_step = Some((middle, b.to_vec(), half));
                    Ok(rates)
                })
                // 中点の解のエラーは場所を持つ。燃焼度の増分（`advance`）のエラーに、ステップを添える
                .map_err(|e| locate(e, &at_step))?;
                // 修正子: 中点の解の反応率（中点の燃焼度の σ、中点の組成で規格化）で、初めの組成から終わりへ
                if let (Some(m), Some(n), Some((middle, half_burnups, Some(half)))) =
                    (&micro, &densities, &middle_of_step)
                {
                    densities = Some(deplete_nodes(
                        m,
                        middle,
                        fraction,
                        half_burnups,
                        half,
                        n,
                        time,
                        &at_step,
                    )?);
                }
                efpd = end;
                days += step_days;
                (fraction, format!("the point at {efpd} EFPD"))
            }
        };
        solve_at_point = solve(point_state(i + 1), &burnups, densities.as_deref(), &|| {
            at_point.clone()
        })?;
        let point = record(
            efpd,
            days,
            Some(power),
            &burnups,
            &solve_at_point,
            densities.as_deref(),
        );
        on_point(&point);
        report(&point);
        points.push(point);
    }
    Ok(CoreHistory {
        specific_power: depletion.specific_power,
        states: used
            .iter()
            .map(|&k| (states[k].name.clone(), core_state(&states[k])))
            .collect(),
        axial_mesh: first.axial_mesh().to_vec(),
        assemblies: placed.keys().map(|id| id.to_string()).collect(),
        heavy_metal,
        nuclides: micro.map(|m| m.chain.nuclides),
        points,
    })
}

/// 1 日 (s)
const DAY: f64 = 86_400.0;

/// 行列指数の次数（格子計算の燃焼と同じ CRAM48。ADR-0173）
const CRAM_ORDER: CramOrder = CramOrder::Order48;

/// CRAM の近似の誤差で説明できる負の数密度の、ノードの最大の数密度に対する大きさ（格子計算の燃焼の `NEGATIVE_TOLERANCE` と同じ。
/// ADR-0161、ADR-0173）。これより大きな負の値はエラーにする
const NEGATIVE_TOLERANCE: f64 = 1.0e-12;

/// 負の数密度が CRAM の近似の誤差の水準を超えないこと
fn check_negative(chain: &MicroChain, n: &[f64]) -> Result<(), String> {
    let largest = n.iter().fold(0.0f64, |a, &b| a.max(b));
    match n
        .iter()
        .zip(&chain.solved)
        .find(|&(&v, _)| v < -NEGATIVE_TOLERANCE * largest)
    {
        Some((v, name)) => Err(format!(
            "{name} has a negative density {v} beyond the CRAM error"
        )),
        None => Ok(()),
    }
}

/// ミクロ燃焼するノード（ADR-0216）
struct MicroNodes {
    chain: MicroChain,
    /// [状態][直方体の位置] のテーブルの使い方（ミクロ燃焼する、核分裂のある位置だけ）。並びを作らない状態は空
    uses: Vec<Vec<Option<NodeUse>>>,
}

/// ミクロ燃焼するノードの微視的断面積の引き方
#[derive(Debug, Clone, Copy)]
enum NodeUse {
    /// 一様なノード: テーブルの使い方の番号
    Uniform(usize),
    /// 部分挿入のノード: `Layout::partial` の番号。区間のテーブルの使い方の微視的断面積を混ぜる（ADR-0217）
    Partial(usize),
}

/// テーブルが核分裂するか
fn fissile(table: &Table) -> bool {
    table
        .states
        .iter()
        .any(|s| s.nu_fission.iter().any(|&v| v > 0.0))
}

/// 並び `layout` の位置ごとのミクロ燃焼の使い方。核分裂するテーブルは、燃焼度で引く使い方で `micro` を持つこと。使った
/// テーブルの `micro` を `files` に集める
fn node_uses<'a>(
    layout: &Layout,
    tables: &'a [Table],
    files: &mut BTreeMap<usize, &'a Micro>,
) -> Result<Vec<Option<NodeUse>>, CaseError> {
    let fail = |m: String| Err(CaseError::Depletion(m));
    let mut uses = vec![None; layout.node_material.len()];
    // 位置 `site` のテーブルの使い方 `table` が核分裂するか。核分裂するなら、燃焼度で引く使い方で、テーブルが `micro` を持つこと
    let mut fissile_use = |site: usize, table: usize| -> Result<bool, CaseError> {
        let file = layout.tables[table].file;
        if !fissile(&tables[file]) {
            return Ok(false);
        }
        // 燃焼度で引かない使い方（反射体の位置など）は、補正も燃焼もしないので、核分裂のあるテーブルを置けない
        if layout.tables[table].node.is_none() {
            return Err(CaseError::Depletion(format!(
                "site {site} uses the fissile table {} but is not a depleted node, so the micro depletion would \
                 neither deplete nor correct it",
                layout.files[file].display()
            )));
        }
        let Some(micro) = &tables[file].micro else {
            return Err(CaseError::Depletion(format!(
                "{}: the table has no micro data for the micro depletion; make it with [history] micro \
                 (ADR-0216)",
                layout.files[file].display()
            )));
        };
        files.insert(file, micro);
        Ok(true)
    };
    for (site, material) in layout.node_material.iter().enumerate() {
        let Some(m) = *material else { continue };
        uses[site] = match layout.materials.get(m) {
            Some(material) => {
                fissile_use(site, material.table)?.then_some(NodeUse::Uniform(material.table))
            }
            None => {
                // 部分挿入のノード: 区間のどれかが核分裂すれば、区間の微視的断面積を混ぜて燃焼させる（ADR-0217）。区間は同じ
                // 層のテーブルで、初めの組成と公称の数密度は区間によらない
                let index = m - layout.materials.len();
                let segments = &layout.partial[index].segments;
                if segments
                    .iter()
                    .any(|&(t, _)| layout.tables[t].file != layout.tables[segments[0].0].file)
                {
                    return fail(format!(
                        "site {site}: the intervals of a partially rodded node use different tables"
                    ));
                }
                let mut any = false;
                for &(table, _) in &layout.partial[index].segments {
                    any |= fissile_use(site, table)?;
                }
                any.then_some(NodeUse::Partial(index))
            }
        };
    }
    Ok(uses)
}

impl MicroNodes {
    /// 状態ごとの並び `layouts`（並びを作らない状態は `None`）の核分裂のある位置をミクロ燃焼する。そのテーブルは `micro` を持ち、
    /// 追う核種はどのテーブルも同じで、チェーンの反応と合うこと。組成は位置ごとに持つので、ミクロ燃焼する位置はどの状態でも同じこと
    fn new(layouts: &[Option<Layout>], tables: &[Table], chain: &Chain) -> Result<Self, CaseError> {
        let fail = |m: String| Err(CaseError::Depletion(m));
        let mut files: BTreeMap<usize, &Micro> = BTreeMap::new();
        let uses = layouts
            .iter()
            .map(|l| {
                l.as_ref()
                    .map_or(Ok(Vec::new()), |l| node_uses(l, tables, &mut files))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let sites =
            |u: &[Option<NodeUse>]| -> Vec<bool> { u.iter().map(Option::is_some).collect() };
        let mut laid = uses.iter().filter(|u| !u.is_empty());
        if let Some(reference) = laid.next()
            && laid.any(|u| sites(u) != sites(reference))
        {
            return fail("the states micro deplete different sites".into());
        }
        let Some((&first, micro)) = files.iter().next() else {
            return fail("no fissile node to micro deplete".into());
        };
        let files_of = layouts
            .iter()
            .flatten()
            .next()
            .map(|l| l.files.clone())
            .unwrap_or_default();
        let name = |file: usize| {
            files_of
                .get(file)
                .map_or_else(|| format!("table {file}"), |p| p.display().to_string())
        };
        let names = &micro.nuclides;
        let mut built = None;
        for (&file, micro) in &files {
            if &micro.nuclides != names {
                return fail(format!(
                    "{} tracks the micro nuclides {:?}, but {} tracks {names:?}; all tables must track the same \
                     nuclides (ADR-0216)",
                    name(file),
                    micro.nuclides,
                    name(first)
                ));
            }
            let reference = micro.states.first().ok_or_else(|| {
                CaseError::Depletion(format!("{}: the micro data have no states", name(file)))
            })?;
            let c = MicroChain::new(chain, names, reference)
                .map_err(|e| CaseError::Depletion(format!("{}: {e}", name(file))))?;
            built.get_or_insert(c);
        }
        let chain = built.ok_or_else(|| CaseError::Depletion("no micro table".into()))?;
        Ok(Self { chain, uses })
    }

    /// 初めの組成: テーブルの燃焼度 0 の点の均質化した数密度（補助の核種は 0）。状態 `k` の並びのテーブルで引く
    fn initial(
        &self,
        k: usize,
        layouts: &[Option<Layout>],
        tables: &[Table],
    ) -> Result<Vec<Option<Vec<f64>>>, CaseError> {
        let layout = layouts[k]
            .as_ref()
            .ok_or_else(|| CaseError::Depletion("the first state is not laid out".into()))?;
        self.uses[k]
            .iter()
            .map(|u| {
                let Some(u) = u else { return Ok(None) };
                let table = match *u {
                    NodeUse::Uniform(table) => table,
                    // 区間は同じ層のテーブル（制御棒の値だけが違う）なので、どの区間も同じ燃焼度 0 の数密度
                    NodeUse::Partial(index) => layout.partial[index].segments[0].0,
                };
                let file = layout.tables[table].file;
                let table = &tables[file];
                let zero = table.axes.iter().any(|a| match a {
                    Axis::Burnup { points } => points.first() == Some(&0.0),
                    _ => false,
                });
                match (&table.micro, zero) {
                    (Some(micro), true) => Ok(Some(self.chain.initial(&micro.densities[0]))),
                    _ => Err(CaseError::Depletion(format!(
                        "{}: the micro depletion starts from the densities at burnup 0, which the table lacks",
                        layout.files[file].display()
                    ))),
                }
            })
            .collect()
    }
}

#[cfg(test)]
mod tests;
