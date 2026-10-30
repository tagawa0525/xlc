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

use std::collections::{BTreeMap, BTreeSet};

use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use xlc_depletion::CramOrder;
use xlc_depletion::chain::Chain;
use xlc_format::msgpack::{self, MsgpackError};
use xlc_format::table::{Axis, Micro, Table};

use crate::case::{
    CaseError, Compositions, Solved, check_search, search_boron, solve_compositions,
};
use crate::input::{
    Boron, BoronSearch, CaseState, DepletionInput, Layout, MidpointBoron, StepKind, Xenon,
};
use crate::micro::{MicroChain, TrackedRates};
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
    /// 探索の状態のステップの中点を解いたホウ素濃度（ADR-0191。ケースの `midpoint-boron`）
    pub midpoint_boron: MidpointBoron,
    /// ミクロ燃焼の Xe-135 と I-135 の扱い（ケースの `micro` の `xenon`）。ミクロ燃焼しなければ `None`
    #[serde(skip_serializing_if = "Option::is_none")]
    pub xenon: Option<Xenon>,
    /// ミクロ燃焼の補助の核種（追う核種へ崩壊で流れ込む核種）の名前（解く核種の順。ADR-0220）。ミクロ燃焼しなければ `None`
    #[serde(skip_serializing_if = "Option::is_none")]
    pub auxiliary: Option<Vec<String>>,
    /// 最後の点のノードの補助の核種の数密度 (1/(barn cm)。[集合体][軸方向のノード][補助の核種]。核分裂のないノードは空)。次の
    /// サイクルの restart が、最後の点の `densities` と合わせて読む（ADR-0220）。ミクロ燃焼しなければ `None`
    #[serde(skip_serializing_if = "Option::is_none")]
    pub auxiliary_densities: Option<Vec<Vec<Vec<f64>>>>,
    /// ミクロ燃焼のノードが持つ履歴の変数の名前（ADR-0215。`HISTORY_VARIABLES`）。ミクロ燃焼しなければ `None`
    #[serde(skip_serializing_if = "Option::is_none")]
    pub history_variables: Option<Vec<String>>,
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
    /// ミクロ燃焼のノードの履歴の変数の値（ADR-0215。[集合体][軸方向のノード][変数]。変数は `history_variables` の順で、核分裂のない
    /// ノードは空）。次のサイクルの restart が最後の点の値を読む。ミクロ燃焼しなければ `None`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub histories: Option<Vec<Vec<Vec<f64>>>>,
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
    /// rmpv の `from_value` は文字列を単位の列挙子として読まないので、名前で読む
    midpoint_boron: String,
    #[serde(default)]
    xenon: Option<String>,
    #[serde(default)]
    auxiliary: Option<Vec<String>>,
    #[serde(default)]
    auxiliary_densities: Option<Vec<Vec<Vec<f64>>>>,
    #[serde(default)]
    history_variables: Option<Vec<String>>,
    points: Vec<CorePoint>,
}

impl CoreHistory {
    /// フィールド名付きの map 形式の MessagePack: `{format: "xlc-core-history", specific_power, states: {<名前>: {boron | search:
    /// {boron: [下端, 上端], k, tolerance, max_solves}, b10, moderator, fuel, banks}}, axial_mesh, assemblies, heavy_metal,
    /// nuclides?, midpoint_boron: "search" | "start", xenon?: "tracked" | "equilibrium", auxiliary?, auxiliary_densities?, points: [{efpd,
    /// days, power?, state, burnup, k, boron, solves, assembly_burnups, burnups, powers, densities?}]}`
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
        let midpoint_boron = match file.midpoint_boron.as_str() {
            "search" => MidpointBoron::Search,
            "start" => MidpointBoron::Start,
            other => {
                return Err(CaseError::Depletion(format!(
                    "midpoint_boron must be search or start, found {other}"
                )));
            }
        };
        // Xe の扱いは、ミクロ燃焼の履歴だけが持つ
        let xenon = match (file.xenon.as_deref(), &file.nuclides) {
            (None, None) => None,
            (Some("tracked"), Some(_)) => Some(Xenon::Tracked),
            (Some("equilibrium"), Some(_)) => Some(Xenon::Equilibrium),
            (other, _) => {
                return Err(CaseError::Depletion(format!(
                    "xenon must be tracked or equilibrium in a micro depletion history and absent otherwise, found {other:?}"
                )));
            }
        };
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
        // 燃焼は燃焼点 0 を必ず書く
        if file.points.is_empty() {
            return Err(CaseError::Depletion("the history has no point".into()));
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
        // ミクロ燃焼の履歴は、最後の点の補助の核種の数密度を、追う核種の数密度を持つノード（核分裂のあるノード）に持つ（ADR-0220）
        let last = file.points.last().and_then(|p| p.densities.as_ref());
        let fits = match (
            &file.nuclides,
            &file.auxiliary,
            &file.auxiliary_densities,
            last,
        ) {
            (None, None, None, _) => true,
            (Some(_), Some(names), Some(rows), Some(tracked)) => {
                rows.len() == file.assemblies.len()
                    && rows.iter().all(|r| r.len() == layers)
                    && rows
                        .iter()
                        .flatten()
                        .zip(tracked.iter().flatten())
                        .all(|(a, t)| {
                            a.len() == if t.is_empty() { 0 } else { names.len() }
                                && a.iter().all(|v| v.is_finite())
                        })
            }
            _ => false,
        };
        if !fits {
            return Err(CaseError::Depletion(
                "the auxiliary densities must match the micro nuclides and the last point of the history (ADR-0220)"
                    .into(),
            ));
        }
        Ok(Self {
            specific_power: file.specific_power,
            states: file.states,
            axial_mesh: file.axial_mesh,
            assemblies: file.assemblies,
            heavy_metal: file.heavy_metal,
            nuclides: file.nuclides,
            midpoint_boron,
            xenon,
            auxiliary: file.auxiliary,
            auxiliary_densities: file.auxiliary_densities,
            history_variables: file.history_variables,
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

/// 履歴の変数の値 `h`（燃焼度 `start` での値）を、燃焼度 `end` まで、ステップの値 `x` で進める（ADR-0215 の燃焼度で重みを付けた平均）:
/// H(end) = (H(start) start + x (end − start)) / end。燃焼度 0 から進めると `x`（燃焼度 0 の値は瞬時の値）、燃焼度が進まなければ `h`
pub fn advance_history(h: &[f64], start: f64, end: f64, x: &[f64]) -> Vec<f64> {
    let _ = (start, end, x);
    h.to_vec()
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

/// ケースの燃焼（ADR-0190、ADR-0218）: 状態 `states`（ケースの名前つきの状態。ステップが番号で指す）で、ノードの燃焼度を
/// `depletion` のステップに沿って進め、燃焼点ごとの値を返す。`layout` は燃焼のケースから作った並び（`radial-nodes` = 1）。初めの
/// 燃焼度と組成は、サイクルが restart するなら前のサイクルの燃焼の履歴 `restart` の最後の点から集合体の ID で引き継ぎ（ADR-0220）、
/// ほかの集合体は 0 とテーブルの燃焼度 0 の組成。点は解くたびに `on_point` にも渡す（後のステップで止まっても、それまでの点を出せるように）
#[allow(
    clippy::too_many_arguments,
    reason = "入力の 4 つのファイルから作る並び、テーブル、状態、燃焼の入力、チェーン、前のサイクルの履歴は出どころが別で、まとめる型を作るより読みやすい"
)]
pub fn deplete(
    layout: &Layout,
    tables: &[Table],
    states: &[CaseState],
    depletion: &DepletionInput,
    chain: Option<&Chain>,
    restart: Option<&CoreHistory>,
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
    let xenon = depletion.micro.as_ref().map(|m| m.xenon);
    // 長い計算を始める前に
    if xenon == Some(Xenon::Equilibrium)
        && micro.as_ref().is_some_and(|m| !m.chain.balances_xenon())
    {
        return fail(
            "equilibrium xenon needs Xe-135 among the tracked nuclides and I-135 among the solved ones (ADR-0219)"
                .into(),
        );
    }
    let mut densities = micro
        .as_ref()
        .map(|m| m.initial(point_state(0), &layouts, tables))
        .transpose()?;
    // 初めの燃焼度。restart するサイクルは、引き継ぐ集合体の燃焼度と組成を前のサイクルの履歴から置く（ADR-0220）
    let mut burnups = vec![0.0; sites];
    match (layout.restart(), restart) {
        (None, None) => {}
        (Some(input), Some(history)) => {
            let ids: Vec<&str> = placed.keys().copied().collect();
            let inherited = inherit(
                history,
                &input.fresh,
                &ids,
                first.axial_mesh(),
                &heavy_metal,
                micro.as_ref().map(|m| &m.chain),
            )?;
            for (&id, state) in placed.keys().zip(inherited) {
                let Some((node_burnups, node_densities)) = state else {
                    continue;
                };
                let [x, y] = placed[id];
                for z in 0..nz {
                    let site = mesh.site([x, y, z]);
                    burnups[site] = node_burnups[z];
                    match (densities.as_mut().map(|n| &mut n[site]), &node_densities) {
                        (None, None) => {}
                        (Some(n), Some(rows)) if n.is_some() == rows[z].is_some() => {
                            if let Some(v) = &rows[z] {
                                *n = Some(v.clone());
                            }
                        }
                        _ => {
                            return fail(format!(
                                "assembly {id}: the restart history densities of node {z} do not match the micro \
                                 depleted nodes of the case (ADR-0220)"
                            ));
                        }
                    }
                }
            }
            let discharged: Vec<&str> = history
                .assemblies
                .iter()
                .map(String::as_str)
                .filter(|id| !placed.contains_key(id))
                .collect();
            tracing::info!(
                inherited = ids.len() - input.fresh.len(),
                fresh = input.fresh.len(),
                ?discharged,
                "restart"
            );
        }
        (Some(_), None) => {
            return fail("the cycle restarts from a history, but no history is given".into());
        }
        (None, Some(_)) => {
            return fail("a restart history is given, but the cycle does not restart".into());
        }
    }
    // 燃焼点 0 の炉心の平均の燃焼度。ステップの炉心の平均の燃焼度の目標は、これに定格の出力密度 × 累積の EFPD を足した値（ADR-0220）
    let start_burnup = weighted(&burnups, &masses);
    if let Some(m) = &micro {
        tracing::info!(
            tracked = m.chain.nuclides.len(),
            auxiliary = m.chain.auxiliary(),
            nodes = m.uses[point_state(0)].iter().flatten().count(),
            "micro depletion"
        );
    }
    // 状態 `k` の燃焼度 `burnups`（ミクロ燃焼なら組成 `n` も）で解く（探索の状態は、k が目標になるホウ素濃度で。`given` が
    // あれば、探索の状態もその濃度で解く）。`at` はエラーに添える解の場所
    let no_histories: Vec<Option<Vec<f64>>> = vec![None; sites];
    let solve = |k: usize,
                 burnups: &[f64],
                 n: Option<&[Option<Vec<f64>>]>,
                 given: Option<f64>,
                 at: &dyn Fn() -> String|
     -> Result<Solve, CaseError> {
        let state = &states[k];
        let layout = layouts[k]
            .as_ref()
            .ok_or_else(|| CaseError::Depletion(format!("state {} is not laid out", state.name)))?;
        let compositions = micro.as_ref().zip(n).map(|(m, densities)| Compositions {
            chain: &m.chain,
            densities,
            histories: &no_histories,
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
        let solved = match (&state.boron, given) {
            (Boron::Given(ppm), _) => solve_ppm(*ppm).map(|solved| (solved, *ppm, 1)),
            (Boron::Search(_), Some(ppm)) => solve_ppm(ppm).map(|solved| (solved, ppm, 1)),
            (Boron::Search(search), None) => search_boron(search, solve_ppm)
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
    // 解 `solve`（燃焼度 `burnups`）の位置 `site` の、組成 `n` で規格化した追う核種の反応率。ノードの出力密度は出力の割合 `power` ×
    // 定格の出力密度 × 解の核分裂率の割合 × 炉心の質量 / ノードの体積（燃焼度を進める `advance` と同じエネルギーを、割合で伸びた
    // 時間で出す。ADR-0218）。ミクロ燃焼しない位置は `None`
    let mass: f64 = masses.iter().sum();
    // `total` は解の核分裂率の炉心の和（位置ごとに足し直さない）
    let node_rates = |m: &MicroNodes,
                      solve: &Solve,
                      total: f64,
                      power: f64,
                      burnups: &[f64],
                      site: usize,
                      n: &[f64]|
     -> Result<Option<Vec<TrackedRates>>, String> {
        let state = states[solve.state].at(solve.ppm);
        let (Some(layout), uses) = (&layouts[solve.state], &m.uses[solve.state]) else {
            return Err(format!(
                "state {} is not laid out",
                states[solve.state].name
            ));
        };
        let Some(node_use) = uses[site] else {
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
                    format!("site {site}: the solve has no weights of the partially rodded node")
                })?;
                mix_micro(&micro, weights, &solve.solved.solution.flux[site])
                    .map_err(|e| format!("site {site}: {e}"))?
            }
        };
        m.chain
            .rates(&sigma, &solve.solved.solution.flux[site], n, power_density)
            .map(Some)
    };
    // 解 `solve`（燃焼度 `burnups`、組成 `n`）の反応率で、組成 `from` を時間 `t` (s) だけ進める
    let deplete_nodes = |m: &MicroNodes,
                         solve: &Solve,
                         power: f64,
                         burnups: &[f64],
                         n: &[Option<Vec<f64>>],
                         from: &[Option<Vec<f64>>],
                         t: f64,
                         at: &dyn Fn() -> String|
     -> Result<Vec<Option<Vec<f64>>>, CaseError> {
        let start = std::time::Instant::now();
        let total: f64 = solve.rates.iter().sum();
        let next = (0..sites)
            .into_par_iter()
            .map(|site| {
                let (Some(n), Some(from)) = (&n[site], &from[site]) else {
                    return Ok(None);
                };
                let Some(rates) = node_rates(m, solve, total, power, burnups, site, n)? else {
                    return Ok(None);
                };
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
    // 出力の割合 `power` で運転している解: 状態 `k`、燃焼度 `burnups`、組成 `n` で解く。Xe を平衡にするケース（`xenon` の
    // `equilibrium`。#0200）は、Xe-135、I-135 とその先行核の補助の核種を、解の反応率で釣り合う値に置き換えて解き直すことを、
    // 置き換える値が変わらなくなるまで繰り返す（Anderson の加速。`equilibrate`）。解と、解いた組成を返す
    let solve_operating = |k: usize,
                           burnups: &[f64],
                           n: Option<Densities>,
                           given: Option<f64>,
                           power: f64,
                           at: &dyn Fn() -> String|
     -> Result<(Solve, Option<Densities>), CaseError> {
        let (Some(m), Some(n)) = (
            micro.as_ref().filter(|_| xenon == Some(Xenon::Equilibrium)),
            n.as_ref(),
        ) else {
            return solve(k, burnups, n.as_deref(), given, at).map(|s| (s, n));
        };
        let n = n.clone();
        // 点の固有値計算の回数は、固定点のすべての解の和
        let mut solves = 0;
        equilibrate(n, |n| {
            let solved = solve(k, burnups, Some(n), given, at)?;
            solves += solved.solves;
            let total: f64 = solved.rates.iter().sum();
            let balanced = (0..sites)
                .into_par_iter()
                .map(|site| {
                    let Some(n) = &n[site] else { return Ok(None) };
                    match node_rates(m, &solved, total, power, burnups, site, n)? {
                        Some(rates) => m.chain.equilibrium(&rates, n).map(Some),
                        None => Ok(None),
                    }
                })
                .collect::<Result<Vec<_>, String>>()
                .map_err(|e| locate(CaseError::Depletion(e), at))?;
            Ok((solved, balanced))
        })
        .map(|(s, n)| (Solve { solves, ..s }, Some(n)))
        .map_err(|e| locate(e, at))
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
            histories: None,
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

    let mut efpd = 0.0;
    let mut days = 0.0;
    let mut solve_at_point = solve(
        point_state(0),
        &burnups,
        densities.as_deref(),
        None,
        &|| "the point at 0 EFPD".into(),
    )?;
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
                // 終わるケースがその点を越えて止まる）。出力の割合によらない（EFPD は割合を含む。ADR-0218）。restart したサイクルは
                // 燃焼点 0 の炉心の平均から数える（ADR-0220）
                let target = start_burnup + depletion.specific_power * end * 1.0e-3;
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
                    // 案 S1 の中点は、ステップの初めの点（同じ状態で解いた点）のホウ素濃度で解く
                    let given = match depletion.midpoint_boron {
                        MidpointBoron::Search => None,
                        MidpointBoron::Start => Some(solve_at_point.ppm),
                    };
                    let (middle, half) =
                        solve_operating(step.state, b, half, given, fraction, &|| {
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
        // 出力のステップの終わりの点は運転している解。点は次のステップの状態で解く（ADR-0218 の案 P2）ので、次も出力のステップなら
        // その出力で、次が崩壊だけのステップか最後の点なら終わったステップの出力で運転している。崩壊だけのステップの後の点は、崩壊で
        // 進めた組成のまま解く（停止の後の Xe は、平衡でなく崩壊の結果）
        let next = depletion.steps.get(i + 1).map(|s| s.kind);
        (solve_at_point, densities) = match step.kind {
            // 崩壊だけのステップの前の点の Xe は、停止の前に運転していた状態（終わったステップの状態）の中性子束で釣り合う。点は
            // 崩壊だけのステップの状態で、その組成で解く
            StepKind::Power { power, .. }
                if xenon == Some(Xenon::Equilibrium)
                    && matches!(next, Some(StepKind::Decay { .. }))
                    && point_state(i + 1) != step.state =>
            {
                let (operating, balanced) =
                    solve_operating(step.state, &burnups, densities.take(), None, power, &|| {
                        at_point.clone()
                    })?;
                let point = solve(
                    point_state(i + 1),
                    &burnups,
                    balanced.as_deref(),
                    None,
                    &|| at_point.clone(),
                )?;
                (
                    Solve {
                        solves: operating.solves + point.solves,
                        ..point
                    },
                    balanced,
                )
            }
            StepKind::Power { power, .. } => solve_operating(
                point_state(i + 1),
                &burnups,
                densities.take(),
                None,
                match next {
                    Some(StepKind::Power { power: next, .. }) => next,
                    _ => power,
                },
                &|| at_point.clone(),
            )?,
            StepKind::Decay { .. } => (
                solve(
                    point_state(i + 1),
                    &burnups,
                    densities.as_deref(),
                    None,
                    &|| at_point.clone(),
                )?,
                densities,
            ),
        };
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
    // 最後の点の補助の核種（次のサイクルの restart が読む。ADR-0220）
    let auxiliary_densities = micro.as_ref().zip(densities.as_ref()).map(|(m, n)| {
        placed
            .values()
            .map(|&[x, y]| {
                (0..nz)
                    .map(|z| {
                        n[mesh.site([x, y, z])]
                            .as_ref()
                            .map_or_else(Vec::new, |v| m.chain.auxiliary_densities(v))
                    })
                    .collect()
            })
            .collect()
    });
    Ok(CoreHistory {
        specific_power: depletion.specific_power,
        states: used
            .iter()
            .map(|&k| (states[k].name.clone(), core_state(&states[k])))
            .collect(),
        axial_mesh: first.axial_mesh().to_vec(),
        assemblies: placed.keys().map(|id| id.to_string()).collect(),
        heavy_metal,
        auxiliary: micro.as_ref().map(|m| m.chain.auxiliary_nuclides()),
        auxiliary_densities,
        history_variables: None,
        nuclides: micro.map(|m| m.chain.nuclides),
        midpoint_boron: depletion.midpoint_boron,
        xenon,
        points,
    })
}

/// 1 日 (s)
const DAY: f64 = 86_400.0;

/// restart で引き継ぐ 1 体の集合体の状態: 軸方向のノードの燃焼度と、ミクロ燃焼なら軸方向のノードの解く核種の数密度（核分裂のない
/// ノードは `None`）
type Inherited = (Vec<f64>, Option<Vec<Option<Vec<f64>>>>);

/// 前のサイクルの燃焼の履歴 `history` の最後の点から、装荷した集合体 `ids`（名前の順）の状態を引き継ぐ（ADR-0220）。`fresh` の
/// 集合体は新しい燃料で `None`。`axial_mesh` はケースの軸方向のメッシュ、`heavy_metal` は `ids` の順のノードの初めの重金属の質量、
/// `chain` はミクロ燃焼の解く核種（ミクロ燃焼しなければ `None`）。履歴と食い違えば推測せずにエラー
fn inherit(
    history: &CoreHistory,
    fresh: &BTreeSet<String>,
    ids: &[&str],
    axial_mesh: &[f64],
    heavy_metal: &[Vec<f64>],
    chain: Option<&MicroChain>,
) -> Result<Vec<Option<Inherited>>, CaseError> {
    fn fail<T>(m: String) -> Result<T, CaseError> {
        Err(CaseError::Depletion(m))
    }
    let Some(last) = history.points.last() else {
        return fail("the restart history has no point".into());
    };
    if history.axial_mesh != axial_mesh {
        return fail(format!(
            "the axial mesh {axial_mesh:?} differs from the axial mesh {:?} of the restart history",
            history.axial_mesh
        ));
    }
    // ミクロ燃焼の有無と、追う核種と補助の核種の並び
    let micro = match (chain, &history.nuclides, &history.auxiliary) {
        (None, None, None) => None,
        (Some(c), Some(tracked), Some(auxiliary)) => {
            if tracked != &c.nuclides {
                return fail(format!(
                    "the restart history tracks the nuclides {tracked:?}, but the tables track {:?}",
                    c.nuclides
                ));
            }
            if auxiliary != &c.auxiliary_nuclides() {
                return fail(format!(
                    "the restart history has the auxiliary nuclides {auxiliary:?}, but the chain solves {:?}",
                    c.auxiliary_nuclides()
                ));
            }
            match (&last.densities, &history.auxiliary_densities) {
                (Some(t), Some(a)) => Some((c, t, a)),
                _ => {
                    return fail(
                        "the restart history has no densities of the tracked and auxiliary nuclides at its last point"
                            .into(),
                    );
                }
            }
        }
        (None, Some(_), _) => {
            return fail(
                "the restart history is a micro depletion, but the case does not micro deplete"
                    .into(),
            );
        }
        _ => {
            return fail(
                "the case micro depletes, but the restart history is not a micro depletion".into(),
            );
        }
    };
    let nz = axial_mesh.len().saturating_sub(1);
    ids.iter()
        .zip(heavy_metal)
        .map(|(&id, masses)| {
            let at = history.assemblies.iter().position(|a| a == id);
            let a = match (at, fresh.contains(id)) {
                (None, true) => return Ok(None),
                (Some(a), false) => a,
                (Some(_), true) => {
                    return fail(format!(
                        "assembly {id} is in the restart history but named fresh"
                    ));
                }
                (None, false) => {
                    return fail(format!(
                        "assembly {id} is neither in the restart history nor named fresh"
                    ));
                }
            };
            if history.heavy_metal.get(a) != Some(masses) {
                return fail(format!(
                    "assembly {id}: the heavy metal {masses:?} g of its nodes differs from {:?} g of the restart history",
                    history.heavy_metal.get(a)
                ));
            }
            let burnups = match last.burnups.get(a) {
                Some(b) if b.len() == nz => b.clone(),
                _ => {
                    return fail(format!(
                        "assembly {id}: the restart history has no burnups of its {nz} nodes"
                    ));
                }
            };
            let densities = micro
                .map(|(c, tracked, auxiliary)| {
                    let (Some(t), Some(x)) = (tracked.get(a), auxiliary.get(a)) else {
                        return fail(format!("assembly {id}: the restart history has no densities"));
                    };
                    t.iter()
                        .zip(x)
                        .map(|(t, x)| {
                            if t.is_empty() && x.is_empty() {
                                Ok(None)
                            } else {
                                c.compose(t, x).map(Some).map_err(|e| {
                                    CaseError::Depletion(format!(
                                        "assembly {id}: the restart history densities: {e}"
                                    ))
                                })
                            }
                        })
                        .collect::<Result<Vec<_>, _>>()
                })
                .transpose()?;
            if densities.as_ref().is_some_and(|d| d.len() != nz) {
                return fail(format!(
                    "assembly {id}: the restart history densities do not have its {nz} nodes"
                ));
            }
            Ok(Some((burnups, densities)))
        })
        .collect()
}

/// 位置ごとの解く核種の数密度（ミクロ燃焼しない位置は `None`）
type Densities = Vec<Option<Vec<f64>>>;

/// 固定点の反復の組: 並べた組成 x と、その解で平衡にした組成 g
type Iterate = (Vec<f64>, Vec<f64>);

/// 平衡 Xe の固定点の許容: 置き換える数密度の、解き直す前後の相対の変化（初めの平衡の値に対する比）の最大
const EQUILIBRIUM_TOLERANCE: f64 = 1.0e-6;
/// 平衡 Xe の固定点の、解く回数の上限
const EQUILIBRIUM_SOLVES: usize = 50;
/// Anderson の加速で使う、前の反復の数
const ANDERSON_DEPTH: usize = 3;

/// 組成 `n` から、`evaluate`（組成 → 解と、その解の反応率で平衡にした組成）の固定点を Anderson の加速（type II）で求め、固定点の
/// 組成で解いた解と、その組成を返す。中性子束が Xe に強く応じる炉心では、平衡の組成をそのまま次の組成にする反復は互い違いに振れて
/// 収束しないことがあるので、前の反復の差の最小二乗で次の組成を選ぶ
fn equilibrate<S>(
    n: Densities,
    mut evaluate: impl FnMut(&[Option<Vec<f64>>]) -> Result<(S, Densities), CaseError>,
) -> Result<(S, Densities), CaseError> {
    let flatten =
        |v: &[Option<Vec<f64>>]| -> Vec<f64> { v.iter().flatten().flatten().copied().collect() };
    let shape: Vec<Option<usize>> = n.iter().map(|v| v.as_ref().map(Vec::len)).collect();
    let unflatten = |x: &[f64]| -> Vec<Option<Vec<f64>>> {
        let mut at = 0;
        shape
            .iter()
            .map(|len| {
                len.map(|len| {
                    at += len;
                    x[at - len..at].to_vec()
                })
            })
            .collect()
    };
    let mut x = n;
    let mut history: Vec<Iterate> = Vec::new();
    let mut scale: Option<Vec<f64>> = None;
    for iteration in 1..=EQUILIBRIUM_SOLVES {
        let (solved, balanced) = evaluate(&x)?;
        let (xf, gf) = (flatten(&x), flatten(&balanced));
        if balanced
            .iter()
            .map(|v| v.as_ref().map(Vec::len))
            .ne(shape.iter().copied())
        {
            return Err(CaseError::Depletion(
                "the equilibrium changed the micro depleted sites".into(),
            ));
        }
        let s = scale.get_or_insert_with(|| {
            gf.iter()
                .map(|v| if *v != 0.0 { v.abs() } else { 1.0 })
                .collect()
        });
        let residual = gf
            .iter()
            .zip(&xf)
            .zip(s.iter())
            .map(|((g, x), s)| (g - x).abs() / s)
            .fold(0.0, f64::max);
        tracing::debug!(iteration, residual, "equilibrium xenon");
        if residual <= EQUILIBRIUM_TOLERANCE {
            tracing::info!(iterations = iteration, residual, "equilibrium xenon");
            return Ok((solved, x));
        }
        let last = (xf, gf);
        x = unflatten(&anderson(&history, &last, s));
        history.push(last);
        if history.len() > ANDERSON_DEPTH {
            history.remove(0);
        }
    }
    Err(CaseError::Depletion(format!(
        "the equilibrium xenon did not converge in {EQUILIBRIUM_SOLVES} solves"
    )))
}

/// Anderson の加速（type II）の次の点: 前の反復の組（x、g = G(x)）`previous`（古い順）と最後の組 `last` から、残差
/// f = (g − x) / `scale` の差の最小二乗で γ を求め、g − ΔG γ を返す（負の数密度は 0 にする）。前の反復がないか最小二乗が特異なら、
/// 最後の g
fn anderson(previous: &[Iterate], last: &Iterate, scale: &[f64]) -> Vec<f64> {
    let all: Vec<&Iterate> = previous.iter().chain([last]).collect();
    let residual = |(x, g): &Iterate| -> Vec<f64> {
        g.iter()
            .zip(x)
            .zip(scale)
            .map(|((g, x), s)| (g - x) / s)
            .collect()
    };
    let f: Vec<Vec<f64>> = all.iter().map(|i| residual(i)).collect();
    let m = previous.len();
    if m == 0 {
        return last.1.clone();
    }
    let df: Vec<Vec<f64>> = (0..m)
        .map(|j| f[j + 1].iter().zip(&f[j]).map(|(a, b)| a - b).collect())
        .collect();
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
    let normal: Vec<Vec<f64>> = (0..m)
        .map(|i| (0..m).map(|j| dot(&df[i], &df[j])).collect())
        .collect();
    let right: Vec<f64> = (0..m).map(|i| dot(&df[i], &f[m])).collect();
    let Some(gamma) = xlc_numeric::lu::solve(normal.concat(), right) else {
        return last.1.clone();
    };
    let mut next = last.1.clone();
    for (j, g) in gamma.iter().enumerate() {
        let (a, b) = (&all[j].1, &all[j + 1].1);
        for ((v, b), a) in next.iter_mut().zip(b).zip(a) {
            *v -= g * (b - a);
        }
    }
    next.iter().map(|v| v.max(0.0)).collect()
}

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
