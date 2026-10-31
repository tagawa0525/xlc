//! ケースの状態ごとの固有値計算（ADR-0146）と、臨界ホウ素濃度の探索（ADR-0145）。

use rayon::prelude::*;
use xlc_format::fewgroup::FewGroup;
use xlc_format::table::{Axis, HistoryRod, MODERATOR_DENSITY, ROD_OUT, Table};

use crate::input::{Boron, BoronSearch, CaseState, InputError, Layout, TableUse};
use crate::micro::{MicroChain, correct, production_index};
use crate::nem::{NemError, Options, Problem, Solution, solve};
use crate::rods::{SegmentWeights, solve_rodded};
use crate::search::{SearchError, find};
use crate::state::{
    MicroAt, State, add_rate_residual, add_residual, interpolate, interpolate_branch_difference,
    interpolate_micro, interpolate_rate_residual, interpolate_residual, non_negative,
};

/// 状態を解いた結果。
#[derive(Debug, Clone)]
pub struct Solved {
    /// 探索したホウ素濃度（探索した状態だけ）
    pub searched: Option<Searched>,
    pub problem: Problem,
    pub solution: Solution,
    /// 部分挿入のノード（`Layout::partial` の順）の定数を混ぜた [区間][群] の重み（ADR-0127、ADR-0217）
    pub segment_weights: SegmentWeights,
}

/// 探索したホウ素濃度と、そのための固有値計算の回数。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Searched {
    pub ppm: f64,
    pub solves: usize,
}

/// 状態を解くときのエラー。
#[derive(Debug, thiserror::Error)]
pub enum CaseError {
    /// テーブルの内挿のエラー（ファイルの名前を含む）
    #[error("{0}")]
    State(String),
    #[error(transparent)]
    Input(#[from] InputError),
    #[error(transparent)]
    Nem(#[from] NemError),
    #[error("boron search: {0}")]
    Search(String),
    /// 燃焼のステップのエラー（ADR-0190）
    #[error("depletion: {0}")]
    Depletion(String),
    /// 燃焼のステップのどの解（点か中点）か、どのステップで起きたか（ADR-0191）
    /// 表示が原因を含むので、原因をエラーの source にしない（連鎖の表示で繰り返さない）
    #[error("{at}: {cause}")]
    At { at: String, cause: Box<CaseError> },
}

/// テーブル（`layout.files` の順）を状態 `state`（`None` なら軸のないテーブル）で、テーブルの使い方（`layout.tables`。
/// 制御棒の区間の値）ごとに内挿し、炉心を解く。部分挿入のノードは `layout.partial_rods` の方法で混ぜた定数で解き、
/// 返す問題の物質はその定数にする（核分裂率とピン出力に使う。ADR-0127）
pub fn solve_at(
    layout: &Layout,
    tables: &[Table],
    state: Option<&State>,
    options: &Options,
) -> Result<Solved, CaseError> {
    solve_with(layout, tables, |_| state.copied(), None, options)
}

/// core がノードごとに持つ履歴の変数（ADR-0215）。ノードの値の並びはこの順
pub const HISTORY_VARIABLES: [&str; 1] = [MODERATOR_DENSITY];

/// ミクロ燃焼のノードの組成（ADR-0216）と履歴の変数（ADR-0215）
pub struct Compositions<'a> {
    pub chain: &'a MicroChain,
    /// [直方体の位置] の解く核種の数密度（`MicroChain::solved` の順。ミクロ燃焼しない位置は `None`）
    pub densities: &'a [Option<Vec<f64>>],
    /// [直方体の位置] の履歴の変数の値（`HISTORY_VARIABLES` の順。ミクロ燃焼しない位置は `None`）
    pub histories: &'a [Option<Vec<f64>>],
}

/// テーブル `table` の履歴の補正のうち、ノードの履歴の変数の値 `histories`（`HISTORY_VARIABLES` の順）の重み
/// w_j = (H_j − H_j⁰) / (H_j^b − H_j⁰) が 0 でない変数の、重みと制御棒の値 `rod`（棒のない区間は `out`）の補正（変数の名前、w_j、補正）。
/// テーブルが core の持たない変数の補正を持つか、w_j が 0 でないのにその制御棒の値の補正を持たなければエラー
fn weighted_histories<'a>(
    table: &'a Table,
    rod: Option<&str>,
    histories: &[f64],
) -> Result<Vec<(&'a str, f64, &'a HistoryRod)>, String> {
    if histories.len() != HISTORY_VARIABLES.len() {
        return Err(format!(
            "{} history variable values for the history variables {HISTORY_VARIABLES:?}",
            histories.len()
        ));
    }
    let mut out = Vec::new();
    for correction in &table.histories {
        let variable = &correction.variable;
        let k = HISTORY_VARIABLES
            .iter()
            .position(|v| v == variable)
            .ok_or_else(|| {
                format!(
                    "the table holds the residuals of the history variable {variable:?}, which core does not \
                     track ({HISTORY_VARIABLES:?}. ADR-0215)"
                )
            })?;
        let w = (histories[k] - correction.nominal) / (correction.branch - correction.nominal);
        // 公称の履歴の変数は補正を使わない（補正を持たないテーブルとビットで同じにする。0 倍を足すと −0 が +0 になる）
        if w == 0.0 {
            continue;
        }
        if !w.is_finite() {
            return Err(format!(
                "history variable {variable}: the value {} gives the weight {w}",
                histories[k]
            ));
        }
        let rod = rod.unwrap_or(ROD_OUT);
        let entry = correction
            .rods
            .iter()
            .find(|r| r.rod == rod)
            .ok_or_else(|| {
                format!(
                    "history variable {variable}: the node has the rod value {rod} and the weight {w}, but the table \
                     holds no residuals for that rod value (ADR-0215)"
                )
            })?;
        out.push((variable.as_str(), w, entry));
    }
    Ok(out)
}

/// ミクロ燃焼のノードの少数群定数 `f`（テーブルを状態で内挿し、線形の補正を足したもの）に、テーブルの履歴の残り
/// Σ_j [w_j U_j,r(u) + (E / E_j,r(u)) L_j,r(u)] を足す（ADR-0215）。w_j と制御棒の値の選び方は `weighted_histories`、U と L は燃焼度
/// `burnup` で内挿した値。E は熱群の νΣf の線形の補正（`production_index`）の、ノードの数密度 `n` と状態で内挿した基準の数密度
/// （`micro`）の差のもの、E_j は同じ係数の、履歴の分岐の数密度の差（`interpolate_branch_difference`）のもの。w_j が 0 の変数は足さない。
/// 履歴の分岐の数密度の差がすべて 0 の燃焼度（燃焼度 0）では、ハイブリッドは公称の組成そのもので L は定義から 0 なので、L を足さない
/// （テーブルの値は 2 つの格子計算の解の差）。燃焼度が残りの点の範囲の外か、数密度の差が 0 でないのに E_j が 0 ならエラー
pub fn add_histories(
    f: &mut FewGroup,
    table: &Table,
    burnup: f64,
    rod: Option<&str>,
    histories: &[f64],
    micro: &MicroAt,
    n: &[f64],
) -> Result<(), String> {
    let weighted = weighted_histories(table, rod, histories)?;
    if weighted.is_empty() {
        return Ok(());
    }
    let Some(Axis::Burnup { points: axis }) =
        table.axes.iter().find(|a| matches!(a, Axis::Burnup { .. }))
    else {
        return Err("history residuals need a burnup axis (ADR-0215)".into());
    };
    let nominal = &table
        .micro
        .as_ref()
        .ok_or("history residuals need micro data (ADR-0215)")?
        .densities;
    let node: Vec<f64> = n
        .iter()
        .zip(&micro.densities)
        .map(|(density, nominal)| density - nominal)
        .collect();
    for (variable, w, entry) in weighted {
        let at = |e: String| format!("history variable {variable}: {e}");
        let untracked = interpolate_residual(entry, &entry.untracked, burnup).map_err(at)?;
        add_residual(f, w, &untracked);
        let branch = interpolate_branch_difference(entry, axis, nominal, burnup).map_err(at)?;
        if branch.iter().all(|&d| d == 0.0) {
            continue;
        }
        let linearization =
            interpolate_residual(entry, &entry.linearization, burnup).map_err(at)?;
        let own = production_index(&micro.nuclides, &node);
        let of_branch = production_index(&micro.nuclides, &branch);
        let ratio = own / of_branch;
        // 比を作れない点で L を w などで足すと、どの近似かが入力から読めない。黙って置き換えない
        if of_branch == 0.0 || !ratio.is_finite() {
            return Err(at(format!(
                "the branch densities differ from the nominal ones at burnup {burnup}, but the branch index (the \
                 thermal nu-fission of their linear correction) is {of_branch}, so its ratio {own} / {of_branch} is \
                 not a weight of the linearization error (ADR-0215)"
            )));
        }
        add_residual(f, ratio, &linearization);
    }
    Ok(())
}

/// `solve_at` と同じく解くが、燃焼のケースのノード（テーブルの使い方が `node` を持つもの）は、直方体の位置ごとの燃焼度
/// `burnups` で引く（ADR-0190）。ほかのテーブルの使い方（反射体の位置）は燃焼しないので、燃焼度なしで引く（燃焼度の軸を持つ
/// テーブルはエラー）
pub fn solve_burnups(
    layout: &Layout,
    tables: &[Table],
    state: &State,
    burnups: &[f64],
    options: &Options,
) -> Result<Solved, CaseError> {
    solve_compositions(layout, tables, state, burnups, None, options)
}

/// `solve_burnups` と同じく解くが、組成 `compositions` を持つ位置の定数に、ミクロ燃焼の線形の補正 Σᵢ σᵢ (Nᵢ − N̄ᵢ⁰) を足す
/// （ADR-0216。σᵢ はその位置のテーブルの使い方の状態と燃焼度で、N̄ᵢ⁰ は燃焼度で内挿した値）。テーブルが履歴の残りを持てば、
/// 位置の履歴の変数で重みを付けた残りも足す（ADR-0215。`node_constants`）。部分挿入のノードは、区間のテーブルの使い方ごとに
/// 足した定数を混ぜるので、残りも少数群定数と同じ重みで混ざる（ADR-0127）
pub fn solve_compositions(
    layout: &Layout,
    tables: &[Table],
    state: &State,
    burnups: &[f64],
    compositions: Option<&Compositions>,
    options: &Options,
) -> Result<Solved, CaseError> {
    if let Some(b) = state.burnup {
        return Err(CaseError::Depletion(format!(
            "the state gives a burnup {b}; the burnups are given by node"
        )));
    }
    if burnups.len() != layout.node_material.len() {
        return Err(CaseError::Depletion(format!(
            "{} burnups for {} positions",
            burnups.len(),
            layout.node_material.len()
        )));
    }
    solve_with(
        layout,
        tables,
        |t| {
            Some(State {
                burnup: t.node.map(|n| burnups[n]),
                ..*state
            })
        },
        compositions,
        options,
    )
}

/// ミクロ燃焼のノードの少数群定数: テーブル `table` を状態 `state`（ノードの燃焼度を持つ）と制御棒の値 `rod` で内挿した定数に、追う
/// 核種の数密度 `n`（テーブルの順）の線形の補正（ADR-0216）と、履歴の変数の値 `histories`（`HISTORY_VARIABLES` の順）の残り（ADR-0215。
/// `add_histories`）を足す。ミクロ燃焼のテーブルは分けない（格子全体の均質化）
pub fn node_constants(
    table: &Table,
    state: &State,
    rod: Option<&str>,
    n: &[f64],
    histories: &[f64],
) -> Result<FewGroup, String> {
    let mut constants = interpolate(table, Some(state), rod)?;
    let micro = interpolate_micro(table, state, rod)?;
    let (Some(mut f), true) = (constants.pop(), constants.is_empty()) else {
        return Err("micro data need an unsplit table".into());
    };
    correct(&mut f, &micro, n);
    let burnup = state
        .burnup
        .ok_or("the micro depleted node needs a burnup (ADR-0190)")?;
    add_histories(&mut f, table, burnup, rod, histories, &micro, n)?;
    Ok(f)
}

/// ミクロ燃焼のノードの反応率の微視的断面積: テーブル `table` を状態 `state` と制御棒の値 `rod` で内挿した値に、履歴の変数の値
/// `histories` の差を足したもの（ADR-0215 の案 S2）
pub fn node_micro(
    table: &Table,
    state: &State,
    rod: Option<&str>,
    histories: &[f64],
) -> Result<MicroAt, String> {
    let mut micro = interpolate_micro(table, state, rod)?;
    let weighted = weighted_histories(table, rod, histories)?;
    if weighted.is_empty() {
        return Ok(micro);
    }
    let burnup = state
        .burnup
        .ok_or("the micro depleted node needs a burnup (ADR-0190)")?;
    let names = &table
        .micro
        .as_ref()
        .ok_or("the table has no micro data for the micro depletion (ADR-0216)")?
        .nuclides;
    for (variable, w, entry) in weighted {
        let r = interpolate_rate_residual(entry, burnup)
            .map_err(|e| format!("history variable {variable}: {e}"))?;
        add_rate_residual(&mut micro, names, w, &r)
            .map_err(|e| format!("history variable {variable}: {e}"))?;
    }
    Ok(micro)
}

/// テーブルの使い方ごとに、`state_of` の状態でテーブルを内挿して解く
fn solve_with(
    layout: &Layout,
    tables: &[Table],
    state_of: impl Fn(&TableUse) -> Option<State>,
    compositions: Option<&Compositions>,
    options: &Options,
) -> Result<Solved, CaseError> {
    if let Some(c) = compositions
        && (c.densities.len() != layout.node_material.len()
            || c.histories.len() != layout.node_material.len())
    {
        return Err(CaseError::Depletion(format!(
            "{} micro compositions and {} node histories for {} positions",
            c.densities.len(),
            c.histories.len(),
            layout.node_material.len()
        )));
    }
    if tables.len() != layout.files.len() {
        return Err(CaseError::State(format!(
            "{} tables for {} few-group constant files",
            tables.len(),
            layout.files.len()
        )));
    }
    let constants = layout
        .tables
        .iter()
        .map(|t| {
            let file = layout.files[t.file].display();
            let table = &tables[t.file];
            if t.cut && table.axes.iter().any(|a| matches!(a, Axis::Rod { .. })) {
                return Err(CaseError::State(format!(
                    "{file}: the layer above the rod top has a rod axis; rods are cut at the top \
                     only where the layer tables have none (ADR-0128)"
                )));
            }
            let state = state_of(t);
            // ミクロ燃焼のノード: 組成と履歴の変数を持つ位置（ADR-0216、ADR-0215）
            let node = compositions.zip(t.node).and_then(|(c, site)| {
                c.densities
                    .get(site)
                    .and_then(Option::as_ref)
                    .map(|n| (c, site, n))
            });
            let Some((c, site, n)) = node else {
                // 履歴の残りは、履歴の変数を進めるミクロ燃焼のノードにだけ足せる。黙って補正なしで使わない
                if !table.histories.is_empty() {
                    return Err(CaseError::State(format!(
                        "{file}: the table holds history residuals (ADR-0215), which core adds only to micro depleted \
                         nodes, whose history variables it advances"
                    )));
                }
                return interpolate(table, state.as_ref(), t.rod.as_deref())
                    .map_err(|e| CaseError::State(format!("{file}: {e}")));
            };
            let s = state.as_ref().ok_or_else(|| {
                CaseError::State(format!("{file}: micro depletion needs a state"))
            })?;
            let histories = c
                .histories
                .get(site)
                .and_then(Option::as_ref)
                .ok_or_else(|| {
                    CaseError::State(format!(
                        "{file}: the micro depleted position {site} has no history variables"
                    ))
                })?;
            let f = node_constants(table, s, t.rod.as_deref(), &c.chain.tracked(n), histories)
                .map_err(|e| CaseError::State(format!("{file}: {e}")))?;
            // 補正で数密度の差が大きいと、吸収や全断面積が負になりうる。内挿の後と同じ検査をして、ソルバに渡す前に止める
            f.validate()
                .map_err(|e| e.to_string())
                .and_then(|()| non_negative(&f))
                .map_err(|e| {
                    CaseError::State(format!(
                        "{file}: the micro and history corrections at position {site} gave invalid constants: {e}"
                    ))
                })?;
            Ok(vec![f])
        })
        .collect::<Result<Vec<_>, _>>()?;
    let problem = layout.problem(&constants)?;
    let partial = layout.partial_nodes(&constants)?;
    if partial.is_empty() {
        let solution = solve(&problem, options)?;
        return Ok(Solved {
            searched: None,
            problem,
            solution,
            segment_weights: Vec::new(),
        });
    }
    let (solution, materials, segment_weights) = solve_rodded(
        &problem,
        partial,
        layout.partial_rods.weighting,
        layout.partial_rods.fine_mesh,
        options,
    )?;
    Ok(Solved {
        searched: None,
        problem: Problem {
            materials,
            ..problem
        },
        solution,
        segment_weights,
    })
}

/// 探索の状態の、解く前にできる検査。
pub fn check_search(state: &CaseState, options: &Options) -> Result<(), CaseError> {
    if let Boron::Search(search) = &state.boron
        // 固有値の反復の判定（見積もった残りの誤差）より細かい許容は区別できない（判定を超えても保証ではない。ADR-0145）
        && search.tolerance <= options.k_tolerance
    {
        return Err(CaseError::Search(format!(
            "tolerance {} must exceed the eigenvalue tolerance {}",
            search.tolerance, options.k_tolerance
        )));
    }
    Ok(())
}

/// ケースの状態を1つ解く。探索の状態は、k が目標になるホウ素濃度を探す（ADR-0145）。
pub fn solve_state(
    layout: &Layout,
    tables: &[Table],
    state: &CaseState,
    options: &Options,
) -> Result<Solved, CaseError> {
    // 状態のバンクの位置で制御棒を入れた並び（ADR-0127）
    let rodded = layout.with_banks(&state.banks)?;
    let layout = &rodded;
    let search = match &state.boron {
        Boron::Given(ppm) => return solve_at(layout, tables, Some(&state.at(*ppm)), options),
        Boron::Search(search) => search,
    };
    check_search(state, options)?;
    let (solved, searched) = search_boron(search, |ppm| {
        solve_at(layout, tables, Some(&state.at(ppm)), options)
    })?;
    Ok(Solved {
        searched: Some(searched),
        ..solved
    })
}

/// `search` の範囲で、`solve`（ホウ素濃度 → 解）の k が目標になるホウ素濃度を Illinois 法で探し（ADR-0145）、最後に解いた
/// 解と、求めたホウ素濃度と回数を返す。検査（`check_search`）は呼ぶ側が済ませる
pub fn search_boron(
    search: &BoronSearch,
    mut solve: impl FnMut(f64) -> Result<Solved, CaseError>,
) -> Result<(Solved, Searched), CaseError> {
    let found = find(
        search.bracket,
        search.k,
        search.tolerance,
        search.max_solves,
        |ppm| {
            let solved = solve(ppm)?;
            tracing::debug!(ppm, k = solved.solution.k, "boron search");
            Ok((solved.solution.k, solved))
        },
    )
    .map_err(|e| match e {
        SearchError::Solve(e) => e,
        e => CaseError::Search(e.to_string()),
    })?;
    let searched = Searched {
        ppm: found.x,
        solves: found.solves,
    };
    Ok((found.at, searched))
}

/// ケースの状態を並列に解く（結果は `states` の順）。各状態の解は1つのスレッドで解くので、スレッドの数によらない。
pub fn solve_states(
    layout: &Layout,
    tables: &[Table],
    states: &[CaseState],
    options: &Options,
) -> Vec<Result<Solved, CaseError>> {
    states
        .par_iter()
        .map(|state| {
            let _span = tracing::info_span!("state", name = %state.name).entered();
            solve_state(layout, tables, state, options)
        })
        .collect()
}

#[cfg(test)]
mod tests;
