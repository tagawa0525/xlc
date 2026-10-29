//! ケースの状態ごとの固有値計算（ADR-0146）と、臨界ホウ素濃度の探索（ADR-0145）。

use rayon::prelude::*;
use xlc_format::table::{Axis, Table};

use crate::input::{Boron, BoronSearch, CaseState, InputError, Layout, TableUse};
use crate::micro::{MicroChain, correct};
use crate::nem::{NemError, Options, Problem, Solution, solve};
use crate::rods::{SegmentWeights, solve_rodded};
use crate::search::{SearchError, find};
use crate::state::{State, interpolate, interpolate_micro, non_negative};

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

/// ミクロ燃焼のノードの組成（ADR-0216）
pub struct Compositions<'a> {
    pub chain: &'a MicroChain,
    /// [直方体の位置] の解く核種の数密度（`MicroChain::solved` の順。ミクロ燃焼しない位置は `None`）
    pub densities: &'a [Option<Vec<f64>>],
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
/// （ADR-0216。σᵢ はその位置のテーブルの使い方の状態と燃焼度で、N̄ᵢ⁰ は燃焼度で内挿した値）
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

/// テーブルの使い方ごとに、`state_of` の状態でテーブルを内挿して解く
fn solve_with(
    layout: &Layout,
    tables: &[Table],
    state_of: impl Fn(&TableUse) -> Option<State>,
    compositions: Option<&Compositions>,
    options: &Options,
) -> Result<Solved, CaseError> {
    if let Some(c) = compositions
        && c.densities.len() != layout.node_material.len()
    {
        return Err(CaseError::Depletion(format!(
            "{} micro compositions for {} positions",
            c.densities.len(),
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
            let mut constants = interpolate(table, state.as_ref(), t.rod.as_deref())
                .map_err(|e| CaseError::State(format!("{file}: {e}")))?;
            if let (Some(c), Some(site)) = (compositions, t.node)
                && let Some(Some(n)) = c.densities.get(site)
            {
                let micro = state
                    .as_ref()
                    .ok_or_else(|| {
                        CaseError::State(format!("{file}: micro depletion needs a state"))
                    })
                    .and_then(|s| {
                        interpolate_micro(table, s, t.rod.as_deref())
                            .map_err(|e| CaseError::State(format!("{file}: {e}")))
                    })?;
                // ミクロ燃焼のテーブルは分けない（格子全体の均質化。ADR-0216）
                let [f] = &mut constants[..] else {
                    return Err(CaseError::State(format!(
                        "{file}: micro data need an unsplit table"
                    )));
                };
                correct(f, &micro, &c.chain.tracked(n));
                // 補正で数密度の差が大きいと、吸収や全断面積が負になりうる。内挿の後と同じ検査をして、ソルバに渡す前に止める
                f.validate()
                    .map_err(|e| e.to_string())
                    .and_then(|()| non_negative(f))
                    .map_err(|e| {
                        CaseError::State(format!(
                            "{file}: the micro correction at position {site} gave invalid constants: {e}"
                        ))
                    })?;
            }
            Ok(constants)
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
