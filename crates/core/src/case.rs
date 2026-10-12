//! ケースの状態ごとの固有値計算（ADR-0146）と、臨界ホウ素濃度の探索（ADR-0145）。

use rayon::prelude::*;
use xlc_format::table::Table;

use crate::input::{Boron, CaseState, InputError, Layout};
use crate::nem::{NemError, Options, Problem, Solution, solve};
use crate::search::{SearchError, find};
use crate::state::{State, interpolate};

/// 状態を解いた結果。
#[derive(Debug, Clone)]
pub struct Solved {
    /// 探索したホウ素濃度（探索した状態だけ）
    pub searched: Option<Searched>,
    pub problem: Problem,
    pub solution: Solution,
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
}

/// テーブル（`layout.files` の順）を状態 `state`（`None` なら軸のないテーブル）で内挿し、炉心を解く。
pub fn solve_at(
    layout: &Layout,
    tables: &[Table],
    state: Option<&State>,
    options: &Options,
) -> Result<Solved, CaseError> {
    if tables.len() != layout.files.len() {
        return Err(CaseError::State(format!(
            "{} tables for {} few-group constant files",
            tables.len(),
            layout.files.len()
        )));
    }
    let constants = tables
        .iter()
        .zip(&layout.files)
        .map(|(table, file)| {
            interpolate(table, state)
                .map_err(|e| CaseError::State(format!("{}: {e}", file.display())))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let problem = layout.problem(&constants)?;
    let solution = solve(&problem, options)?;
    Ok(Solved {
        searched: None,
        problem,
        solution,
    })
}

/// 探索の状態の、解く前にできる検査。
pub fn check_search(state: &CaseState, options: &Options) -> Result<(), CaseError> {
    if let Boron::Search(search) = &state.boron
        // 固有値の反復の判定より細かい許容は区別できない（判定を超えても k の誤差を抑える保証ではない。ADR-0145）
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
    let search = match &state.boron {
        Boron::Given(ppm) => return solve_at(layout, tables, Some(&state.at(*ppm)), options),
        Boron::Search(search) => search,
    };
    check_search(state, options)?;
    let found = find(
        search.bracket,
        search.k,
        search.tolerance,
        search.max_solves,
        |ppm| {
            let solved = solve_at(layout, tables, Some(&state.at(ppm)), options)?;
            tracing::debug!(ppm, k = solved.solution.k, "boron search");
            Ok((solved.solution.k, solved))
        },
    )
    .map_err(|e| match e {
        SearchError::Solve(e) => e,
        e => CaseError::Search(e.to_string()),
    })?;
    Ok(Solved {
        searched: Some(Searched {
            ppm: found.x,
            solves: found.solves,
        }),
        ..found.at
    })
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
