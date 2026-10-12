//! ケースの状態ごとの固有値計算（ADR-0146）。

use xlc_format::table::Table;

use crate::input::{CaseState, InputError, Layout};
use crate::nem::{NemError, Options, Problem, Solution};
use crate::state::State;

/// 状態を解いた結果。
#[derive(Debug, Clone)]
pub struct Solved {
    pub problem: Problem,
    pub solution: Solution,
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
}

/// テーブル（`layout.files` の順）を状態 `state`（`None` なら軸のないテーブル）で内挿し、炉心を解く。
pub fn solve_at(
    layout: &Layout,
    tables: &[Table],
    state: Option<&State>,
    options: &Options,
) -> Result<Solved, CaseError> {
    let _ = (layout, tables, state, options);
    todo!()
}

/// ケースの状態を並列に解く（結果は `states` の順）。各状態の解は1つのスレッドで解くので、スレッドの数によらない。
pub fn solve_states(
    layout: &Layout,
    tables: &[Table],
    states: &[CaseState],
    options: &Options,
) -> Vec<Result<Solved, CaseError>> {
    let _ = (layout, tables, states, options);
    todo!()
}

#[cfg(test)]
mod tests;
