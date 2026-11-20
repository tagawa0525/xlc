//! `xs check` の検査: ENDF の評価を、`process_nuclide` が読むものの読み込みと形式の検査だけで確かめる（ADR-0266）。
//!
//! 線形化、共鳴の再構成、ドップラー拡がり、非分離共鳴のラダー、多群化は行わないので、それらの中で初めて止まる評価は通る
//! （MF8 の経路の和が反応を超える検査、共鳴の式の数値のエラーなど）。

use serde::Serialize;
use xlc_format::library::CHAIN_REACTIONS;

use crate::endf::{Mf2, Mf3, RangeKind, Tape, final_state_table, final_states, total_nubar};
use crate::multigroup::check_reactions;
use crate::pointwise::ColdParts;
use crate::process::{ProcessError, read_distributions, resolved};
use crate::resonance::{UnresolvedXs, unresolved_interpolation};

/// テープ中の（ただ1つの）物質を、輸送する核種として `process_nuclide` が読むものの読み込みと形式の検査で確かめる。
pub fn check_evaluation(tape: &Tape) -> Result<(), ProcessError> {
    let materials = tape.materials();
    let [mat] = materials[..] else {
        return Err(ProcessError::Materials(materials));
    };
    let first = *tape.mts(mat, 3).first().ok_or(ProcessError::NoMf3)?;
    Mf3::parse(&tape.section(mat, 3, first).ok_or(ProcessError::NoMf3)?)?;
    let nubar = total_nubar(tape, mat)?;
    let mts = ColdParts::check(tape, mat)?;
    let distributions = read_distributions(tape, mat, nubar.as_ref())?;
    resolved(tape, mat)?;
    check_unresolved(tape, mat)?;
    // `process_nuclide` の温度ごとの処理（`at_temperature`）が、χ と遅発中性子の前に確かめる条件
    let fission = mts.contains(&18);
    if fission && (nubar.is_none() || distributions.fission_parts.is_empty()) {
        return Err(ProcessError::FissionData);
    }
    if distributions.delayed.is_some() && !fission {
        return Err(ProcessError::DelayedWithoutFission);
    }
    check_reactions(
        &mts,
        &distributions.q_values,
        &distributions.emissions,
        &distributions.angular,
    )?;
    check_final_states(tape, mat, &mts)
}

/// 非分離共鳴の自己遮蔽（`unresolved_shielding`）が、ラダーを抽出する前に確かめる条件: 範囲が 1 つ以下、半径の指定（NAPS）、
/// 系列の補間則が揃うこと。
fn check_unresolved(tape: &Tape, mat: i32) -> Result<(), ProcessError> {
    let Some(section) = tape.section(mat, 2, 151) else {
        return Ok(());
    };
    let mf2 = Mf2::parse(&section)?;
    let mut unresolved = mf2.ranges.iter().filter_map(|r| match &r.kind {
        RangeKind::Unresolved(urr) => Some(urr),
        _ => None,
    });
    let Some(urr) = unresolved.next() else {
        return Ok(());
    };
    if unresolved.next().is_some() {
        return Err(ProcessError::MultipleUnresolved);
    }
    UnresolvedXs::new(urr)?;
    unresolved_interpolation(urr)?;
    Ok(())
}

/// チェーンの反応（`chain_reactions`）の MF8 と、それが挙げる MF9、MF10 の表を読む。表のある反応 `mts` だけを読む。
fn check_final_states(tape: &Tape, mat: i32, mts: &[u32]) -> Result<(), ProcessError> {
    for spec in CHAIN_REACTIONS.iter().filter(|spec| mts.contains(&spec.mt)) {
        let entries = final_states(tape, mat, spec.mt)?.unwrap_or_default();
        for entry in &entries {
            final_state_table(tape, mat, entry.mf, spec.mt, entry.lfs)?;
        }
        if let Some(first) = entries.first()
            && entries.iter().all(|e| e.excitation > 0.0)
            && entries.iter().any(|e| e.mf != first.mf)
        {
            return Err(ProcessError::MixedFinalStates { mt: spec.mt });
        }
    }
    Ok(())
}

/// ファイルの中身を、テープとして読んで `check_evaluation` で確かめる。止まれば理由の文。
pub fn check_contents(_bytes: &[u8]) -> Result<(), String> {
    Ok(())
}

/// 1 つの ENDF のファイルの検査の結果。
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Checked {
    /// 与えられたパス
    pub endf: String,
    /// 止まった理由（読めれば `None`）
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

/// `report_toml` が書く報告の中身。
#[derive(Serialize)]
struct Report<'a> {
    evaluation: &'a [Checked],
}

/// 検査の結果を、与えられた順に TOML の `[[evaluation]]` の並びで書く。
pub fn report_toml(checked: &[Checked]) -> Result<String, toml::ser::Error> {
    toml::to_string(&Report {
        evaluation: checked,
    })
}

#[cfg(test)]
mod tests;
