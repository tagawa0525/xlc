//! ポイントワイズの断面積からの多群定数の作成（NJOY の GROUPR に相当）。

mod fission;
mod free_gas;
mod integrate;
mod outgoing;
mod sab;
mod scatter;
mod slowing;
mod weight;

use std::collections::BTreeMap;

use xlc_format::group::GroupStructure;

use crate::endf::{Interpolation, Tab1};

pub use fission::{FissionNeutrons, fission_spectrum};
pub use free_gas::{FreeGas, free_gas_group_integrals, free_gas_kernel};
pub use integrate::{group_integrals, group_integrals_multi};
pub use outgoing::{AnalyticSpectrum, GroupedLaw, GroupedOutgoing, outgoing_group_integrals};
#[cfg(test)]
pub(crate) use sab::tests as sab_tests;
pub use sab::{SabError, SabOutside, ThermalSab};
pub use scatter::{
    Emission, Reactions, ScatterError, Thermal, emits_neutrons, multiplicity, scattering_matrix,
    scattering_moments, two_body_range, uses_angular,
};
pub use slowing::{SlowingDownError, slowing_down_shape, slowing_down_shape_with};
pub use weight::Weight;

/// 群平均断面積。配列は群の順（高エネルギーから）。
#[derive(Debug, Clone, PartialEq)]
pub struct ReactionAverages {
    /// 群ごとの中性子束の積分 ∫_g φ dE（無限希釈なら φ は重み関数 w）
    pub flux: Vec<f64>,
    /// MT ごとの σ_g = ∫_g σ φ dE / ∫_g φ dE (barn)
    pub sigma: BTreeMap<u32, Vec<f64>>,
    /// νσf の群平均（核分裂がなければ `None`）
    pub nu_fission: Option<Vec<f64>>,
}

/// ポイントワイズの断面積 `tables`（lin-lin、キーは MT）を、重み関数 `weight` で群平均する。
/// νσf は ν(E)σf(E) を群平均する。σf には MT18 を使う。
pub fn reaction_averages(
    groups: &GroupStructure,
    weight: &Weight,
    tables: &BTreeMap<u32, Tab1>,
    nubar: Option<&Tab1>,
) -> ReactionAverages {
    averages_with_flux(groups, weight, tables, nubar, &[], |_| 1.0, |_, _, v| v)
}

/// 非分離共鳴の自己遮蔽によるポイントワイズの断面積の補正（NJOY GROUPR getunr と同じ）。
///
/// 範囲 [el, eh] で、全断面積・弾性散乱・核分裂・捕獲（MT1, 2, 18/19, 102）を、LSSF=1 なら
/// s(σ0)/s(∞) 倍し、LSSF=0 なら s(σ0) - s(∞) を足す。s(σ0) は σ0 ごとの実効断面積、s(∞) は
/// 無限希釈の値で、どちらも格子点の間を `law` で補間する。
#[derive(Debug, Clone, PartialEq)]
pub struct UnresolvedShielding {
    el: f64,
    eh: f64,
    multiplicative: bool,
    grid: Vec<f64>,
    law: Interpolation,
    /// [σ0][格子点] の [全, 弾性, 核分裂, 捕獲]
    shielded: Vec<Vec<[f64; 4]>>,
    /// [格子点] の無限希釈の値
    infinite: Vec<[f64; 4]>,
}

impl UnresolvedShielding {
    /// `shielded` は σ0 ごとに格子点の数だけ、`infinite` は格子点の数だけ値を持つこと。
    pub fn new(
        el: f64,
        eh: f64,
        lssf: i64,
        grid: Vec<f64>,
        law: Interpolation,
        shielded: Vec<Vec<[f64; 4]>>,
        infinite: Vec<[f64; 4]>,
    ) -> Result<Self, ShieldingError> {
        let multiplicative = match lssf {
            0 => false,
            1 => true,
            _ => return Err(ShieldingError::Lssf(lssf)),
        };
        let n = grid.len();
        if n < 2 || grid.windows(2).any(|w| w[1] <= w[0]) || el >= eh {
            return Err(ShieldingError::Grid);
        }
        if infinite.len() != n || shielded.iter().any(|row| row.len() != n) {
            return Err(ShieldingError::GridValues);
        }
        Ok(Self {
            el,
            eh,
            multiplicative,
            grid,
            law,
            shielded,
            infinite,
        })
    }

    /// 補正する範囲 [el, eh] (eV)。
    pub fn range(&self) -> (f64, f64) {
        (self.el, self.eh)
    }

    /// σ0 の番号 `s` で、MT `mt` のポイントワイズの断面積 `value`（エネルギー e）を補正した値。
    fn adjust(&self, s: usize, mt: u32, e: f64, value: f64) -> f64 {
        let component = match mt {
            1 => 0,
            2 => 1,
            18 | 19 => 2,
            102 => 3,
            _ => return value,
        };
        if e < self.el || e > self.eh {
            return value;
        }
        let grid = &self.grid;
        let k = grid.partition_point(|&g| g <= e).clamp(1, grid.len() - 1);
        let interpolate = |v: &[[f64; 4]]| {
            self.law.interpolate(
                e,
                (grid[k - 1], v[k - 1][component]),
                (grid[k], v[k][component]),
            )
        };
        let (shielded, infinite) = (interpolate(&self.shielded[s]), interpolate(&self.infinite));
        if !self.multiplicative {
            value + shielded - infinite
        } else if infinite != 0.0 {
            value * shielded / infinite
        } else {
            value
        }
    }
}

/// 自己遮蔽した群定数を作るときのエラー。
#[derive(Debug, thiserror::Error)]
pub enum ShieldingError {
    #[error("unresolved shielding needs LSSF 0 or 1, got {0}")]
    Lssf(i64),
    #[error("unresolved grid must have increasing energies in [el, eh]")]
    Grid,
    #[error("unresolved values must match the grid")]
    GridValues,
    #[error("self-shielding needs the total cross section (MT1)")]
    NoTotal,
    #[error("sigma zero must be positive and finite, got {0}")]
    Sigma0(f64),
    #[error("unresolved shielding has {found} sigma zeros, expected {expected}")]
    UnresolvedSigma0 { found: usize, expected: usize },
    #[error("the slowing-down weight needs elastic scattering (MT2)")]
    NoElastic,
    #[error("the slowing-down weight is not defined at {0} eV")]
    Undefined(f64),
    #[error(transparent)]
    SlowingDown(#[from] SlowingDownError),
}

/// 減速計算で重みを作る範囲。範囲の中では、超微細群の減速計算の中性子束で重みづけする（`slowing`）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlowingDown {
    /// 吸収体の質量比 AWR
    pub awr: f64,
    /// 吸収体のポテンシャル散乱断面積 (b)。範囲より上の散乱源（NR 近似）に使う
    pub potential: f64,
    /// 範囲 (eV)。上端は群の境界であること（群の中で重みの式が切り替わらないように）
    pub range: (f64, f64),
}

/// バックグラウンド断面積 σ0 ごとの群平均断面積。`sigma0` と同じ順に返す。
///
/// 中性子束は、狭い共鳴の近似（Bondarenko 法）の w(E)/(σt(E) + σ0) とする（σt は MT1）。NJOY GROUPR の
/// 中性子束はこれに定数 (σpot + σ0) を掛けたものだが、群平均では打ち消し合う。`slowing` を与えると、
/// その範囲では超微細群の減速計算の中性子束を使い、範囲の上では Bondarenko の中性子束に定数
/// (σp + σ0) を掛けて上端でつなぐ（定数なので範囲の上の群平均は変わらない）。範囲の下
/// （熱中性子の領域）では重み関数 w の定数倍（σ0 によらない。燃料の中の熱中性子束のへこみは輸送計算で扱う）。
#[tracing::instrument(level = "debug", skip_all)]
pub fn shielded_averages(
    groups: &GroupStructure,
    weight: &Weight,
    tables: &BTreeMap<u32, Tab1>,
    nubar: Option<&Tab1>,
    sigma0: &[f64],
    unresolved: Option<&UnresolvedShielding>,
    slowing: Option<&SlowingDown>,
) -> Result<Vec<ReactionAverages>, ShieldingError> {
    let total = tables.get(&1).ok_or(ShieldingError::NoTotal)?;
    if let Some(&bad) = sigma0.iter().find(|s| !(s.is_finite() && **s > 0.0)) {
        return Err(ShieldingError::Sigma0(bad));
    }
    if let Some(u) = unresolved
        && u.shielded.len() != sigma0.len()
    {
        return Err(ShieldingError::UnresolvedSigma0 {
            found: u.shielded.len(),
            expected: sigma0.len(),
        });
    }
    let mut breakpoints = total.x.clone();
    if let Some(u) = unresolved {
        breakpoints.extend(&u.grid);
    }
    sigma0
        .iter()
        .enumerate()
        .map(|(s, &s0)| {
            let adjust = |mt: u32, e: f64, value: f64| match unresolved {
                Some(u) => u.adjust(s, mt, e, value),
                None => value,
            };
            let bondarenko = |e: f64| 1.0 / (adjust(1, e, total.eval(e).unwrap_or(0.0)) + s0);
            let Some(slowing) = slowing else {
                return Ok(averages_with_flux(
                    groups,
                    weight,
                    tables,
                    nubar,
                    &breakpoints,
                    bondarenko,
                    adjust,
                ));
            };
            let elastic = tables.get(&2).ok_or(ShieldingError::NoElastic)?;
            let shape = slowing_down_shape(total, elastic, slowing, weight, s0)?;
            let flux = slowing_down_flux(&shape, slowing, s0, bondarenko)?;
            let points: Vec<f64> = breakpoints.iter().chain(&shape.x).copied().collect();
            Ok(averages_with_flux(
                groups, weight, tables, nubar, &points, flux, adjust,
            ))
        })
        .collect()
}

/// 減速計算の形 `shape` と、範囲の上の Bondarenko の中性子束 `bondarenko` = 1/(σt + σ0) をつないだ、重み関数に
/// 対する中性子束の比。範囲の上は (σp + σ0) を掛けて上端でつなぐ（減速計算の上端の値は
/// W (σp + σ0)/(σt + σ0)）。範囲の下は重み関数の定数倍（下端でつながるよう、下端の値を掛ける）。
fn slowing_down_flux<'a>(
    shape: &'a Tab1,
    slowing: &SlowingDown,
    sigma0: f64,
    bondarenko: impl Fn(f64) -> f64 + 'a,
) -> Result<impl Fn(f64) -> f64 + 'a, ShieldingError> {
    let (low, high) = slowing.range;
    let scale = slowing.potential + sigma0;
    let at_low = shape.eval(low).ok_or(ShieldingError::Undefined(low))?;
    Ok(move |e: f64| {
        if e < low {
            at_low
        } else if e <= high {
            // 表は範囲をちょうど覆うので、引けないのは誤り。NaN はライブラリを書く前の検査で止まる
            shape.eval(e).unwrap_or(f64::NAN)
        } else {
            scale * bondarenko(e)
        }
    })
}

/// 中性子束 w(E)·`shape`(E) による群平均。`shape` は `breakpoints` の間で滑らかであること。
/// 断面積は `adjust`(MT, E, ポイントワイズの断面積) で補正した値を使う。
fn averages_with_flux(
    groups: &GroupStructure,
    weight: &Weight,
    tables: &BTreeMap<u32, Tab1>,
    nubar: Option<&Tab1>,
    breakpoints: &[f64],
    shape: impl Fn(f64) -> f64,
    adjust: impl Fn(u32, f64, f64) -> f64,
) -> ReactionAverages {
    let flux = group_integrals(groups, weight, breakpoints, &shape);
    let average = |points: &[f64], f: &dyn Fn(f64) -> f64| -> Vec<f64> {
        let points: Vec<f64> = points.iter().chain(breakpoints).copied().collect();
        group_integrals(groups, weight, &points, |e| f(e) * shape(e))
            .iter()
            .zip(&flux)
            .map(|(integral, phi)| integral / phi)
            .collect()
    };
    let sigma = tables
        .iter()
        .map(|(&mt, t)| {
            (
                mt,
                average(&t.x, &|e| adjust(mt, e, t.eval(e).unwrap_or(0.0))),
            )
        })
        .collect();
    let nu_fission = match (nubar, tables.get(&18)) {
        (Some(nu), Some(fission)) => {
            let points: Vec<f64> = fission.x.iter().chain(&nu.x).copied().collect();
            Some(average(&points, &|e| {
                nu.eval(e).unwrap_or(0.0) * adjust(18, e, fission.eval(e).unwrap_or(0.0))
            }))
        }
        _ => None,
    };
    ReactionAverages {
        flux,
        sigma,
        nu_fission,
    }
}

#[cfg(test)]
mod tests;
