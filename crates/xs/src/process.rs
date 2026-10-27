//! ENDF の1核種から多群定数を作る処理の全体。

use std::collections::BTreeMap;

use crate::endf::{
    Mf2, Mf3, Mf4, Mf5, Mf6, Mf7Incoherent, OrbitalL, RangeKind, Tab1, Tape, final_state_table,
    final_states, parse_delayed_nubar, parse_nubar,
};
use crate::multigroup::{
    Emission, FissionNeutrons, FreeGas, ReactionAverages, Reactions, SabError, SabOutside,
    ScatterError, ShieldingError, SlowingDown, Thermal, ThermalSab, UnresolvedShielding, Weight,
    emits_neutrons, fission_spectrum, reaction_averages, scattering_moments, shielded_averages,
    uses_angular,
};
use crate::parallel::in_parallel;
use crate::physics::BK;
use crate::pointwise::{
    ColdParts, DEFAULT_TOLERANCE, left_limit, linearize, resonance_tables, right_limit,
};
use crate::resonance::{
    Sampling, Sigma, UnresolvedXs, UrrError, unresolved_grid, unresolved_interpolation,
};
use xlc_format::group::GroupStructure;
use xlc_format::library::{
    CHAIN_REACTIONS, ChannelState, GroupChannel, GroupReaction, NuclideData, Pointwise, Resolved,
    Shielded, TemperatureData,
};

/// 処理の条件。
#[derive(Debug, Clone)]
pub struct Settings {
    pub groups: GroupStructure,
    pub weight: Weight,
    /// バックグラウンド断面積 (barn)、降順。最初の値を無限希釈とみなす
    pub sigma0: Vec<f64>,
    /// S(α, β) を持たない核種で、この入射エネルギー (eV) 未満の弾性散乱を自由ガス模型で扱う
    pub thermal_cutoff: f64,
    /// ポイントワイズの断面積の許容誤差
    pub tolerance: f64,
    /// ライブラリに持たせるポイントワイズの断面積の範囲 (eV)。`None` なら持たせない
    pub pointwise_range: Option<(f64, f64)>,
    /// 散乱行列のルジャンドルのモーメントの次数
    pub scattering_order: usize,
}

/// S(α, β) を使う入射エネルギーの上限。
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SabUpper {
    /// 評価済みデータの上限（MF7 MT4 の B(4)）
    Evaluated,
    /// 評価済みデータの上限以下のエネルギー (eV)
    Energy(f64),
}

/// S(α, β) を持つ核種の熱領域の扱い（ADR-0133）。
#[derive(Debug, Clone)]
pub struct ThermalScattering {
    /// MF7 MT4 の S(α, β)
    pub data: Mf7Incoherent,
    pub upper: SabUpper,
    pub outside: SabOutside,
    /// S(α, β) の上限から、この入射エネルギー (eV) 未満の弾性散乱を自由ガス模型で扱う
    pub free_gas_upper: f64,
}

impl ThermalScattering {
    /// S(α, β) と自由ガスの上限 (eV)。S(α, β) の上限は評価済みデータの上限以下で正、自由ガスの上限は有限で
    /// S(α, β) の上限以上であること。
    pub fn limits(&self) -> Result<(f64, f64), ProcessError> {
        let evaluated = self.data.e_max;
        let upper = match self.upper {
            SabUpper::Evaluated => evaluated,
            SabUpper::Energy(e) if e > 0.0 && e <= evaluated => e,
            SabUpper::Energy(e) => {
                return Err(ProcessError::SabUpper {
                    given: e,
                    evaluated,
                });
            }
        };
        if !(self.free_gas_upper.is_finite() && self.free_gas_upper >= upper) {
            return Err(ProcessError::FreeGasUpper {
                given: self.free_gas_upper,
                sab: upper,
            });
        }
        Ok((upper, self.free_gas_upper))
    }
}

/// 核種の処理のエラー。
#[derive(Debug, thiserror::Error)]
pub enum ProcessError {
    #[error(transparent)]
    Endf(#[from] crate::endf::Error),
    #[error(transparent)]
    Unresolved(#[from] UrrError),
    #[error(transparent)]
    Shielding(#[from] ShieldingError),
    #[error(transparent)]
    Sab(#[from] SabError),
    #[error(transparent)]
    Scatter(#[from] ScatterError),
    #[error("empty thermal scattering tape")]
    EmptyThermal,
    #[error("thermal elastic scattering (MF7 MT2) is not supported")]
    ThermalElastic,
    #[error("thermal scattering tape has no MF7 MT4")]
    NoIncoherent,
    #[error("S(α, β) limit {given} eV must be positive and at most the evaluated {evaluated} eV")]
    SabUpper { given: f64, evaluated: f64 },
    #[error("free gas limit {given} eV must be finite and at least the S(α, β) limit {sab} eV")]
    FreeGasUpper { given: f64, sab: f64 },
    #[error("more than one unresolved range")]
    MultipleUnresolved,
    #[error(
        "the R-matrix limited range has no s-wave elastic channel for the potential scattering"
    )]
    NoSWaveRadius,
    #[error("the scattering radius table does not cover {energy} eV")]
    RadiusTable { energy: f64 },
    #[error("s-wave scattering radius {0} must be positive")]
    SWaveRadius(f64),
    #[error("s-wave elastic channels have different radii {0:?}")]
    SWaveRadii(Vec<f64>),
    #[error("expected one material, found {0:?}")]
    Materials(Vec<i32>),
    #[error("MF3 MT{0} not found")]
    MissingMf3(u32),
    #[error("no MF3 data")]
    NoMf3,
    #[error("fission needs MF1 MT452 and MF5 MT18")]
    FissionData,
    #[error("the delayed neutron spectrum (MF5 MT455) needs MF1 {missing}")]
    DelayedFission { missing: &'static str },
    #[error("MF8/MT{mt} mixes MF9 and MF10 in its final states")]
    MixedFinalStates { mt: u32 },
    #[error(
        "the final states listed for MT{mt} exceed the reaction by {excess} b in group {group}"
    )]
    NegativeRemainder { mt: u32, group: usize, excess: f64 },
    #[error("the final state table of MT{mt} has no reaction cross section to follow")]
    EmptyProduct { mt: u32 },
    #[error("{temperature} K")]
    Temperature {
        temperature: f64,
        source: Box<ProcessError>,
    },
}

/// 熱中性子散乱のテープ（MF7）から、非干渉性非弾性散乱の S(α, β)（MT4）を読む。
///
/// 弾性散乱（MT2）を持つデータはエラーにする。S(α, β) の散乱行列は cutoff 未満で MF3 の弾性散乱を
/// 置き換えるので、MT2 を扱わないと弾性散乱が黙って抜け落ちる（黒鉛、ZrH など）。
pub fn read_thermal_scattering(tape: &Tape) -> Result<Mf7Incoherent, ProcessError> {
    let mat = *tape.materials().first().ok_or(ProcessError::EmptyThermal)?;
    if tape.section(mat, 7, 2).is_some() {
        return Err(ProcessError::ThermalElastic);
    }
    let section = tape.section(mat, 7, 4).ok_or(ProcessError::NoIncoherent)?;
    Ok(Mf7Incoherent::parse(&section)?)
}

/// 非分離共鳴の自己遮蔽の補正を、温度 `temperature` (K) とバックグラウンド断面積 `sigma0` について作る
/// （NJOY UNRESR に相当）。非分離共鳴がなければ `None`。
#[tracing::instrument(level = "debug", skip_all)]
pub fn unresolved_shielding(
    tape: &Tape,
    mat: i32,
    temperature: f64,
    sigma0: &[f64],
) -> Result<Option<UnresolvedShielding>, ProcessError> {
    let Some(section) = tape.section(mat, 2, 151) else {
        return Ok(None);
    };
    let mf2 = Mf2::parse(&section)?;
    let mut unresolved = None;
    for range in &mf2.ranges {
        match &range.kind {
            RangeKind::Unresolved(urr) if unresolved.is_none() => unresolved = Some((range, urr)),
            RangeKind::Unresolved(_) => return Err(ProcessError::MultipleUnresolved),
            _ => {}
        }
    }
    let Some((range, urr)) = unresolved else {
        return Ok(None);
    };
    let xs = UnresolvedXs::new(urr)?;
    let law = unresolved_interpolation(urr)?;
    let grid = unresolved_grid(urr, range.el, range.eh);

    // MF3 のバックグラウンド断面積（LSSF=0 のとき）。範囲の上端では内側（左側）の値を使う
    let background_table = |mt: u32| -> Result<Option<Tab1>, ProcessError> {
        if urr.lssf != 0 {
            return Ok(None);
        }
        match tape.section(mat, 3, mt) {
            Some(section) => Ok(Some(linearize(
                &Mf3::parse(&section)?.sigma,
                DEFAULT_TOLERANCE,
            ))),
            None => Ok(None),
        }
    };
    let tables = [1, 2, 18, 102].map(background_table);
    let [total, elastic, fission, capture] = tables;
    let (total, elastic, fission, capture) = (total?, elastic?, fission?, capture?);
    let at = |t: &Option<Tab1>, e: f64| {
        t.as_ref().map_or(0.0, |t| {
            if e < range.eh {
                right_limit(t, e)
            } else {
                left_limit(t, e)
            }
        })
    };

    // 格子点ごとに独立なので並列に求める。種は格子点の番号（温度と σ0 によらず同じラダーを使い、温度の間の相関を保つ。
    // スレッドの数によらない。ADR-0177）
    let sampled = in_parallel(grid.len(), |i| {
        let e = grid[i];
        let background = Sigma {
            total: at(&total, e),
            elastic: at(&elastic, e),
            fission: at(&fission, e),
            capture: at(&capture, e),
        };
        xs.shielded(
            e,
            temperature,
            sigma0,
            background,
            Sampling::LIBRARY,
            i as u64,
        )
    })?;
    let floored: usize = sampled.iter().map(|b| b.floored).sum();
    if floored > 0 {
        let samples = Sampling::LIBRARY.ladders * Sampling::LIBRARY.samples * grid.len();
        tracing::warn!(
            "{temperature} K: the sampled elastic (with the MF3 background) is negative at {floored} of \
             {samples} sampled energies in the unresolved range and is raised to 1 microbarn, as NJOY PURR does"
        );
    }
    let mut shielded = vec![Vec::with_capacity(grid.len()); sigma0.len()];
    let mut infinite = Vec::with_capacity(grid.len());
    for b in &sampled {
        for (row, s) in shielded.iter_mut().zip(&b.shielded) {
            row.push([s.total, s.elastic, s.fission, s.capture]);
        }
        // 無限希釈の欄は同じサンプルの平均にして、補正（LSSF=0 の差、LSSF=1 の比）がサンプルの揺らぎの 1 次を打ち消すようにする
        let mean = b.mean;
        infinite.push([mean.total, mean.elastic, mean.fission, mean.capture]);
    }
    Ok(Some(UnresolvedShielding::new(
        range.el, range.eh, urr.lssf, grid, law, shielded, infinite,
    )?))
}

/// 分離共鳴の範囲を持つ核種の、σ0 ごとの群定数の重みを減速計算で作る範囲（ADR-0072）。
/// 重み関数の熱中性子の区切り以上で最も低い群の境界から、分離共鳴の上端以下で最も高い群の境界まで
/// （熱中性子の領域では標的の熱運動が効き、静止した標的の減速計算は成り立たない）。範囲より上の散乱源の
/// ポテンシャル散乱は、分離共鳴の範囲の s 波の散乱半径から求める（`resolved`）。
/// 分離共鳴の範囲がなければ `None`。
pub fn slowing_down_range(
    tape: &Tape,
    mat: i32,
    groups: &GroupStructure,
    weight: &Weight,
) -> Result<Option<SlowingDown>, ProcessError> {
    let Some(section) = tape.section(mat, 2, 151) else {
        return Ok(None);
    };
    let mf2 = Mf2::parse(&section)?;
    let Some(Resolved { top, potential }) = resolved_in(&mf2)? else {
        return Ok(None);
    };
    let awr = mf2.awr;
    let bounds = groups.bounds();
    let [thermal, _] = weight.breakpoints();
    let (Some(&low), Some(&high)) = (
        bounds.iter().rev().find(|&&b| b >= thermal),
        bounds.iter().find(|&&b| b <= top),
    ) else {
        return Ok(None);
    };
    Ok((high > low).then_some(SlowingDown {
        awr,
        potential,
        range: (low, high),
    }))
}

/// 範囲 [low, high] (eV) のポイントワイズの断面積を、全断面積の格子で取り出す。全断面積は部分反応の和で作り直してあるので、
/// その格子は全反応の格子点と不連続点（x の重複）を含む。不連続点では左と右の極限を並べ、範囲の端では範囲の内側の
/// 極限を使う。νσf は ν(E)σf(E)（MT18 がなければ核分裂は 0）。
pub fn pointwise_within(
    tables: &BTreeMap<u32, Tab1>,
    nubar: Option<&Tab1>,
    (low, high): (f64, f64),
) -> Result<Pointwise, ProcessError> {
    let total = tables.get(&1).ok_or(ProcessError::MissingMf3(1))?;
    let elastic = tables.get(&2).ok_or(ProcessError::MissingMf3(2))?;
    let fission = tables.get(&18);
    let capture = tables.get(&102);
    // (エネルギー, 右の極限か)
    let mut points = vec![(low, true)];
    let inside: Vec<f64> = total
        .x
        .iter()
        .copied()
        .filter(|&e| e > low && e < high)
        .collect();
    for (i, &e) in inside.iter().enumerate() {
        // 同じエネルギーが続けば1つ目は左、2つ目は右の極限
        let repeated_before = i > 0 && inside[i - 1] == e;
        points.push((e, repeated_before));
    }
    points.push((high, false));
    let value = |t: &Tab1, (e, right): (f64, bool)| {
        if right {
            right_limit(t, e)
        } else {
            left_limit(t, e)
        }
    };
    let column = |t: &Tab1| points.iter().map(|&p| value(t, p)).collect::<Vec<f64>>();
    let fission_values = fission.map_or(vec![0.0; points.len()], column);
    let nu_fission = match nubar {
        Some(nu) => points
            .iter()
            .zip(&fission_values)
            .map(|(&(e, _), f)| nu.eval(e).unwrap_or(0.0) * f)
            .collect(),
        None => vec![0.0; points.len()],
    };
    Ok(Pointwise {
        energy: points.iter().map(|&(e, _)| e).collect(),
        total: column(total),
        elastic: column(elastic),
        fission: fission_values,
        nu_fission,
        capture: capture.map_or(vec![0.0; points.len()], column),
    })
}

/// 分離共鳴の上端と、ポテンシャル散乱断面積 4π a² (barn)。上端は分離共鳴の範囲の上端のうち最も高いもの。
/// a はその範囲の、s 波の位相シフトに使う散乱半径（Reich-Moore は l=0 の APL、0 なら AP。MLBW は AP、エネルギーの表が
/// あれば範囲の上端の値（範囲より上へつなぐ。ADR-0167）。R-matrix limited は弾性散乱の l=0 のチャネルの APE）。
/// 分離共鳴がなければ `None`。
pub fn resolved(tape: &Tape, mat: i32) -> Result<Option<Resolved>, ProcessError> {
    let Some(section) = tape.section(mat, 2, 151) else {
        return Ok(None);
    };
    resolved_in(&Mf2::parse(&section)?)
}

/// 読んだ MF2 から `resolved` を求める。
fn resolved_in(mf2: &Mf2) -> Result<Option<Resolved>, ProcessError> {
    let Some(range) = mf2
        .ranges
        .iter()
        .filter(|r| r.kind.is_resolved())
        .max_by(|a, b| a.eh.total_cmp(&b.eh))
    else {
        return Ok(None);
    };
    let radius = s_wave_radius(&range.kind, range.eh)?;
    Ok(Some(Resolved {
        top: range.eh,
        potential: 4.0 * std::f64::consts::PI * radius * radius,
    }))
}

/// 上端が `top` (eV) の分離共鳴の範囲の、s 波の位相シフトに使う散乱半径 (10^-12 cm)。エネルギーによるなら上端の値。
fn s_wave_radius(kind: &RangeKind, top: f64) -> Result<f64, ProcessError> {
    let radius = match kind {
        RangeKind::ReichMoore(rm) => rm
            .l_values
            .iter()
            .find(|lv| lv.l == OrbitalL::S && lv.apl != 0.0)
            .map_or(rm.ap, |lv| lv.apl),
        // 表は範囲を覆い、不連続を持たない（MF2 を読むときに確かめる）
        RangeKind::Mlbw(bw) => bw
            .scattering_radius(top)
            .ok_or(ProcessError::RadiusTable { energy: top })?,
        RangeKind::RMatrixLimited(rml) => {
            // l=0 の弾性散乱のチャネルはどれも同じ半径であること（違えば、どれを使うかの根拠がない）
            let mut radii: Vec<f64> = rml
                .spin_groups
                .iter()
                .flat_map(|g| &g.channels)
                .filter(|c| c.l == OrbitalL::S && rml.pairs[c.pair].mt == 2)
                .map(|c| c.ape)
                .collect();
            radii.sort_by(f64::total_cmp);
            radii.dedup();
            match radii[..] {
                [] => return Err(ProcessError::NoSWaveRadius),
                [radius] => radius,
                _ => return Err(ProcessError::SWaveRadii(radii)),
            }
        }
        RangeKind::ScatteringRadius { .. } | RangeKind::Unresolved(_) => {
            return Err(ProcessError::NoSWaveRadius);
        }
    };
    if !(radius > 0.0 && radius.is_finite()) {
        return Err(ProcessError::SWaveRadius(radius));
    }
    Ok(radius)
}

/// 熱領域の散乱の扱い（自由ガス、S(α, β)）による弾性散乱の群断面積の変化。
///
/// 弾性散乱の散乱行列の行の和は、熱領域の扱いを含む実効的な散乱断面積になる。無限希釈の弾性散乱との
/// 差を、全 σ0 の全断面積と弾性散乱に足して、全断面積 = 吸収 + 散乱 の関係を保つ。
/// 熱領域の外では、行の和は無限希釈の弾性散乱に等しいので、差は 0 になる。
pub fn thermal_correction(elastic: &[f64], elastic_matrix: &[Vec<f64>]) -> Vec<f64> {
    elastic
        .iter()
        .zip(elastic_matrix)
        .map(|(el, row)| row.iter().sum::<f64>() - el)
        .collect()
}

/// テープ中の（ただ1つの）物質を、温度 `temperatures` (K) ごとに処理する。
/// `thermal` があれば、その上限まで弾性散乱を S(α, β) と自由ガス模型で、なければ `settings.thermal_cutoff` 未満を
/// 自由ガス模型で扱う。
#[tracing::instrument(name = "process", level = "debug", skip_all)]
pub fn process_nuclide(
    tape: &Tape,
    name: &str,
    temperatures: &[f64],
    thermal: Option<&ThermalScattering>,
    settings: &Settings,
) -> Result<NuclideData, ProcessError> {
    if let Some(thermal) = thermal {
        thermal.limits()?;
        // NJOY THERMR（参照解の OpenMC の S(α, β) を作る）は ln S を 3 点の 2 次式で補間する（thermr.f90 の terpq）。
        // 評価済みデータの補間則に従い、差は知らせる（ADR-0133）
        tracing::warn!(
            nuclide = name,
            "S(α, β) follows the ENDF interpolation law (ln S linear); NJOY THERMR interpolates ln S \
             quadratically, so its scattering cross section differs (ADR-0133)"
        );
    }
    let evaluation = read_evaluation(tape, settings.tolerance)?;
    let common = Common {
        tape,
        evaluation: &evaluation,
        thermal,
        settings,
    };
    // 温度ごとに独立なので並列に求める。温度の処理のうち格子点の並列（`unresolved_shielding`）の外の部分を、ほかの温度と
    // 重ねる。結果は温度の順で、エラーは最初の温度のもの（順に求めたときと同じ）
    let temperatures = in_parallel(temperatures.len(), |k| {
        let t = temperatures[k];
        common
            .at_temperature(t)
            .map_err(|source| ProcessError::Temperature {
                temperature: t,
                source: Box::new(source),
            })
    })?;
    Ok(NuclideData {
        name: name.to_string(),
        za: evaluation.za,
        awr: evaluation.awr,
        resolved: resolved(tape, evaluation.mat)?,
        thermal_scattering: thermal.is_some(),
        temperatures,
    })
}

/// 燃焼チェーンの反応（`CHAIN_REACTIONS`。ADR-0164）の、経路ごとの無限希釈の群の断面積（重み関数だけで平均する。ADR-0172）。
/// 反応の断面積は `tables` の代表の MT（部分の MT があれば、再構成が部分の和で作り直したもの。`ColdParts`。部分の MT はチェーンの
/// `ChainReaction` と同じ）。MF8 が状態を MF9（多重度）か MF10（生成断面積）で
/// 挙げれば状態ごとに分け（テープの順）、基底状態を挙げなければ、残りを基底状態への経路（LFS が nil）にする（チェーンと同じ）。
/// MF9 の表の外の多重度は 0（NJOY GROUPR と同じ。反応がそこにあれば警告する）。挙げた経路の和が、経路の表の線形化の
/// 許容誤差 `tolerance`（相対）より反応を超える群があればエラー、基底状態を挙げて和が反応に足りなければ警告する。`tables` は
/// 温度に広げた断面積で、MF10 は 0 K の値を、反応の断面積が 0 でないエネルギーでだけ使う（しきい値より下を数えない）。
fn chain_reactions(
    tape: &Tape,
    mat: i32,
    groups: &GroupStructure,
    weight: &Weight,
    tables: &BTreeMap<u32, Tab1>,
    tolerance: f64,
) -> Result<Vec<GroupReaction>, ProcessError> {
    // 反応と状態の経路の表を番号で 1 つの表にまとめ、1 度に平均する
    let mut parts: BTreeMap<u32, Tab1> = BTreeMap::new();
    let mut layout = Vec::new();
    for spec in &CHAIN_REACTIONS {
        let Some(sigma) = tables.get(&spec.mt) else {
            continue;
        };
        let entries = final_states(tape, mat, spec.mt)?.unwrap_or_default();
        let mut keys = Vec::with_capacity(entries.len());
        for entry in &entries {
            let table = final_state_table(tape, mat, entry.mf, spec.mt, entry.lfs)?;
            let state = state_table(&table, entry.mf, sigma, tolerance, spec.mt, entry.lfs)?;
            let key = u32::try_from(parts.len()).unwrap_or(u32::MAX);
            parts.insert(key, state);
            keys.push(key);
        }
        let key = u32::try_from(parts.len()).unwrap_or(u32::MAX);
        parts.insert(key, sigma.clone());
        layout.push((spec.mt, entries, key, keys));
    }
    let averaged = reaction_averages(groups, weight, &parts, None).sigma;
    let mut reactions = Vec::with_capacity(layout.len());
    for (mt, entries, key, keys) in layout {
        let total = averaged[&key].clone();
        if entries.is_empty() {
            reactions.push(GroupReaction {
                mt,
                channels: vec![GroupChannel {
                    final_state: None,
                    cross_section: total,
                }],
            });
            continue;
        }
        let mut channels: Vec<GroupChannel> = entries
            .iter()
            .zip(&keys)
            .map(|(entry, key)| GroupChannel {
                final_state: Some(ChannelState {
                    mf: entry.mf,
                    lfs: Some(entry.lfs),
                }),
                cross_section: averaged[key].clone(),
            })
            .collect();
        let remainder: Vec<f64> = (0..groups.len())
            .map(|g| total[g] - channels.iter().map(|c| c.cross_section[g]).sum::<f64>())
            .collect();
        if let Some((group, excess)) = remainder
            .iter()
            .zip(&total)
            .enumerate()
            .find(|(_, (r, t))| **r < -tolerance * t.abs())
            .map(|(g, (r, _))| (g, -r))
        {
            return Err(ProcessError::NegativeRemainder { mt, group, excess });
        }
        if entries.iter().all(|e| e.excitation > 0.0) {
            let mf = entries[0].mf;
            if entries.iter().any(|e| e.mf != mf) {
                return Err(ProcessError::MixedFinalStates { mt });
            }
            // 許容誤差の内の負は、線形化の誤差なので 0 にする
            channels.push(GroupChannel {
                final_state: Some(ChannelState { mf, lfs: None }),
                cross_section: remainder.iter().map(|r| r.max(0.0)).collect(),
            });
        } else if let Some((group, missing)) = remainder
            .iter()
            .zip(&total)
            .enumerate()
            .find(|(_, (r, t))| **r > tolerance * t.abs())
            .map(|(g, (r, _))| (g, *r))
        {
            tracing::warn!(
                mt,
                group,
                missing,
                "the final states listed in MF8 leave part of the reaction to no channel (b)"
            );
        }
        reactions.push(GroupReaction { mt, channels });
    }
    Ok(reactions)
}

/// 状態の経路の断面積を、許容誤差 `tolerance` の lin-lin の表にする。MF9 なら多重度 `table` × 反応の断面積 `sigma`、MF10 なら
/// 生成断面積 `table` そのもので、どちらも `sigma` が 0 のエネルギーでは 0（MF3 と同じく運動学のしきい値から数える。ADR-0100）。
/// 表の外は 0 とし、MF9 の表の外に反応があれば警告する（NJOY GROUPR と同じく、その分は状態に入らない）。
fn state_table(
    table: &Tab1,
    mf: u32,
    sigma: &Tab1,
    tolerance: f64,
    mt: u32,
    lfs: u32,
) -> Result<Tab1, ProcessError> {
    let (Some(&low), Some(&high)) = (sigma.x.first(), sigma.x.last()) else {
        return Err(ProcessError::EmptyProduct { mt });
    };
    let uncovered = std::cell::Cell::new(None::<f64>);
    let value = |e: f64| -> Result<[f64; 1], std::convert::Infallible> {
        let s = sigma.eval(e).unwrap_or(0.0);
        if s == 0.0 {
            return Ok([0.0]);
        }
        match table.eval(e) {
            Some(v) if mf == 9 => Ok([v * s]),
            Some(v) => Ok([v]),
            None => {
                if mf == 9 && uncovered.get().is_none() {
                    uncovered.set(Some(e));
                }
                Ok([0.0])
            }
        }
    };
    let nodes: Vec<f64> = table.x.iter().chain(&sigma.x).copied().collect();
    let Ok(tables) = resonance_tables(value, &nodes, low, high, tolerance);
    if let Some(energy) = uncovered.get() {
        tracing::warn!(
            mt,
            lfs,
            energy,
            "the MF9 multiplicity does not cover all energies where the reaction is nonzero; the state \
             gets none of the reaction there (as NJOY GROUPR)"
        );
    }
    let [table] = &tables[..] else {
        return Err(ProcessError::EmptyProduct { mt });
    };
    Ok(table.clone())
}

/// 温度によらない、1核種の評価の読み込み（多群の処理の前に `process_nuclide` が読むもの）。
#[derive(Debug, Clone)]
pub struct Evaluation {
    pub mat: i32,
    /// 最初の MF3 の見出しの ZA と AWR
    pub za: f64,
    pub awr: f64,
    pub nubar: Option<Tab1>,
    /// 0 K の部分反応（許容誤差 `tolerance` で再構成したもの）
    pub cold: ColdParts,
    /// 核分裂で生まれる中性子の成分（χ の材料。`fission_parts`）
    pub fission_parts: Vec<FissionNeutrons>,
    /// 離散準位の非弾性散乱（MT51〜90）の Q 値
    pub q_values: BTreeMap<u32, f64>,
    /// 中性子を放出する反応（弾性散乱と核分裂を除く）の出射の分布
    pub emissions: BTreeMap<u32, Emission>,
    /// 散乱行列に使う角度分布（MF4）
    pub angular: BTreeMap<u32, Mf4>,
}

/// テープ中の（ただ1つの）物質を読み、0 K のポイントワイズの断面積を許容誤差 `tolerance` で再構成する。
pub fn read_evaluation(tape: &Tape, tolerance: f64) -> Result<Evaluation, ProcessError> {
    let materials = tape.materials();
    let [mat] = materials[..] else {
        return Err(ProcessError::Materials(materials));
    };
    let mf3 = |mt: u32| -> Result<Mf3, ProcessError> {
        let section = tape
            .section(mat, 3, mt)
            .ok_or(ProcessError::MissingMf3(mt))?;
        Ok(Mf3::parse(&section)?)
    };
    let first = *tape.mts(mat, 3).first().ok_or(ProcessError::NoMf3)?;
    let head = mf3(first)?;
    let nubar = match tape.section(mat, 1, 452) {
        Some(section) => Some(parse_nubar(&section)?),
        None => None,
    };
    let cold = ColdParts::read(tape, mat, tolerance)?;
    let fission_parts = fission_parts(tape, mat, nubar.as_ref())?;
    let mut q_values = BTreeMap::new();
    for mt in tape
        .mts(mat, 3)
        .into_iter()
        .filter(|mt| (51..=90).contains(mt))
    {
        q_values.insert(mt, mf3(mt)?.qi);
    }
    let mut emissions = BTreeMap::new();
    for (mt, section) in tape.sections(mat, 6) {
        if mt != 2 && emits_neutrons(mt) {
            emissions.insert(mt, Emission::Mf6(Mf6::parse(&section)?));
        }
    }
    // MF6 のない反応は、MF4/MF5 で与えられている（MF5 のエネルギー分布だけを使う）
    for (mt, section) in tape.sections(mat, 5) {
        if mt != 2 && emits_neutrons(mt) && !emissions.contains_key(&mt) {
            emissions.insert(mt, Emission::Mf5(Mf5::parse(&section)?));
        }
    }
    // 角度分布（MF4）: 弾性散乱、離散準位（2体反応）と、MF5 で与えられた連続的な放出の分
    let mut angular = BTreeMap::new();
    for (mt, section) in tape.sections(mat, 4) {
        if uses_angular(mt, emissions.get(&mt)) {
            angular.insert(mt, Mf4::parse(&section)?);
        }
    }
    Ok(Evaluation {
        mat,
        za: head.za,
        awr: head.awr,
        nubar,
        cold,
        fission_parts,
        q_values,
        emissions,
        angular,
    })
}

/// χ を作る中性子の成分を、物質 `mat` の MF1 と MF5 から選ぶ。`nubar` は MF1 MT452。MF5 MT18 がなければ空。
///
/// - 遅発中性子の分布（MF5 MT455）があれば、即発（MF1 MT456 の ν_p と MF5 MT18）と遅発（MF1 MT455 の ν_d と
///   MF5 MT455。p_k(E) は族の割合）の 2 つ。k の固有値の計算では遅発中性子も核分裂源に入るので、その分布で生まれる
/// - 遅発中性子の分布がなければ、全数（MT452）が MF5 MT18 の分布で生まれる。MF1 MT455 があれば warn を出す（ENDF/B-VIII.0 の
///   Np-239、Am-242m、Cm-244。ほかに分布の出どころがない）。この場合 ν_d は使わないので MT455 は読まない
pub fn fission_parts(
    tape: &Tape,
    mat: i32,
    nubar: Option<&Tab1>,
) -> Result<Vec<FissionNeutrons>, ProcessError> {
    let Some(section) = tape.section(mat, 5, 18) else {
        return Ok(Vec::new());
    };
    let prompt_spectrum = Mf5::parse(&section)?;
    let delayed = tape.section(mat, 1, 455);
    let Some(delayed_spectrum) = tape.section(mat, 5, 455) else {
        if delayed.is_some() {
            tracing::warn!(
                mat,
                "MF1 MT455 has no spectrum in MF5 MT455; the delayed neutrons get the MF5 MT18 spectrum"
            );
        }
        return Ok(nubar
            .map(|total| FissionNeutrons {
                nubar: total.clone(),
                spectrum: prompt_spectrum,
            })
            .into_iter()
            .collect());
    };
    let delayed = delayed.ok_or(ProcessError::DelayedFission { missing: "MT455" })?;
    let prompt = tape
        .section(mat, 1, 456)
        .ok_or(ProcessError::DelayedFission { missing: "MT456" })?;
    Ok(vec![
        FissionNeutrons {
            nubar: parse_nubar(&prompt)?,
            spectrum: prompt_spectrum,
        },
        FissionNeutrons {
            nubar: parse_delayed_nubar(&delayed)?,
            spectrum: Mf5::parse(&delayed_spectrum)?,
        },
    ])
}

/// 温度によらない、1核種の処理の材料。
struct Common<'a> {
    tape: &'a Tape,
    evaluation: &'a Evaluation,
    thermal: Option<&'a ThermalScattering>,
    settings: &'a Settings,
}

impl Common<'_> {
    #[tracing::instrument(level = "debug", skip(self))]
    fn at_temperature(&self, temperature: f64) -> Result<TemperatureData, ProcessError> {
        let settings = self.settings;
        let (groups, weight) = (&settings.groups, &settings.weight);
        let ng = groups.len();
        let tables = self.evaluation.cold.at(temperature);
        let shielding = unresolved_shielding(
            self.tape,
            self.evaluation.mat,
            temperature,
            &settings.sigma0,
        )?;
        let slowing = slowing_down_range(self.tape, self.evaluation.mat, groups, weight)?;
        let averages = shielded_averages(
            groups,
            weight,
            &tables,
            self.evaluation.nubar.as_ref(),
            &settings.sigma0,
            shielding.as_ref(),
            slowing.as_ref(),
        )?;
        let parts = &self.evaluation.fission_parts;
        let chi = match (tables.get(&18), self.evaluation.nubar.as_ref()) {
            (None, _) => vec![0.0; ng],
            (Some(fission), Some(_)) if !parts.is_empty() => {
                fission_spectrum(groups, weight, fission, parts)
            }
            (Some(_), _) => return Err(ProcessError::FissionData),
        };

        // 弾性散乱だけ熱領域の扱いを使い、その実効断面積で全断面積と弾性散乱を補正する
        // S(α, β) の核種は、UFG の減速計算に使う S(α, β) の移行の表も作る（ポイントワイズの範囲があるときだけ。ADR-0205）
        let (thermal, sab_transfer) = match self.thermal {
            Some(thermal) => {
                let (upper, free_gas_upper) = thermal.limits()?;
                let sab = ThermalSab::new(&thermal.data, temperature, upper, thermal.outside)?;
                let transfer = settings
                    .pointwise_range
                    .map(|(low, _)| sab.ufg_transfer(low))
                    .transpose()?
                    .flatten();
                let thermal = Thermal::Sab {
                    sab,
                    free_gas: FreeGas {
                        kt: BK * temperature,
                        cutoff: free_gas_upper,
                    },
                };
                (thermal, transfer)
            }
            None => (
                Thermal::FreeGas(FreeGas {
                    kt: BK * temperature,
                    cutoff: settings.thermal_cutoff,
                }),
                None,
            ),
        };
        let elastic_table = tables.get(&2).ok_or(ProcessError::MissingMf3(2))?;
        let elastic_only = BTreeMap::from([(2, elastic_table.clone())]);
        let others: BTreeMap<_, _> = tables
            .iter()
            .filter(|(mt, _)| **mt != 2)
            .map(|(mt, t)| (*mt, t.clone()))
            .collect();
        let order = settings.scattering_order;
        let elastic = scattering_moments(
            groups,
            weight,
            self.evaluation.awr,
            &Reactions {
                tables: &elastic_only,
                q_values: &self.evaluation.q_values,
                emissions: &self.evaluation.emissions,
                angular: &self.evaluation.angular,
            },
            Some(&thermal),
            order,
        )?;
        let mut scatter = scattering_moments(
            groups,
            weight,
            self.evaluation.awr,
            &Reactions {
                tables: &others,
                q_values: &self.evaluation.q_values,
                emissions: &self.evaluation.emissions,
                angular: &self.evaluation.angular,
            },
            None,
            order,
        )?;
        for (matrix, el) in scatter.iter_mut().zip(&elastic) {
            for (row, el_row) in matrix.iter_mut().zip(el) {
                for (v, e) in row.iter_mut().zip(el_row) {
                    *v += e;
                }
            }
        }
        let zeros = vec![0.0; ng];
        let sigma =
            |avg: &ReactionAverages, mt: u32| avg.sigma.get(&mt).cloned().unwrap_or(zeros.clone());
        let correction = thermal_correction(&sigma(&averages[0], 2), &elastic[0]);
        let corrected = |values: Vec<f64>| -> Vec<f64> {
            values.iter().zip(&correction).map(|(v, c)| v + c).collect()
        };
        let shielded = averages
            .iter()
            .map(|avg| Shielded {
                total: corrected(sigma(avg, 1)),
                elastic: corrected(sigma(avg, 2)),
                fission: sigma(avg, 18),
                nu_fission: avg.nu_fission.clone().unwrap_or(zeros.clone()),
                capture: sigma(avg, 102),
            })
            .collect();
        let pointwise = match settings.pointwise_range {
            Some(range) => pointwise_within(&tables, self.evaluation.nubar.as_ref(), range)?,
            None => Pointwise::default(),
        };
        let reactions = chain_reactions(
            self.tape,
            self.evaluation.mat,
            groups,
            weight,
            &tables,
            settings.tolerance,
        )?;
        Ok(TemperatureData {
            temperature,
            chi,
            scatter,
            shielded,
            pointwise,
            reactions,
            sab_transfer,
        })
    }
}

#[cfg(test)]
mod tests;
