//! ENDF の1核種から多群定数を作る処理の全体。

use std::collections::BTreeMap;

use crate::endf::{
    Mf2, Mf3, Mf4, Mf5, Mf6, Mf7Incoherent, OrbitalL, RangeKind, Tab1, Tape, parse_nubar,
};
use crate::multigroup::{
    Emission, FreeGas, ReactionAverages, Reactions, SabError, SabOutside, ScatterError,
    ShieldingError, SlowingDown, Thermal, ThermalSab, UnresolvedShielding, Weight, emits_neutrons,
    fission_spectrum, scattering_moments, shielded_averages, uses_angular,
};
use crate::physics::BK;
use crate::pointwise::{ColdParts, DEFAULT_TOLERANCE, left_limit, linearize, right_limit};
use crate::resonance::{Sigma, UnresolvedXs, UrrError, unresolved_grid, unresolved_interpolation};
use xlc_format::group::GroupStructure;
use xlc_format::library::{NuclideData, Pointwise, Resolved, Shielded, TemperatureData};

/// 処理の条件。
#[derive(Debug, Clone)]
pub struct Settings {
    pub groups: GroupStructure,
    pub weight: Weight,
    /// 背景断面積 (barn)、降順。最初の値を無限希釈とみなす
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

/// 非分離共鳴の自己遮蔽の補正を、温度 `temperature` (K) と背景断面積 `sigma0` について作る
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

    // MF3 の背景断面積（LSSF=0 のとき）。範囲の上端では内側（左側）の値を使う
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

    let mut shielded = vec![Vec::with_capacity(grid.len()); sigma0.len()];
    let mut infinite = Vec::with_capacity(grid.len());
    for &e in &grid {
        let background = Sigma {
            total: at(&total, e),
            elastic: at(&elastic, e),
            fission: at(&fission, e),
            capture: at(&capture, e),
        };
        for (row, s) in shielded
            .iter_mut()
            .zip(xs.shielded(e, temperature, sigma0, background)?)
        {
            row.push([s.total, s.elastic, s.fission, s.capture]);
        }
        let average = xs.cross_sections(e)?;
        infinite.push([
            average.total + background.total,
            average.elastic + background.elastic,
            average.fission + background.fission,
            average.capture + background.capture,
        ]);
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
        RangeKind::Mlbw(bw) => bw.ap_energy.as_ref().map_or(bw.ap, |t| t.eval(top)),
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
    let cold = ColdParts::read(tape, mat, settings.tolerance)?;
    let mf5 = match tape.section(mat, 5, 18) {
        Some(section) => Some(Mf5::parse(&section)?),
        None => None,
    };
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
    let common = Common {
        tape,
        mat,
        awr: head.awr,
        nubar: nubar.as_ref(),
        mf5: mf5.as_ref(),
        q_values: &q_values,
        emissions: &emissions,
        angular: &angular,
        thermal,
        settings,
        cold: &cold,
    };
    let temperatures = temperatures
        .iter()
        .map(|&t| {
            common
                .at_temperature(t)
                .map_err(|source| ProcessError::Temperature {
                    temperature: t,
                    source: Box::new(source),
                })
        })
        .collect::<Result<_, _>>()?;
    Ok(NuclideData {
        name: name.to_string(),
        za: head.za,
        awr: head.awr,
        resolved: resolved(tape, mat)?,
        thermal_scattering: thermal.is_some(),
        temperatures,
    })
}

/// 温度によらない、1核種の処理の材料。
struct Common<'a> {
    tape: &'a Tape,
    mat: i32,
    awr: f64,
    nubar: Option<&'a Tab1>,
    mf5: Option<&'a Mf5>,
    q_values: &'a BTreeMap<u32, f64>,
    emissions: &'a BTreeMap<u32, Emission>,
    angular: &'a BTreeMap<u32, Mf4>,
    thermal: Option<&'a ThermalScattering>,
    settings: &'a Settings,
    /// 0 K の部分反応（温度ごとに広げる）
    cold: &'a ColdParts,
}

impl Common<'_> {
    #[tracing::instrument(level = "debug", skip(self))]
    fn at_temperature(&self, temperature: f64) -> Result<TemperatureData, ProcessError> {
        let settings = self.settings;
        let (groups, weight) = (&settings.groups, &settings.weight);
        let ng = groups.len();
        let tables = self.cold.at(temperature);
        let shielding = unresolved_shielding(self.tape, self.mat, temperature, &settings.sigma0)?;
        let slowing = slowing_down_range(self.tape, self.mat, groups, weight)?;
        let averages = shielded_averages(
            groups,
            weight,
            &tables,
            self.nubar,
            &settings.sigma0,
            shielding.as_ref(),
            slowing.as_ref(),
        )?;
        let chi = match (tables.get(&18), self.nubar, self.mf5) {
            (None, ..) => vec![0.0; ng],
            (Some(fission), Some(nu), Some(mf5)) => {
                fission_spectrum(groups, weight, fission, nu, mf5)
            }
            (Some(_), ..) => return Err(ProcessError::FissionData),
        };

        // 弾性散乱だけ熱領域の扱いを使い、その実効断面積で全断面積と弾性散乱を補正する
        let thermal = match self.thermal {
            Some(thermal) => {
                let (upper, free_gas_upper) = thermal.limits()?;
                Thermal::Sab {
                    sab: ThermalSab::new(&thermal.data, temperature, upper, thermal.outside)?,
                    free_gas: FreeGas {
                        kt: BK * temperature,
                        cutoff: free_gas_upper,
                    },
                }
            }
            None => Thermal::FreeGas(FreeGas {
                kt: BK * temperature,
                cutoff: settings.thermal_cutoff,
            }),
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
            self.awr,
            &Reactions {
                tables: &elastic_only,
                q_values: self.q_values,
                emissions: self.emissions,
                angular: self.angular,
            },
            Some(&thermal),
            order,
        )?;
        let mut scatter = scattering_moments(
            groups,
            weight,
            self.awr,
            &Reactions {
                tables: &others,
                q_values: self.q_values,
                emissions: self.emissions,
                angular: self.angular,
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
            })
            .collect();
        let pointwise = match settings.pointwise_range {
            Some(range) => pointwise_within(&tables, self.nubar, range)?,
            None => Pointwise::default(),
        };
        Ok(TemperatureData {
            temperature,
            chi,
            scatter,
            shielded,
            pointwise,
        })
    }
}

#[cfg(test)]
mod tests;
