//! ENDF の1核種から多群定数を作る処理の全体。

use std::collections::BTreeMap;

use crate::endf::{
    Mf2, Mf3, Mf5, Mf6, Mf7Incoherent, OrbitalL, RangeKind, Tab1, Tape, parse_nubar,
};
use crate::library::{NuclideData, Shielded, TemperatureData};
use crate::multigroup::{
    Emission, FreeGas, GroupStructure, ReactionAverages, SabError, ScatterError, ShieldingError,
    SlowingDown, Thermal, ThermalSab, UnresolvedShielding, Weight, emits_neutrons,
    fission_spectrum, scattering_matrix, shielded_averages,
};
use crate::physics::BK;
use crate::pointwise::{ColdParts, DEFAULT_TOLERANCE, left_limit, linearize, right_limit};
use crate::resonance::{Sigma, UnresolvedXs, UrrError, unresolved_grid, unresolved_interpolation};

/// 処理の条件。
#[derive(Debug, Clone)]
pub struct Settings {
    pub groups: GroupStructure,
    pub weight: Weight,
    /// 背景断面積 (barn)、降順。最初の値を無限希釈とみなす
    pub sigma0: Vec<f64>,
    /// この入射エネルギー (eV) 未満の弾性散乱を、自由ガス模型か S(α, β) で扱う
    pub thermal_cutoff: f64,
    /// ポイントワイズの断面積の許容誤差
    pub tolerance: f64,
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
    #[error("more than one unresolved range")]
    MultipleUnresolved,
    #[error(
        "the R-matrix limited range has no s-wave elastic channel for the potential scattering"
    )]
    NoSWaveRadius,
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

/// 分離共鳴の範囲を持つ核種の、σ0 ごとの群定数の重みを減速計算で作る範囲（docs/plans/m2-assembly.md）。
/// 重み関数の熱中性子の区切り以上で最も低い群の境界から、分離共鳴の上端以下で最も高い群の境界まで
/// （熱中性子の領域では標的の熱運動が効き、静止した標的の減速計算は成り立たない）。
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
    let Some(top) = mf2
        .ranges
        .iter()
        .filter(|r| r.kind.is_resolved())
        .map(|r| r.eh)
        .reduce(f64::max)
    else {
        return Ok(None);
    };
    let bounds = groups.bounds();
    let [thermal, _] = weight.breakpoints();
    let (Some(&low), Some(&high)) = (
        bounds.iter().rev().find(|&&b| b >= thermal),
        bounds.iter().find(|&&b| b <= top),
    ) else {
        return Ok(None);
    };
    Ok((high > low).then_some(SlowingDown {
        awr: mf2.awr,
        range: (low, high),
    }))
}

/// 分離共鳴を持つ核種のポテンシャル散乱断面積 4π a² (barn)。a は上端がいちばん高い分離共鳴の範囲の、
/// s 波の位相シフトに使う散乱半径（Reich-Moore は l=0 の APL、0 なら AP。MLBW は AP。R-matrix limited は
/// 弾性散乱の l=0 のチャネルの APE）。分離共鳴がなければ `None`。ESSM の線源と全断面積に使う。
pub fn potential_scattering(tape: &Tape, mat: i32) -> Result<Option<f64>, ProcessError> {
    let Some(section) = tape.section(mat, 2, 151) else {
        return Ok(None);
    };
    let mf2 = Mf2::parse(&section)?;
    let Some(range) = mf2
        .ranges
        .iter()
        .filter(|r| r.kind.is_resolved())
        .max_by(|a, b| a.eh.total_cmp(&b.eh))
    else {
        return Ok(None);
    };
    let radius = match &range.kind {
        RangeKind::ReichMoore(rm) => rm
            .l_values
            .iter()
            .find(|lv| lv.l == OrbitalL::S && lv.apl != 0.0)
            .map_or(rm.ap, |lv| lv.apl),
        RangeKind::Mlbw(bw) => bw.ap,
        RangeKind::RMatrixLimited(rml) => rml
            .spin_groups
            .iter()
            .flat_map(|g| &g.channels)
            .find(|c| c.l == OrbitalL::S && rml.pairs[c.pair].mt == 2)
            .map(|c| c.ape)
            .ok_or(ProcessError::NoSWaveRadius)?,
        RangeKind::ScatteringRadius { .. } | RangeKind::Unresolved(_) => return Ok(None),
    };
    Ok(Some(4.0 * std::f64::consts::PI * radius * radius))
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
/// `thermal` があれば cutoff 未満の弾性散乱を S(α, β) で、なければ自由ガス模型で扱う。
#[tracing::instrument(name = "process", level = "debug", skip_all)]
pub fn process_nuclide(
    tape: &Tape,
    name: &str,
    temperatures: &[f64],
    thermal: Option<&Mf7Incoherent>,
    settings: &Settings,
) -> Result<NuclideData, ProcessError> {
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
    let common = Common {
        tape,
        mat,
        awr: head.awr,
        nubar: nubar.as_ref(),
        mf5: mf5.as_ref(),
        q_values: &q_values,
        emissions: &emissions,
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
        potential: potential_scattering(tape, mat)?,
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
    thermal: Option<&'a Mf7Incoherent>,
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
            Some(data) => {
                Thermal::Sab(ThermalSab::new(data, temperature, settings.thermal_cutoff)?)
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
        let elastic = scattering_matrix(
            groups,
            weight,
            self.awr,
            &elastic_only,
            self.q_values,
            self.emissions,
            Some(&thermal),
        )?;
        let mut scatter = scattering_matrix(
            groups,
            weight,
            self.awr,
            &others,
            self.q_values,
            self.emissions,
            None,
        )?;
        for (row, el) in scatter.iter_mut().zip(&elastic) {
            for (v, e) in row.iter_mut().zip(el) {
                *v += e;
            }
        }
        let zeros = vec![0.0; ng];
        let sigma =
            |avg: &ReactionAverages, mt: u32| avg.sigma.get(&mt).cloned().unwrap_or(zeros.clone());
        let correction = thermal_correction(&sigma(&averages[0], 2), &elastic);
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
        Ok(TemperatureData {
            temperature,
            chi,
            scatter,
            shielded,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thermal_correction_is_row_sum_minus_elastic() {
        let elastic = vec![4.0, 20.0];
        let matrix = vec![vec![3.9, 0.1], vec![2.0, 27.0]];
        assert_eq!(thermal_correction(&elastic, &matrix), vec![0.0, 9.0]);
    }

    /// 非分離共鳴（LSSF=1、s 波の1系列）と、MF3 の一定の断面積だけのテープ。
    fn unresolved_tape() -> Tape {
        use crate::endf::rec;
        let mf3 = |mt: u32, value: &str| {
            [
                rec(["1.001000+3", "2.360060+2", "0", "0", "0", "0"], 1, 3, mt),
                rec(["0.0", "0.0", "0", "0", "1", "2"], 1, 3, mt),
                rec(["2", "2", "", "", "", ""], 1, 3, mt),
                rec(["1.0-5", value, "2.0+7", value, "", ""], 1, 3, mt),
                rec(["", "", "", "", "", ""], 1, 3, 0),
            ]
            .concat()
        };
        let text = [
            rec(["1.001000+3", "2.360060+2", "0", "0", "1", "0"], 1, 2, 151),
            rec(["1.001000+3", "1.0", "0", "0", "1", "0"], 1, 2, 151),
            rec(["1.0+3", "1.0+5", "2", "2", "0", "1"], 1, 2, 151),
            rec(["0.0", "0.948", "1", "0", "1", "0"], 1, 2, 151),
            rec(["2.360060+2", "0.0", "0", "0", "1", "0"], 1, 2, 151),
            rec(["0.5", "0.0", "2", "0", "18", "2"], 1, 2, 151),
            rec(["0.0", "0.0", "0.0", "1.0", "0.0", "0.0"], 1, 2, 151),
            rec(["1.0+3", "20.0", "0.0", "1.0-3", "2.3-2", "0.0"], 1, 2, 151),
            rec(["1.0+5", "20.0", "0.0", "1.0-3", "2.3-2", "0.0"], 1, 2, 151),
            rec(["", "", "", "", "", ""], 1, 2, 0),
            mf3(1, "11.0"),
            mf3(2, "10.0"),
            mf3(102, "1.0"),
        ]
        .concat();
        Tape::parse(&text).unwrap()
    }

    /// 分離共鳴（MLBW、1e-5 eV〜20 keV、6.67 eV に1つ）を持つテープ
    fn resolved_tape() -> Tape {
        use crate::endf::rec;
        let text = [
            rec(["9.223800+4", "2.360060+2", "0", "0", "1", "0"], 1, 2, 151),
            rec(["9.223800+4", "1.0", "0", "0", "1", "0"], 1, 2, 151),
            rec(["1.0-5", "2.0+4", "1", "2", "0", "0"], 1, 2, 151),
            rec(["0.0", "0.948", "0", "0", "1", "0"], 1, 2, 151),
            rec(["2.360060+2", "0.0", "0", "0", "6", "1"], 1, 2, 151),
            rec(
                ["6.674", "0.5", "2.6e-2", "1.5e-3", "2.3e-2", "0.0"],
                1,
                2,
                151,
            ),
            rec(["", "", "", "", "", ""], 1, 2, 0),
        ]
        .concat();
        Tape::parse(&text).unwrap()
    }

    #[test]
    fn slowing_down_covers_the_resolved_range_in_whole_groups() {
        let groups = GroupStructure::new(vec![1.0e5, 3.0e4, 1.0e4, 10.0, 1.0e-5]).unwrap();
        let slowing = slowing_down_range(&resolved_tape(), 1, &groups, &Weight::lwr())
            .unwrap()
            .unwrap();
        assert_eq!(slowing.awr, 236.006);
        // 重み関数の熱中性子の区切り 0.1 eV 以上で最も低い群の境界から、分離共鳴の上端 20 keV 以下で
        // 最も高い群の境界まで
        assert_eq!(slowing.range, (10.0, 1.0e4));
    }

    #[test]
    fn potential_scattering_comes_from_the_s_wave_radius() {
        use crate::endf::rec;
        let expected = 4.0 * std::f64::consts::PI * 0.948 * 0.948;
        // MLBW: AP
        let mlbw = potential_scattering(&resolved_tape(), 1).unwrap().unwrap();
        assert!((mlbw - expected).abs() < 1e-12 * expected);
        // Reich-Moore で AP が 0 なら、l=0 の APL（U-238 と同じ）
        let text = [
            rec(["9.223800+4", "2.360060+2", "0", "0", "1", "0"], 1, 2, 151),
            rec(["9.223800+4", "1.0", "0", "0", "1", "0"], 1, 2, 151),
            rec(["1.0-5", "2.0+4", "1", "3", "0", "0"], 1, 2, 151),
            rec(["0.0", "0.0", "0", "0", "1", "0"], 1, 2, 151),
            rec(["2.360060+2", "0.948", "0", "0", "6", "1"], 1, 2, 151),
            rec(
                ["6.674", "0.5", "1.5e-3", "2.3e-2", "0.0", "0.0"],
                1,
                2,
                151,
            ),
            rec(["", "", "", "", "", ""], 1, 2, 0),
        ]
        .concat();
        let rm = potential_scattering(&Tape::parse(&text).unwrap(), 1)
            .unwrap()
            .unwrap();
        assert!((rm - expected).abs() < 1e-12 * expected);
        // 分離共鳴がなければ持たない
        assert!(
            potential_scattering(&unresolved_tape(), 1)
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn no_resolved_range_needs_no_slowing_down() {
        let groups = GroupStructure::new(vec![1.0e5, 1.0e4, 1.0e-5]).unwrap();
        assert!(
            slowing_down_range(&unresolved_tape(), 1, &groups, &Weight::lwr())
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn unresolved_shielding_lowers_capture_at_small_background() {
        use crate::multigroup::shielded_averages;
        let tape = unresolved_tape();
        let sigma0 = [1.0e10, 1.0];
        let shielding = unresolved_shielding(&tape, 1, 300.0, &sigma0)
            .unwrap()
            .unwrap();
        let tables = crate::pointwise::reconstruct(&tape, 1, 1.0e-3).unwrap();
        let groups = GroupStructure::new(vec![5.0e4, 1.0e4]).unwrap();
        let shielded = shielded_averages(
            &groups,
            &Weight::lwr(),
            &tables,
            None,
            &sigma0,
            Some(&shielding),
            None,
        )
        .unwrap();
        let plain =
            shielded_averages(&groups, &Weight::lwr(), &tables, None, &sigma0, None, None).unwrap();
        // 無限希釈では遮蔽因子は 1、σ0 = 1 b では捕獲が減る
        let capture = |avg: &crate::multigroup::ReactionAverages| avg.sigma[&102][0];
        assert!((capture(&shielded[0]) / capture(&plain[0]) - 1.0).abs() < 1e-5);
        assert!(capture(&shielded[1]) < 0.99 * capture(&plain[1]));
    }

    #[test]
    fn no_unresolved_range_gives_no_shielding() {
        use crate::endf::rec;
        let text = [
            rec(["1.001000+3", "1.0", "0", "0", "0", "0"], 1, 3, 1),
            rec(["0.0", "0.0", "0", "0", "1", "2"], 1, 3, 1),
            rec(["2", "2", "", "", "", ""], 1, 3, 1),
            rec(["1.0-5", "1.0", "2.0+7", "1.0", "", ""], 1, 3, 1),
            rec(["", "", "", "", "", ""], 1, 3, 0),
        ]
        .concat();
        let tape = Tape::parse(&text).unwrap();
        assert_eq!(
            unresolved_shielding(&tape, 1, 300.0, &[10.0]).unwrap(),
            None
        );
    }

    /// MF7 の MT4（定数と β が1点の小さな表）と、指定があれば MT2 を持つ熱中性子散乱のテープ。
    fn thermal_tape(with_elastic: bool) -> Tape {
        use crate::endf::rec;
        let mut lines: Vec<String> = vec![
            rec(["1.001000+3", "9.991673-1", "0", "0", "0", "0"], 1, 7, 4),
            rec(["0.0", "0.0", "0", "0", "6", "0"], 1, 7, 4),
            rec(["4.0", "10.0", "9.9917-1", "5.0", "0.0", "2.0"], 1, 7, 4),
            rec(["0.0", "0.0", "0", "0", "1", "2"], 1, 7, 4),
            rec(["2", "4", "", "", "", ""], 1, 7, 4),
        ];
        for beta in ["0.0", "2.0"] {
            lines.push(rec(["2.936000+2", beta, "0", "0", "1", "2"], 1, 7, 4));
            lines.push(rec(["2", "4", "", "", "", ""], 1, 7, 4));
            lines.push(rec(["0.1", "1.0", "10.0", "0.01", "", ""], 1, 7, 4));
        }
        lines.extend([
            rec(["0.0", "0.0", "0", "0", "1", "1"], 1, 7, 4),
            rec(["1", "2", "", "", "", ""], 1, 7, 4),
            rec(["2.936000+2", "1400.0", "", "", "", ""], 1, 7, 4),
            rec(["", "", "", "", "", ""], 1, 7, 0),
        ]);
        if with_elastic {
            // 非干渉性弾性散乱（LTHR=2）: SB と Debye-Waller 積分
            lines.splice(
                0..0,
                [
                    rec(["1.001000+3", "9.991673-1", "2", "0", "0", "0"], 1, 7, 2),
                    rec(["2.0", "0.0", "0", "0", "1", "1"], 1, 7, 2),
                    rec(["1", "2", "", "", "", ""], 1, 7, 2),
                    rec(["2.936000+2", "10.0", "", "", "", ""], 1, 7, 2),
                    rec(["", "", "", "", "", ""], 1, 7, 0),
                ],
            );
        }
        Tape::parse(&lines.concat()).unwrap()
    }

    #[test]
    fn reads_incoherent_inelastic_scattering() {
        let data = read_thermal_scattering(&thermal_tape(false)).unwrap();
        assert_eq!(data.beta, vec![0.0, 2.0]);
    }

    #[test]
    fn thermal_elastic_scattering_is_rejected() {
        assert!(read_thermal_scattering(&thermal_tape(true)).is_err());
    }

    #[test]
    fn processes_each_temperature_and_sigma0() {
        let tape = unresolved_tape();
        let settings = Settings {
            groups: GroupStructure::wims69(),
            weight: Weight::lwr(),
            sigma0: vec![1.0e10, 1.0],
            thermal_cutoff: 4.0,
            tolerance: 1.0e-3,
        };
        let data = process_nuclide(&tape, "X", &[293.6, 600.0], None, &settings).unwrap();
        assert_eq!(data.name, "X");
        let temperatures: Vec<f64> = data.temperatures.iter().map(|t| t.temperature).collect();
        assert_eq!(temperatures, vec![293.6, 600.0]);
        for t in &data.temperatures {
            assert_eq!(t.shielded.len(), 2);
            let [dilute, shielded] = &t.shielded[..] else {
                panic!("two σ0 values")
            };
            let last = dilute.total.len() - 1;
            // 熱群の補正は全断面積と弾性散乱の両方に入るので、差は同じ温度に拡げた捕獲の群平均に等しい
            // （一定の 1 b でも、拡げると低いエネルギーで大きくなる）
            let tables = crate::pointwise::reconstruct_at(&tape, 1, 1.0e-3, t.temperature).unwrap();
            let capture = crate::multigroup::reaction_averages(
                &settings.groups,
                &settings.weight,
                &tables,
                None,
            )
            .sigma[&102][last];
            assert!(capture > 1.0);
            let absorption = dilute.total[last] - dilute.elastic[last];
            // σ0 = 1e10 の重み 1/(σt + σ0) と無限希釈の重みは、σt/σ0 ≈ 1e-9 だけ違う
            assert!(
                (absorption - capture).abs() <= 1e-9 * capture,
                "{absorption} vs {capture}"
            );
            // 自由ガスの行の和は温度に拡げた弾性散乱の群平均になり、熱領域の補正は 0 になる
            let elastic = crate::multigroup::reaction_averages(
                &settings.groups,
                &settings.weight,
                &tables,
                None,
            )
            .sigma[&2][last];
            assert!(
                (dilute.elastic[last] - elastic).abs() <= 1e-9 * elastic,
                "{} vs {elastic}",
                dilute.elastic[last]
            );
            // 一定の 10 b を拡げると 10 N(E)（N(E) = (1 + 1/(2a²)) erf(a) + e^(-a²)/(a√π)、a² = A E/kT）。
            // 許容誤差はポイントワイズの断面積の許容誤差 1e-3
            let (awr, kt) = (236.006, BK * t.temperature);
            let effective = |e: f64| {
                let a = (awr * e / kt).sqrt();
                let n = (1.0 + 0.5 / (a * a)) * crate::special::erf(a)
                    + (-a * a).exp() / (a * std::f64::consts::PI.sqrt());
                10.0 * n
            };
            let integrate = |f: &dyn Fn(f64) -> f64| {
                crate::multigroup::group_integrals(&settings.groups, &settings.weight, &[], f)[last]
            };
            let expected = integrate(&effective) / integrate(&|_| 1.0);
            assert!(
                (elastic - expected).abs() <= 1e-3 * expected,
                "{elastic} vs {expected}"
            );
            // 非分離共鳴の群（[1.5e4, 2.5e4] eV）では、σ0 が小さいと吸収が自己遮蔽で減る
            let g = settings.groups.find(2.0e4).unwrap();
            let absorption = |s: &Shielded| s.total[g] - s.elastic[g];
            assert!(absorption(shielded) < absorption(dilute));
            assert!(t.chi.iter().all(|&c| c == 0.0));
        }
    }
}
