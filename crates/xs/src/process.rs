//! ENDF の1核種から多群定数を作る処理の全体。

use std::collections::BTreeMap;

use crate::endf::{Mf2, Mf3, Mf5, Mf6, Mf7Incoherent, RangeKind, Tab1, Tape, parse_nubar};
use crate::library::{NuclideData, Shielded, TemperatureData};
use crate::multigroup::{
    Emission, FreeGas, GroupStructure, ReactionAverages, Thermal, ThermalSab, UnresolvedShielding,
    Weight, emits_neutrons, fission_spectrum, scattering_matrix, shielded_averages,
};
use crate::physics::BK;
use crate::pointwise::{ColdParts, DEFAULT_TOLERANCE, left_limit, linearize, right_limit};
use crate::resonance::{Sigma, UnresolvedXs, unresolved_grid, unresolved_interpolation};

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

/// 熱中性子散乱のテープ（MF7）から、非干渉性非弾性散乱の S(α, β)（MT4）を読む。
///
/// 弾性散乱（MT2）を持つデータはエラーにする。S(α, β) の散乱行列は cutoff 未満で MF3 の弾性散乱を
/// 置き換えるので、MT2 を扱わないと弾性散乱が黙って抜け落ちる（黒鉛、ZrH など）。
pub fn read_thermal_scattering(tape: &Tape) -> Result<Mf7Incoherent, String> {
    let mat = *tape
        .materials()
        .first()
        .ok_or("empty thermal scattering tape")?;
    if tape.section(mat, 7, 2).is_some() {
        return Err("thermal elastic scattering (MF7 MT2) is not supported".into());
    }
    let section = tape
        .section(mat, 7, 4)
        .ok_or("thermal scattering tape has no MF7 MT4")?;
    Mf7Incoherent::parse(&section).map_err(|e| e.to_string())
}

/// 非分離共鳴の自己遮蔽の補正を、温度 `temperature` (K) と背景断面積 `sigma0` について作る
/// （NJOY UNRESR に相当）。非分離共鳴がなければ `None`。
pub fn unresolved_shielding(
    tape: &Tape,
    mat: i32,
    temperature: f64,
    sigma0: &[f64],
) -> Result<Option<UnresolvedShielding>, String> {
    let text = |e: crate::endf::Error| e.to_string();
    let Some(section) = tape.section(mat, 2, 151) else {
        return Ok(None);
    };
    let mf2 = Mf2::parse(&section).map_err(text)?;
    let mut unresolved = None;
    for range in &mf2.ranges {
        match &range.kind {
            RangeKind::Unresolved(urr) if unresolved.is_none() => unresolved = Some((range, urr)),
            RangeKind::Unresolved(_) => return Err("more than one unresolved range".into()),
            RangeKind::UnresolvedLrf1 { .. } => {
                return Err("self-shielding of unresolved LRF=1 ranges is not supported".into());
            }
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
    let background_table = |mt: u32| -> Result<Option<Tab1>, String> {
        if urr.lssf != 0 {
            return Ok(None);
        }
        match tape.section(mat, 3, mt) {
            Some(section) => Ok(Some(linearize(
                &Mf3::parse(&section).map_err(text)?.sigma,
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
        let average = xs.cross_sections(e);
        infinite.push([
            average.total + background.total,
            average.elastic + background.elastic,
            average.fission + background.fission,
            average.capture + background.capture,
        ]);
    }
    UnresolvedShielding::new(range.el, range.eh, urr.lssf, grid, law, shielded, infinite).map(Some)
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

/// 自由ガス模型で散乱行列を作るための弾性散乱断面積。`cutoff` 未満は自由原子の（0 K の）断面積
/// `free_atom`、以上は温度 T に拡げた断面積 `broadened` とする。
///
/// 自由ガスの核は規格化せず、その積分 N(E) が熱運動による実効断面積の増加を表す（詳細バランスを
/// 保つため。`free_gas` を参照）。温度 T に拡げた断面積はすでに N(E) を含むので、核に掛けると
/// N(E) を二重に数える。一定の断面積なら、拡げた断面積はちょうど N(E) 倍になる。
/// NJOY THERMR は代わりに、拡げた断面積に規格化した出射分布を掛ける。
fn free_gas_elastic(free_atom: &Tab1, broadened: &Tab1, cutoff: f64) -> Tab1 {
    let (mut x, mut y) = (Vec::new(), Vec::new());
    for (&e, &v) in free_atom
        .x
        .iter()
        .zip(&free_atom.y)
        .filter(|(e, _)| **e < cutoff)
    {
        x.push(e);
        y.push(v);
    }
    x.extend([cutoff, cutoff]);
    y.extend([
        left_limit(free_atom, cutoff),
        right_limit(broadened, cutoff),
    ]);
    for (&e, &v) in broadened
        .x
        .iter()
        .zip(&broadened.y)
        .filter(|(e, _)| **e > cutoff)
    {
        x.push(e);
        y.push(v);
    }
    Tab1::linear(x, y)
}

/// テープ中の（ただ1つの）物質を、温度 `temperatures` (K) ごとに処理する。
/// `thermal` があれば cutoff 未満の弾性散乱を S(α, β) で、なければ自由ガス模型で扱う。
pub fn process_nuclide(
    tape: &Tape,
    name: &str,
    temperatures: &[f64],
    thermal: Option<&Mf7Incoherent>,
    settings: &Settings,
) -> Result<NuclideData, String> {
    let materials = tape.materials();
    let [mat] = materials[..] else {
        return Err(format!("expected one material, found {materials:?}"));
    };
    let text = |e: crate::endf::Error| e.to_string();
    let mf3 = |mt: u32| -> Result<Mf3, String> {
        let section = tape
            .section(mat, 3, mt)
            .ok_or_else(|| format!("MF3 MT{mt} not found"))?;
        Mf3::parse(&section).map_err(text)
    };
    let first = *tape.mts(mat, 3).first().ok_or("no MF3 data")?;
    let head = mf3(first)?;
    let nubar = match tape.section(mat, 1, 452) {
        Some(section) => Some(parse_nubar(&section).map_err(text)?),
        None => None,
    };
    let cold = ColdParts::read(tape, mat, settings.tolerance).map_err(text)?;
    let mf5 = match tape.section(mat, 5, 18) {
        Some(section) => Some(Mf5::parse(&section).map_err(text)?),
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
    for mt in tape.mts(mat, 6) {
        if mt != 2 && emits_neutrons(mt) {
            let section = tape.section(mat, 6, mt).expect("listed by mts");
            emissions.insert(mt, Emission::Mf6(Mf6::parse(&section).map_err(text)?));
        }
    }
    // MF6 のない反応は、MF4/MF5 で与えられている（MF5 のエネルギー分布だけを使う）
    for mt in tape.mts(mat, 5) {
        if mt != 2 && emits_neutrons(mt) && !emissions.contains_key(&mt) {
            let section = tape.section(mat, 5, mt).expect("listed by mts");
            emissions.insert(mt, Emission::Mf5(Mf5::parse(&section).map_err(text)?));
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
        .map(|&t| common.at_temperature(t).map_err(|e| format!("{t} K: {e}")))
        .collect::<Result<_, _>>()?;
    Ok(NuclideData {
        name: name.to_string(),
        za: head.za,
        awr: head.awr,
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
    fn at_temperature(&self, temperature: f64) -> Result<TemperatureData, String> {
        let settings = self.settings;
        let (groups, weight) = (&settings.groups, &settings.weight);
        let ng = groups.len();
        let tables = self.cold.at(temperature);
        let shielding = unresolved_shielding(self.tape, self.mat, temperature, &settings.sigma0)?;
        let averages = shielded_averages(
            groups,
            weight,
            &tables,
            self.nubar,
            &settings.sigma0,
            shielding.as_ref(),
        )?;
        let chi = match (tables.contains_key(&18), self.nubar, self.mf5) {
            (false, ..) => vec![0.0; ng],
            (true, Some(nu), Some(mf5)) => fission_spectrum(groups, weight, &tables, nu, mf5),
            (true, ..) => return Err("fission needs MF1 MT452 and MF5 MT18".into()),
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
        let elastic_table = tables.get(&2).ok_or("MF3 MT2 not found")?;
        // 自由ガス模型は自由原子の断面積（0 K の弾性散乱）に核を掛ける（`free_gas_elastic` を参照）
        let elastic_table = match thermal {
            Thermal::FreeGas(_) => {
                let free_atom = self.cold.part(2).ok_or("MF3 MT2 not found")?;
                free_gas_elastic(free_atom, elastic_table, settings.thermal_cutoff)
            }
            Thermal::Sab(_) => elastic_table.clone(),
        };
        let elastic_only = BTreeMap::from([(2, elastic_table)]);
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

    #[test]
    fn unresolved_shielding_lowers_capture_at_small_background() {
        use crate::multigroup::shielded_averages;
        let tape = unresolved_tape();
        let sigma0 = [1.0e10, 1.0];
        let shielding = unresolved_shielding(&tape, 1, 300.0, &sigma0)
            .unwrap()
            .unwrap();
        let tables = crate::pointwise::reconstruct(&tape, 1, 1.0e-3).unwrap();
        let groups = GroupStructure::new(vec![5.0e4, 1.0e4]);
        let shielded = shielded_averages(
            &groups,
            &Weight::lwr(),
            &tables,
            None,
            &sigma0,
            Some(&shielding),
        )
        .unwrap();
        let plain =
            shielded_averages(&groups, &Weight::lwr(), &tables, None, &sigma0, None).unwrap();
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
                unreachable!()
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
