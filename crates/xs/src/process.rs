//! ENDF の1核種から多群定数を作る処理の全体。

use std::collections::BTreeMap;

use crate::endf::{Mf2, Mf3, Mf5, Mf6, RangeKind, Tab1, Tape, parse_nubar};
use crate::library::Nuclide;
use crate::multigroup::{
    FreeGas, GroupStructure, Thermal, UnresolvedShielding, Weight, emits_neutrons,
    fission_spectrum, reaction_averages, scattering_matrix,
};
use crate::physics::BK;
use crate::pointwise::{DEFAULT_TOLERANCE, left_limit, linearize, reconstruct, right_limit};
use crate::resonance::{Sigma, UnresolvedXs, unresolved_grid, unresolved_interpolation};

/// 処理の条件。
#[derive(Debug, Clone)]
pub struct Settings {
    pub groups: GroupStructure,
    pub weight: Weight,
    /// 自由ガス模型の温度 (K)
    pub temperature: f64,
    /// この入射エネルギー (eV) 未満の弾性散乱を自由ガス模型で扱う
    pub thermal_cutoff: f64,
    /// ポイントワイズの断面積の許容誤差
    pub tolerance: f64,
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

/// 熱領域で弾性散乱の実効断面積が 0 K の値と違う分だけ、全断面積を補正する。
///
/// 自由ガスの散乱行列の行の和は、標的の熱運動を含む実効的な弾性散乱断面積になる。
/// 全断面積をそれに合わせて、全断面積 = 吸収 + 散乱 の関係を保つ。
pub fn effective_total(total: &[f64], elastic: &[f64], elastic_matrix: &[Vec<f64>]) -> Vec<f64> {
    total
        .iter()
        .zip(elastic)
        .zip(elastic_matrix)
        .map(|((t, el), row)| t - el + row.iter().sum::<f64>())
        .collect()
}

/// テープ中の（ただ1つの）物質を処理する。
pub fn process_nuclide(tape: &Tape, name: &str, settings: &Settings) -> Result<Nuclide, String> {
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
    let (groups, weight) = (&settings.groups, &settings.weight);
    let ng = groups.len();

    let tables = reconstruct(tape, mat, settings.tolerance).map_err(text)?;
    let nubar = match tape.section(mat, 1, 452) {
        Some(section) => Some(parse_nubar(&section).map_err(text)?),
        None => None,
    };
    let averages = reaction_averages(groups, weight, &tables, nubar.as_ref());
    let fissile = tables.contains_key(&18);
    let chi = match (fissile, nubar.as_ref(), tape.section(mat, 5, 18)) {
        (false, ..) => vec![0.0; ng],
        (true, Some(nu), Some(section)) => {
            let mf5 = Mf5::parse(&section).map_err(text)?;
            fission_spectrum(groups, weight, &tables, nu, &mf5)
        }
        (true, ..) => return Err("fission needs MF1 MT452 and MF5 MT18".into()),
    };

    let mut q_values = BTreeMap::new();
    for &mt in tables.keys().filter(|mt| (51..=90).contains(*mt)) {
        q_values.insert(mt, mf3(mt)?.qi);
    }
    let mut mf6 = BTreeMap::new();
    for mt in tape.mts(mat, 6) {
        if mt != 2 && emits_neutrons(mt) && tables.contains_key(&mt) {
            let section = tape.section(mat, 6, mt).expect("listed by mts");
            mf6.insert(mt, Mf6::parse(&section).map_err(text)?);
        }
    }

    // 弾性散乱だけ自由ガスを使い、その実効断面積で全断面積を補正する
    let thermal = FreeGas {
        kt: BK * settings.temperature,
        cutoff: settings.thermal_cutoff,
    };
    let elastic_table = tables.get(&2).ok_or("MF3 MT2 not found")?;
    let elastic_only = BTreeMap::from([(2, elastic_table.clone())]);
    let others: BTreeMap<_, _> = tables
        .iter()
        .filter(|(mt, _)| **mt != 2)
        .map(|(mt, t)| (*mt, t.clone()))
        .collect();
    let awr = head.awr;
    let elastic = scattering_matrix(
        groups,
        weight,
        awr,
        &elastic_only,
        &q_values,
        &mf6,
        Some(&Thermal::FreeGas(thermal)),
    )?;
    let mut scatter = scattering_matrix(groups, weight, awr, &others, &q_values, &mf6, None)?;
    for (row, el) in scatter.iter_mut().zip(&elastic) {
        for (v, e) in row.iter_mut().zip(el) {
            *v += e;
        }
    }
    let sigma = |mt: u32| averages.sigma.get(&mt).cloned().unwrap_or(vec![0.0; ng]);
    Ok(Nuclide {
        name: name.to_string(),
        za: head.za,
        awr,
        temperature: settings.temperature,
        total: effective_total(&sigma(1), &sigma(2), &elastic),
        fission: sigma(18),
        nu_fission: averages.nu_fission.unwrap_or(vec![0.0; ng]),
        chi,
        scatter,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn total_follows_effective_elastic() {
        let total = vec![10.0, 30.0];
        let elastic = vec![4.0, 20.0];
        let matrix = vec![vec![3.9, 0.1], vec![2.0, 27.0]];
        assert_eq!(effective_total(&total, &elastic, &matrix), vec![10.0, 39.0]);
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
}
