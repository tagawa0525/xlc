//! ENDF の1核種から多群定数を作る処理の全体。

use std::collections::BTreeMap;

use crate::endf::{Mf3, Mf5, Mf6, Tape, parse_nubar};
use crate::library::Nuclide;
use crate::multigroup::{
    FreeGas, GroupStructure, Weight, emits_neutrons, fission_spectrum, reaction_averages,
    scattering_matrix,
};
use crate::physics::BK;
use crate::pointwise::reconstruct;

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
        Some(&thermal),
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
}
