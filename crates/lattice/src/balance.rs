//! 外周が鏡面反射の格子の中性子の釣り合い。物質ごと・群ごとの物理的な中性子束、正味の吸収、生成と、群の間の散乱による
//! 移行の反応率を、生成の和が 1 になるように規格化して求める（参照解と比べる用）。

use std::collections::BTreeMap;

use crate::material::Material;
use crate::moc::Solution;

/// 1 つの物質（名前）の、群ごとの値（群は高いエネルギーから）。
#[derive(Debug, Clone, PartialEq)]
pub struct Rates {
    /// 物理的な中性子束 × 体積
    pub flux: Vec<f64>,
    /// 正味の吸収（全断面積 − 群から出る散乱の行の和。(n,2n) などで増える中性子を引いた吸収）
    pub absorption: Vec<f64>,
    /// 生成（νΣf）
    pub production: Vec<f64>,
    /// 群の間の散乱による移行 [入射の群][出射の群]（放出の数を含む）
    pub transfer: Vec<Vec<f64>>,
}

/// 釣り合いを求めるときのエラー。
#[derive(Debug, thiserror::Error)]
pub enum BalanceError {
    #[error("{given} values for {regions} regions")]
    Regions { given: usize, regions: usize },
    #[error("region {0}: the flux and the material differ in groups")]
    Groups(usize),
    #[error("no production to normalize by")]
    NoProduction,
}

/// 領域ごとの物質 `materials`（名前は `names`）、MOC の解 `solution` と、物理的な中性子束 `physical`（[領域][群]）から、
/// 物質ごとの釣り合いを求める。真空の領域は数えない。
pub fn by_material(
    solution: &Solution,
    materials: &[&Material],
    names: &[String],
    physical: &[Vec<f64>],
) -> Result<BTreeMap<String, Rates>, BalanceError> {
    let regions = solution.flux.len();
    for given in [materials.len(), names.len(), physical.len()] {
        if given != regions {
            return Err(BalanceError::Regions { given, regions });
        }
    }
    let mut by_name: BTreeMap<String, Rates> = BTreeMap::new();
    for (region, (name, m)) in names.iter().zip(materials).enumerate() {
        if m.is_void() {
            continue;
        }
        let groups = m.groups();
        if solution.flux[region].len() != groups || physical[region].len() != groups {
            return Err(BalanceError::Groups(region));
        }
        let rates = by_name.entry(name.clone()).or_insert_with(|| Rates {
            flux: vec![0.0; groups],
            absorption: vec![0.0; groups],
            production: vec![0.0; groups],
            transfer: vec![vec![0.0; groups]; groups],
        });
        let volume = solution.volumes[region];
        for (g, phi_physical) in physical[region].iter().enumerate() {
            // 反応率は MOC の中性子束と領域の物質（SPH 因子を掛けたもの）から求め、中性子束だけを物理的な値にする
            let phi = solution.flux[region][g] * volume;
            let out: f64 = m.scatter[g].iter().sum();
            rates.flux[g] += phi_physical * volume;
            rates.absorption[g] += (m.total[g] - out) * phi;
            rates.production[g] += m.nu_fission[g] * phi;
            for (to, s) in m.scatter[g].iter().enumerate() {
                rates.transfer[g][to] += s * phi;
            }
        }
    }
    let production: f64 = by_name.values().flat_map(|r| &r.production).sum();
    if production.is_nan() || production <= 0.0 {
        return Err(BalanceError::NoProduction);
    }
    for rates in by_name.values_mut() {
        for v in rates
            .flux
            .iter_mut()
            .chain(&mut rates.absorption)
            .chain(&mut rates.production)
            .chain(rates.transfer.iter_mut().flatten())
        {
            *v /= production;
        }
    }
    Ok(by_name)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Lattice, PinCell};
    use crate::material::FissionTerm;
    use crate::moc::{Options, solve};
    use crate::tracking::Tracking;

    fn two_group(total: [f64; 2], scatter: [[f64; 2]; 2], nu_fission: [f64; 2]) -> Material {
        Material {
            total: total.to_vec(),
            scatter: scatter.iter().map(|r| r.to_vec()).collect(),
            anisotropic: Vec::new(),
            fission: nu_fission.iter().map(|v| v / 2.5).collect(),
            nu_fission: nu_fission.to_vec(),
            fission_terms: if nu_fission.iter().any(|&v| v > 0.0) {
                vec![FissionTerm {
                    chi: vec![1.0, 0.0],
                    nu_fission: nu_fission.to_vec(),
                }]
            } else {
                vec![]
            },
        }
    }

    #[test]
    fn the_balance_of_a_reflective_cell_closes() {
        let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
        // (n,2n) のように、散乱の行の和が全断面積を超える群のある減速材
        let moderator = two_group([0.6, 2.0], [[0.52, 0.1], [0.0, 1.98]], [0.0, 0.0]);
        let cell = PinCell::new(1.26, vec![0.25, 0.4096, 0.55]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 16, 0.05).unwrap();
        let materials = [&fuel, &fuel, &moderator, &moderator];
        let names: Vec<String> = ["fuel", "fuel", "water", "water"]
            .into_iter()
            .map(String::from)
            .collect();
        let solution = solve(&tracking, &materials, &Options::default()).unwrap();
        // 物理的な中性子束は MOC の値の 2 倍とする（中性子束だけが変わり、反応率は変わらない）
        let physical: Vec<Vec<f64>> = solution
            .flux
            .iter()
            .map(|f| f.iter().map(|v| 2.0 * v).collect())
            .collect();
        let rates = by_material(&solution, &materials, &names, &physical).unwrap();
        // 生成の和は 1、正味の吸収の和は 1/k（外周が鏡面反射）
        let production: f64 = rates.values().flat_map(|r| &r.production).sum();
        let absorption: f64 = rates.values().flat_map(|r| &r.absorption).sum();
        assert!((production - 1.0).abs() < 1e-12, "{production}");
        assert!(
            (absorption * solution.k - 1.0).abs() < 1e-6,
            "{absorption} {}",
            solution.k
        );
        // 群から出る移行の和と正味の吸収の和は、全断面積の反応率
        let water = &rates["water"];
        let phi_moc: f64 = water.flux[0] / 2.0;
        let out: f64 = water.transfer[0].iter().sum();
        assert!((out + water.absorption[0] - 0.6 * phi_moc).abs() < 1e-12);
        // 中性子束だけが物理的な値（MOC の 2 倍）になる
        let fuel_flux: f64 = (0..2)
            .map(|r| solution.flux[r][1] * solution.volumes[r])
            .sum::<f64>();
        assert!(
            (rates["fuel"].flux[1] - 2.0 * fuel_flux / production_of(&solution, &materials)).abs()
                < 1e-12
        );
    }

    #[test]
    fn mismatched_inputs_are_errors() {
        let fuel = two_group([0.5, 1.2], [[0.45, 0.02], [0.0, 0.9]], [0.008, 0.25]);
        let cell = PinCell::new(1.26, vec![0.4096]).unwrap();
        let tracking = Tracking::new(&Lattice::pin_cell(cell), 8, 0.1).unwrap();
        let materials = [&fuel, &fuel];
        let names: Vec<String> = vec!["fuel".into(); 2];
        let solution = solve(&tracking, &materials, &Options::default()).unwrap();
        assert!(matches!(
            by_material(&solution, &materials[..1], &names, &solution.flux),
            Err(BalanceError::Regions { .. })
        ));
        let short: Vec<Vec<f64>> = solution.flux.iter().map(|f| f[..1].to_vec()).collect();
        assert!(matches!(
            by_material(&solution, &materials, &names, &short),
            Err(BalanceError::Groups(0))
        ));
    }

    /// 規格化の前の生成の和
    fn production_of(solution: &Solution, materials: &[&Material]) -> f64 {
        materials
            .iter()
            .enumerate()
            .map(|(r, m)| {
                m.nu_fission
                    .iter()
                    .zip(&solution.flux[r])
                    .map(|(n, f)| n * f)
                    .sum::<f64>()
                    * solution.volumes[r]
            })
            .sum()
    }
}
