//! 外周が鏡面反射の格子の中性子収支。物質ごと・群ごとの物理的な中性子束、正味の吸収、生成と、群の間の散乱による
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

/// バランスを求めるときのエラー。
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
/// 物質ごとのバランスを求める。真空の領域は数えない。
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
mod tests;
