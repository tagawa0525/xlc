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
    _solution: &Solution,
    _materials: &[&Material],
    _names: &[String],
    _physical: &[Vec<f64>],
) -> Result<BTreeMap<String, Rates>, BalanceError> {
    Err(BalanceError::NoProduction)
}
