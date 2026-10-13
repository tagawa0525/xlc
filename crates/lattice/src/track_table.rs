//! 円筒のセルの飛程の表（ADR-0156）。

use crate::cp::{CpError, Cylinder};

/// 表の細かさの既定（ln Σ の 10 倍あたりの格子点の数。ADR-0156）
pub const DEFAULT_DENSITY: usize = 16;

/// 表を作るときのエラー。
#[derive(Debug, thiserror::Error)]
pub enum TableError {
    #[error("the track table needs at least one point per decade")]
    Density,
    #[error("material {material}: total cross sections {range:?} must be positive and finite")]
    Range { material: usize, range: [f64; 2] },
    #[error("{given} ring materials for {rings} rings")]
    Rings { given: usize, rings: usize },
    #[error(transparent)]
    Cp(#[from] CpError),
}

/// 物質ごとの全断面積の格子での飛程の表。
pub struct TrackTable {}

impl TrackTable {
    /// `classes` はリングごとの物質の番号（真空は `None`）、`ranges` は物質ごとの全断面積の [最小, 最大] (1/cm)、
    /// `density` は ln Σ の 10 倍あたりの格子点の数。
    pub fn new(
        _cell: &Cylinder,
        _classes: Vec<Option<usize>>,
        _ranges: &[[f64; 2]],
        _density: usize,
    ) -> Result<Self, TableError> {
        todo!()
    }

    /// 物質ごとの全断面積 `sigma` での飛程 T[j][i]（リング i の単位の線源密度あたりの、リング j の中を飛ぶ長さ）。
    pub fn tracks(&self, _sigma: &[f64]) -> Vec<Vec<f64>> {
        todo!()
    }
}

#[cfg(test)]
mod tests;
