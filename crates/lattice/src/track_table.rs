//! 円筒のセルの飛程の表（ADR-0156）。
//!
//! 同じ物質のリングは同じ全断面積を持つので、セルの飛程 W は物質ごとの全断面積 Σ_m の関数になる。物質ごとに ln Σ の
//! 等間隔の格子をとり、格子点の組で W を求めておき、4 点の Lagrange 補間（物質ごとの積）で内挿する。
//!
//! 内挿するのは R_ij = W_ij Σ_j（真空のリング j は W_ij のまま）。真空でない列についての R の行の和は体積 V_i
//! （白色境界ではどこかで必ず衝突する）で、補間の重みの和は 1 なので、内挿した R も中性子の保存を丸めの範囲で満たす。
//! W は R_ij/Σ_j で戻す。

use rayon::prelude::*;

use crate::cp::{CpError, Cylinder};

/// 表の細かさの既定（ln Σ の 10 倍あたりの格子点の数。ADR-0156）
pub const DEFAULT_DENSITY: usize = 16;

/// 1つの物質の内挿に使う格子点の数（4 点の Lagrange 補間）
const STENCIL: usize = 4;

/// 表を作るときのエラー。
#[derive(Debug, thiserror::Error)]
pub enum TableError {
    #[error("the track table needs at least one point per decade")]
    Density,
    #[error("material {material}: total cross sections {range:?} must be positive and finite")]
    Range { material: usize, range: [f64; 2] },
    #[error("{given} ring materials for {rings} rings")]
    Rings { given: usize, rings: usize },
    #[error("the number of table nodes overflows")]
    Size,
    #[error("{given} total cross sections for {materials} materials")]
    Materials { given: usize, materials: usize },
    #[error("material {material}: total cross section {sigma} is outside the table {range:?}")]
    Outside {
        material: usize,
        sigma: f64,
        range: [f64; 2],
    },
    #[error(transparent)]
    Cp(#[from] CpError),
}

/// 1つの物質の ln Σ の等間隔の格子
struct Axis {
    start: f64,
    step: f64,
    count: usize,
    /// 覆う全断面積の [最小, 最大]
    range: [f64; 2],
}

impl Axis {
    /// [min, max] を覆う格子。幅が 0 なら 1 点、それ以外は 10 倍あたり `density` 点以上で、最低 [`STENCIL`] 点
    fn new(min: f64, max: f64, density: usize) -> Self {
        let (start, span) = (min.ln(), max.ln() - min.ln());
        let range = [min, max];
        if span == 0.0 {
            return Self {
                start,
                step: 0.0,
                count: 1,
                range,
            };
        }
        let step = std::f64::consts::LN_10 / density as f64;
        // 細かさが大きすぎると点の数が usize を超える。そのときは usize::MAX にして、組の数の検査で止める
        let intervals = (span / step).ceil();
        let count = if intervals < usize::MAX as f64 {
            (intervals as usize).saturating_add(1).max(STENCIL)
        } else {
            usize::MAX
        };
        Self {
            start,
            step: span / (count - 1) as f64,
            count,
            range,
        }
    }

    fn sigma(&self, k: usize) -> f64 {
        (self.start + k as f64 * self.step).exp()
    }

    /// `sigma` の内挿に使う最初の格子点と、点ごとの重み（使う点の数は `count` と [`STENCIL`] の小さいほう）
    fn stencil(&self, sigma: f64) -> (usize, [f64; STENCIL]) {
        if self.count == 1 {
            return (0, [1.0, 0.0, 0.0, 0.0]);
        }
        let s = (sigma.ln() - self.start) / self.step;
        // 点 s を囲む 4 点（端では内側に寄せる）
        let last = (self.count - STENCIL) as f64;
        let first = (s.floor() - 1.0).clamp(0.0, last) as usize;
        let x = s - first as f64;
        // 点 0, 1, 2, 3 の Lagrange の基底
        let weights = [
            -(x - 1.0) * (x - 2.0) * (x - 3.0) / 6.0,
            x * (x - 2.0) * (x - 3.0) / 2.0,
            -x * (x - 1.0) * (x - 3.0) / 2.0,
            x * (x - 1.0) * (x - 2.0) / 6.0,
        ];
        (first, weights)
    }

    fn width(&self) -> usize {
        self.count.min(STENCIL)
    }
}

/// 物質ごとの全断面積の格子での飛程の表。
pub struct TrackTable {
    /// [リング] の物質の番号（真空は `None`）
    classes: Vec<Option<usize>>,
    axes: Vec<Axis>,
    /// [格子点の組 × リング × リング] の R_ij（格子点の組は最初の物質がいちばん速く変わる）
    values: Vec<f64>,
}

impl TrackTable {
    /// `classes` はリングごとの物質の番号（真空は `None`）、`ranges` は物質ごとの全断面積の [最小, 最大] (1/cm)、
    /// `density` は ln Σ の 10 倍あたりの格子点の数。
    pub fn new(
        cell: &Cylinder,
        classes: Vec<Option<usize>>,
        ranges: &[[f64; 2]],
        density: usize,
    ) -> Result<Self, TableError> {
        if density == 0 {
            return Err(TableError::Density);
        }
        let rings = cell.volumes().len();
        if classes.len() != rings {
            return Err(TableError::Rings {
                given: classes.len(),
                rings,
            });
        }
        let axes = axes(ranges, density)?;
        let nodes = count(&axes)?;
        let size = nodes.checked_mul(rings * rings).ok_or(TableError::Size)?;
        // 格子点の組ごとの R をその場の並びに書く（組ごとの小さな確保と、集めてからつなぐコピーを避ける）
        let mut values = vec![0.0; size];
        values
            .par_chunks_mut(rings * rings)
            .enumerate()
            .try_for_each(|(node, r)| {
                let mut rest = node;
                let sigma: Vec<f64> = axes
                    .iter()
                    .map(|a| {
                        let k = rest % a.count;
                        rest /= a.count;
                        a.sigma(k)
                    })
                    .collect();
                let totals = totals(&classes, &sigma);
                let w = cell.tracks(&totals)?;
                for (row, out) in w.iter().zip(r.chunks_mut(rings)) {
                    for ((x, w), &t) in out.iter_mut().zip(row).zip(&totals) {
                        *x = if t > 0.0 { w * t } else { *w };
                    }
                }
                Ok::<(), CpError>(())
            })?;
        tracing::debug!(nodes, "track table");
        Ok(Self {
            classes,
            axes,
            values,
        })
    }

    /// 格子点の組の数
    pub fn nodes(&self) -> usize {
        self.axes.iter().map(|a| a.count).product()
    }

    /// `ranges` と `density` で作る表の格子点の組の数（作る前に、表にする価値があるかを見る）
    pub fn count(ranges: &[[f64; 2]], density: usize) -> Result<usize, TableError> {
        count(&axes(ranges, density)?)
    }

    /// 物質ごとの全断面積 `sigma` での飛程 T[j][i]（リング i の単位の線源密度あたりの、リング j の中を飛ぶ長さ）。
    /// 表の範囲の外（外挿になる）と、物質の数の違う `sigma` はエラー。
    pub fn tracks(&self, sigma: &[f64]) -> Result<Vec<Vec<f64>>, TableError> {
        if sigma.len() != self.axes.len() {
            return Err(TableError::Materials {
                given: sigma.len(),
                materials: self.axes.len(),
            });
        }
        for (material, (axis, &s)) in self.axes.iter().zip(sigma).enumerate() {
            if !(axis.range[0] <= s && s <= axis.range[1]) {
                return Err(TableError::Outside {
                    material,
                    sigma: s,
                    range: axis.range,
                });
            }
        }
        let n = self.classes.len();
        let stencils: Vec<(usize, [f64; STENCIL])> = self
            .axes
            .iter()
            .zip(sigma)
            .map(|(a, &s)| a.stencil(s))
            .collect();
        let mut r = vec![0.0; n * n];
        // 物質ごとの内挿の点の組を、最初の物質がいちばん速く変わる順に回る
        let mut offsets = vec![0usize; self.axes.len()];
        'combos: loop {
            let (mut node, mut stride, mut weight) = (0, 1, 1.0);
            for ((axis, &(first, w)), &k) in self.axes.iter().zip(&stencils).zip(&offsets) {
                node += (first + k) * stride;
                stride *= axis.count;
                weight *= w[k];
            }
            for (x, v) in r
                .iter_mut()
                .zip(&self.values[node * n * n..(node + 1) * n * n])
            {
                *x += weight * v;
            }
            for (k, axis) in offsets.iter_mut().zip(&self.axes) {
                *k += 1;
                if *k < axis.width() {
                    continue 'combos;
                }
                *k = 0;
            }
            break;
        }
        let totals = totals(&self.classes, sigma);
        Ok((0..n)
            .map(|j| {
                let t = totals[j];
                (0..n)
                    .map(|i| {
                        if t > 0.0 {
                            r[i * n + j] / t
                        } else {
                            r[i * n + j]
                        }
                    })
                    .collect()
            })
            .collect())
    }
}

/// 物質ごとの ln Σ の格子
fn axes(ranges: &[[f64; 2]], density: usize) -> Result<Vec<Axis>, TableError> {
    if density == 0 {
        return Err(TableError::Density);
    }
    ranges
        .iter()
        .enumerate()
        .map(|(material, &[min, max])| {
            if !(min > 0.0 && min <= max && max.is_finite()) {
                return Err(TableError::Range {
                    material,
                    range: [min, max],
                });
            }
            Ok(Axis::new(min, max, density))
        })
        .collect()
}

/// 格子点の組の数（usize に収まらなければエラー）
fn count(axes: &[Axis]) -> Result<usize, TableError> {
    axes.iter()
        .try_fold(1usize, |n, a| n.checked_mul(a.count))
        .ok_or(TableError::Size)
}

/// [リング] の全断面積（真空は 0）
fn totals(classes: &[Option<usize>], sigma: &[f64]) -> Vec<f64> {
    classes
        .iter()
        .map(|m| m.map_or(0.0, |m| sigma[m]))
        .collect()
}

#[cfg(test)]
mod tests;
