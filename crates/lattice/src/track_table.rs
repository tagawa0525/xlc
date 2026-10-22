//! 円筒のセルの飛程の表（ADR-0156）。
//!
//! 同じ物質のリングは同じ全断面積を持つので、セルの飛程 W は物質ごとの全断面積 Σ_m の関数になる。物質ごとに ln Σ の
//! 等間隔の格子をとり、格子点の組で W を求めておき、4 点の Lagrange 補間（物質ごとの積）で内挿する。
//!
//! 内挿するのは R_ij = W_ij Σ_j（真空のリング j は W_ij のまま）。真空でない列についての R の行の和は体積 V_i
//! （白色境界ではどこかで必ず衝突する）で、補間の重みの和は 1 なので、内挿した R も中性子の保存を丸めの範囲で満たす。
//! W は R_ij/Σ_j で戻す。
//!
//! 核種の並びが同じで数密度の違う物質（リングごとに燃焼させる燃料）は、1 つの軸を共有する（ADR-0182）。軸の値は、その物質の
//! リングの ln Σ の、体積の重みの平均（中心）。表は物質ごとの ln Σ についての R の微分も持ち、格子点では中心からの差 δ で
//! R + Σ δ ∂R/∂ln Σ を内挿する。δ が [`MAX_DEVIATION`] を超える組は内挿せず、呼び出し側が格子点ごとに求める。

use rayon::prelude::*;

use crate::cp::{CpError, Cylinder};

/// 表の細かさの既定（ln Σ の 10 倍あたりの格子点の数。ADR-0156）
pub const DEFAULT_DENSITY: usize = 16;

/// 1つの物質の内挿に使う格子点の数（4 点の Lagrange 補間）
const STENCIL: usize = 4;

/// 軸を共有する物質の ln Σ の、中心からの差の上限（ADR-0182）。1 次の展開の誤差は差の 2 乗に比例し、この上限で表の内挿の誤差
/// （10 倍あたり 16 点で 1e-5 未満）の水準になる
pub const MAX_DEVIATION: f64 = 0.01;

/// 微分の中心差分の ln Σ の幅（ADR-0182）。打ち切りの誤差は幅の 2 乗で、展開の項（差 ≤ [`MAX_DEVIATION`] を掛ける）では丸めの
/// 水準になる
const STEP: f64 = 1e-3;

/// 表を作るときのエラー。
#[derive(Debug, thiserror::Error)]
pub enum TableError {
    #[error("the track table needs at least one point per decade")]
    Density,
    #[error("material {material}: total cross sections {range:?} must be positive and finite")]
    Range { material: usize, range: [f64; 2] },
    #[error("{given} ring materials for {rings} rings")]
    Rings { given: usize, rings: usize },
    #[error("every material needs an axis, and every axis a material")]
    Axes,
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

/// リングの物質と表の軸の対応。
pub struct Layout {
    /// [リング] の物質の番号（真空は `None`）
    classes: Vec<Option<usize>>,
    /// [物質] の軸の番号
    axes: Vec<usize>,
    /// [物質] の、軸の中心（ln Σ の平均）での体積の重み
    weights: Vec<f64>,
    /// [物質] の微分の番号（軸を共有しない物質は `None`）
    slots: Vec<Option<usize>>,
    /// [軸] の物質の数
    members: Vec<usize>,
}

impl Layout {
    /// `classes` はリングごとの物質の番号（真空は `None`）、`axes` は物質ごとの軸の番号（0 から隙間なく）。
    pub fn new(
        cell: &Cylinder,
        classes: Vec<Option<usize>>,
        axes: Vec<usize>,
    ) -> Result<Self, TableError> {
        let volumes = cell.volumes();
        if classes.len() != volumes.len() {
            return Err(TableError::Rings {
                given: classes.len(),
                rings: volumes.len(),
            });
        }
        let count = axes.iter().max().map_or(0, |&a| a + 1);
        let mut members = vec![0; count];
        for &a in &axes {
            members[a] += 1;
        }
        if members.contains(&0) || classes.iter().flatten().any(|&m| m >= axes.len()) {
            return Err(TableError::Axes);
        }
        let mut material_volume = vec![0.0; axes.len()];
        for (m, v) in classes.iter().zip(&volumes) {
            if let Some(m) = m {
                material_volume[*m] += v;
            }
        }
        let mut axis_volume = vec![0.0; count];
        for (&a, v) in axes.iter().zip(&material_volume) {
            axis_volume[a] += v;
        }
        let weights = axes
            .iter()
            .zip(&material_volume)
            .map(|(&a, v)| v / axis_volume[a])
            .collect();
        let mut next = 0;
        let slots = axes
            .iter()
            .map(|&a| {
                (members[a] > 1).then(|| {
                    next += 1;
                    next - 1
                })
            })
            .collect();
        Ok(Self {
            classes,
            axes,
            weights,
            slots,
            members,
        })
    }

    /// 物質ごとの全断面積 `sigma` での、[軸] の値（中心）と、軸を共有する物質の ln Σ の中心からの差の絶対値の最大。
    /// 物質が 1 つの軸は、その物質の Σ そのもの（ln と exp を通さないので、ADR-0156 の表と同じ値になる）。
    pub fn centers(&self, sigma: &[f64]) -> (Vec<f64>, f64) {
        let mut log = vec![0.0; self.members.len()];
        let mut centers = vec![0.0; self.members.len()];
        for (m, (&a, &s)) in self.axes.iter().zip(sigma).enumerate() {
            if self.members[a] == 1 {
                centers[a] = s;
            } else {
                log[a] += self.weights[m] * s.ln();
            }
        }
        for (a, c) in centers.iter_mut().enumerate() {
            if self.members[a] > 1 {
                *c = log[a].exp();
            }
        }
        let deviation = self
            .axes
            .iter()
            .zip(sigma)
            .filter(|&(&a, _)| self.members[a] > 1)
            .map(|(&a, &s)| (s.ln() - log[a]).abs())
            // NaN を飛ばさない（中心も NaN になり、範囲の検査で止まる）
            .fold(
                0.0,
                |worst: f64, d| if d.is_nan() || d > worst { d } else { worst },
            );
        (centers, deviation)
    }

    /// 表の格子点の組ごとに求める衝突確率の数: 表の値と、軸を共有する物質ごとの中心差分の 2 つ
    pub fn evaluations(&self) -> usize {
        1 + 2 * self.derivatives()
    }

    /// 微分の数（軸を共有する物質の数）
    fn derivatives(&self) -> usize {
        self.slots.iter().flatten().count()
    }

    /// [リング] の全断面積（真空は 0）。`sigma` は物質ごと
    fn totals(&self, sigma: &[f64]) -> Vec<f64> {
        self.classes
            .iter()
            .map(|m| m.map_or(0.0, |m| sigma[m]))
            .collect()
    }
}

/// 物質ごとの全断面積の格子での飛程の表。
pub struct TrackTable {
    layout: Layout,
    axes: Vec<Axis>,
    /// [格子点の組 × (1 + 微分) × リング × リング] の R_ij と、微分ごとの ∂R_ij/∂ln Σ（格子点の組は最初の軸がいちばん速く変わる）
    values: Vec<f64>,
}

impl TrackTable {
    /// `layout` はリングの物質と軸の対応、`ranges` は軸ごとの全断面積（中心）の [最小, 最大] (1/cm)、
    /// `density` は ln Σ の 10 倍あたりの格子点の数。
    pub fn new(
        cell: &Cylinder,
        layout: Layout,
        ranges: &[[f64; 2]],
        density: usize,
    ) -> Result<Self, TableError> {
        if density == 0 {
            return Err(TableError::Density);
        }
        let rings = cell.volumes().len();
        if layout.classes.len() != rings {
            return Err(TableError::Rings {
                given: layout.classes.len(),
                rings,
            });
        }
        if ranges.len() != layout.members.len() {
            return Err(TableError::Materials {
                given: ranges.len(),
                materials: layout.members.len(),
            });
        }
        let axes = axes(ranges, density)?;
        let nodes = count(&axes)?;
        let block = (1 + layout.derivatives()) * rings * rings;
        let size = nodes.checked_mul(block).ok_or(TableError::Size)?;
        // 格子点の組ごとの R をその場の並びに書く（組ごとの小さな確保と、集めてからつなぐコピーを避ける）
        let mut values = vec![0.0; size];
        values
            .par_chunks_mut(block)
            .enumerate()
            .try_for_each(|(node, out)| {
                let mut rest = node;
                let centers: Vec<f64> = axes
                    .iter()
                    .map(|a| {
                        let k = rest % a.count;
                        rest /= a.count;
                        a.sigma(k)
                    })
                    .collect();
                let sigma: Vec<f64> = layout.axes.iter().map(|&a| centers[a]).collect();
                let (base, derivatives) = out.split_at_mut(rings * rings);
                collisions(cell, &layout.totals(&sigma), base)?;
                // 軸を共有する物質ごとに、その物質のリングの ln Σ だけを ±STEP ずらした R の中心差分
                let mut plus = vec![0.0; rings * rings];
                let mut minus = vec![0.0; rings * rings];
                for (m, slot) in layout.slots.iter().enumerate() {
                    let Some(slot) = slot else { continue };
                    let mut shifted = sigma.clone();
                    shifted[m] = sigma[m] * STEP.exp();
                    collisions(cell, &layout.totals(&shifted), &mut plus)?;
                    shifted[m] = sigma[m] * (-STEP).exp();
                    collisions(cell, &layout.totals(&shifted), &mut minus)?;
                    let d = &mut derivatives[slot * rings * rings..(slot + 1) * rings * rings];
                    for ((d, p), q) in d.iter_mut().zip(&plus).zip(&minus) {
                        *d = (p - q) / (2.0 * STEP);
                    }
                }
                Ok::<(), CpError>(())
            })?;
        tracing::debug!(nodes, derivatives = layout.derivatives(), "track table");
        Ok(Self {
            layout,
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
    /// 軸を共有する物質の ln Σ の中心からの差が [`MAX_DEVIATION`] を超えるときは `None`（格子点ごとに求める）。
    /// 表の範囲の外（外挿になる）と、物質の数の違う `sigma` はエラー。
    pub fn tracks(&self, sigma: &[f64]) -> Result<Option<Vec<Vec<f64>>>, TableError> {
        let materials = self.layout.axes.len();
        if sigma.len() != materials {
            return Err(TableError::Materials {
                given: sigma.len(),
                materials,
            });
        }
        let (centers, deviation) = self.layout.centers(sigma);
        if deviation > MAX_DEVIATION {
            return Ok(None);
        }
        for (material, (axis, &s)) in self.axes.iter().zip(&centers).enumerate() {
            if !(axis.range[0] <= s && s <= axis.range[1]) {
                return Err(TableError::Outside {
                    material,
                    sigma: s,
                    range: axis.range,
                });
            }
        }
        let n = self.layout.classes.len();
        // 表の値と微分ごとの係数: 1 と、微分の物質の ln Σ の中心からの差
        let mut coefficients = vec![1.0];
        for (m, slot) in self.layout.slots.iter().enumerate() {
            if slot.is_some() {
                let a = self.layout.axes[m];
                coefficients.push(sigma[m].ln() - centers[a].ln());
            }
        }
        let block = coefficients.len() * n * n;
        let stencils: Vec<(usize, [f64; STENCIL])> = self
            .axes
            .iter()
            .zip(&centers)
            .map(|(a, &s)| a.stencil(s))
            .collect();
        let mut r = vec![0.0; n * n];
        // 軸ごとの内挿の点の組を、最初の軸がいちばん速く変わる順に回る
        let mut offsets = vec![0usize; self.axes.len()];
        'combos: loop {
            let (mut node, mut stride, mut weight) = (0, 1, 1.0);
            for ((axis, &(first, w)), &k) in self.axes.iter().zip(&stencils).zip(&offsets) {
                node += (first + k) * stride;
                stride *= axis.count;
                weight *= w[k];
            }
            let values = &self.values[node * block..(node + 1) * block];
            for (c, v) in coefficients.iter().zip(values.chunks(n * n)) {
                let factor = weight * c;
                for (x, v) in r.iter_mut().zip(v) {
                    *x += factor * v;
                }
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
        let totals = self.layout.totals(sigma);
        Ok(Some(
            (0..n)
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
                .collect(),
        ))
    }
}

/// リングの全断面積 `totals` での R_ij = W_ij Σ_j（真空の列 j は W_ij）を `out`（[リング × リング]）に書く
fn collisions(cell: &Cylinder, totals: &[f64], out: &mut [f64]) -> Result<(), CpError> {
    let w = cell.tracks(totals)?;
    for (row, out) in w.iter().zip(out.chunks_mut(totals.len())) {
        for ((x, w), &t) in out.iter_mut().zip(row).zip(totals) {
            *x = if t > 0.0 { w * t } else { *w };
        }
    }
    Ok(())
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

#[cfg(test)]
mod tests;
