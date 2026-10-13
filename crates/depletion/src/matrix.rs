//! 燃焼行列（dn/dt = A n）の疎な形。

use crate::DepletionError;

/// 同じ位置の寄与の和を 0 とみなす、寄与の絶対値の和に対する割合。m 個の和の丸めは (m − 1) ε × 絶対値の和で
/// 抑えられるので、16 ε は 17 個までの寄与を覆う（1 つの位置の寄与は、崩壊と反応の数種）
const ROUNDING: f64 = 16.0 * f64::EPSILON;

/// 行の順に並んだ (行, 値) を行ごとに足し、(行, 和, 寄与の絶対値の和) にする
fn sum_by_row(column: &[(usize, f64)]) -> Vec<(usize, f64, f64)> {
    let mut sums: Vec<(usize, f64, f64)> = Vec::new();
    for &(row, value) in column {
        match sums.last_mut() {
            Some(last) if last.0 == row => {
                last.1 += value;
                last.2 += value.abs();
            }
            _ => sums.push((row, value, value.abs())),
        }
    }
    sums
}

/// 燃焼行列 A。列 j が親の核種、行 i が生まれる核種で、A[i][j] は核種 j の 1 個が単位時間に核種 i を生む割合
/// （1/s）。対角は核種の消える割合の符号を変えたもの。OpenMC の `A[j, i]`（i から j へ）と同じ向き。
///
/// 圧縮列の形（CSC）で持ち、各列の行は昇順で、対角は値が 0 でも必ず持つ（CRAM の A t − θ I が対角を使うため）。
#[derive(Debug, Clone)]
pub struct DepletionMatrix {
    n: usize,
    col_ptr: Vec<usize>,
    row_idx: Vec<usize>,
    values: Vec<f64>,
    /// 列 j の対角の、`row_idx` と `values` の中の位置
    diagonal: Vec<usize>,
}

impl DepletionMatrix {
    /// (行, 列, 値) の並びから作る。同じ位置の値は足す（同じ親から同じ核種への反応と崩壊を別々に書ける）。
    ///
    /// 足したあとの非対角が負（負の生成）、対角が正（負の消滅）ならエラーにする。燃焼行列はこの符号を必ず
    /// 満たすので、外れるのは向き（転置）か符号の取り違えである。
    pub fn from_triplets(
        n: usize,
        entries: &[(usize, usize, f64)],
    ) -> Result<Self, DepletionError> {
        let mut columns: Vec<Vec<(usize, f64)>> = (0..n).map(|j| vec![(j, 0.0)]).collect();
        for &(row, col, value) in entries {
            if row >= n || col >= n {
                return Err(DepletionError::OutOfRange { row, col, n });
            }
            if !value.is_finite() {
                return Err(DepletionError::NonFiniteEntry { row, col, value });
            }
            columns[col].push((row, value));
        }
        let mut col_ptr = Vec::with_capacity(n + 1);
        let mut row_idx = Vec::new();
        let mut values = Vec::new();
        let mut diagonal = Vec::with_capacity(n);
        col_ptr.push(0);
        for (col, mut column) in columns.into_iter().enumerate() {
            column.sort_by_key(|&(row, _)| row);
            for (row, value, magnitude) in sum_by_row(&column) {
                // 1 つずつは有限でも、同じ位置の和があふれうる
                if !value.is_finite() {
                    return Err(DepletionError::NonFiniteEntry { row, col, value });
                }
                // 打ち消し合う寄与の和は、寄与の大きさの丸め（数 ε）だけ符号を外れうるので、その分は 0 とみなす
                let value = if value.abs() <= ROUNDING * magnitude {
                    0.0
                } else {
                    value
                };
                if row == col {
                    if value > 0.0 {
                        return Err(DepletionError::PositiveLoss {
                            nuclide: col,
                            value,
                        });
                    }
                    diagonal.push(row_idx.len());
                } else if value < 0.0 {
                    return Err(DepletionError::NegativeProduction { row, col, value });
                }
                row_idx.push(row);
                values.push(value);
            }
            col_ptr.push(row_idx.len());
        }
        Ok(Self {
            n,
            col_ptr,
            row_idx,
            values,
            diagonal,
        })
    }

    /// 核種の数
    pub fn len(&self) -> usize {
        self.n
    }

    /// 核種がない
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    pub(crate) fn col_ptr(&self) -> &[usize] {
        &self.col_ptr
    }

    pub(crate) fn row_idx(&self) -> &[usize] {
        &self.row_idx
    }

    pub(crate) fn values(&self) -> &[f64] {
        &self.values
    }

    pub(crate) fn diagonal(&self) -> &[usize] {
        &self.diagonal
    }
}

#[cfg(test)]
mod tests;
