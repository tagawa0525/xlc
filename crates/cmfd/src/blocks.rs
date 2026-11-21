//! セルごとの対角のブロックの LU 分解を、すべてのセルで 1 つの平らな配列に持つ。
//!
//! 演算は `xlc_numeric::lu::Lu` と同じ（部分ピボット選択、行ごとの内積の前進・後退代入）で、解は同じ値になる。行の入れ替えは、
//! 分解した行に当たる元の行の並びにして、右辺を読むときに適用する。`Lu` をセルごとに持つ形は、core の全炉心の BOC の固有値計算が
//! この形より 7 % 長かった（ADR-0275）。`xlc-numeric` に平らな形を足さないのは、それに依存する xs の実行ファイルのバイトを
//! 変えないため（ADR-0232）。

use rayon::prelude::*;

/// 並列に分解するときの、1 つのスレッドがまとめて分解するセルの数の下限
const MIN_CELLS_PER_TASK: usize = 64;

/// [セル] の n × n の LU 分解
#[derive(Debug, Clone)]
pub(crate) struct Blocks {
    n: usize,
    /// [セル × n × n] の行優先の L（対角は 1、格納しない）と U
    lu: Vec<f64>,
    /// [セル × n] の、分解した行 row に当たる元の行（部分ピボット選択の入れ替えを順に適用した並び）
    order: Vec<usize>,
}

impl Blocks {
    /// 行優先の n × n 行列 `matrices` を分解する。ピボットが 0 か NaN のセルがあれば、セルの順で最初のものの番号を `Err` で返す
    pub(crate) fn new(matrices: Vec<Vec<f64>>, n: usize) -> Result<Self, usize> {
        let cells = matrices.len();
        let mut lu = Vec::with_capacity(cells * n * n);
        for a in matrices {
            lu.extend(a);
        }
        let mut order: Vec<usize> = (0..cells).flat_map(|_| 0..n).collect();
        // セルごとの分解は互いに独立なので並列に行う（どのスレッドで分解しても同じ値）
        let singular = lu
            .par_chunks_mut(n * n)
            .zip(order.par_chunks_mut(n))
            .enumerate()
            .with_min_len(MIN_CELLS_PER_TASK)
            .filter_map(|(cell, (a, order))| (!factor(a, order, n)).then_some(cell))
            .min();
        match singular {
            Some(cell) => Err(cell),
            None => Ok(Self { n, lu, order }),
        }
    }

    /// セル `cell` の A x = `b` の解を `x` に書く
    pub(crate) fn solve(&self, cell: usize, b: &[f64], x: &mut [f64]) {
        let n = self.n;
        let lu = &self.lu[cell * n * n..][..n * n];
        let order = &self.order[cell * n..][..n];
        for row in 0..n {
            let l = &lu[row * n..][..row];
            let sum: f64 = l.iter().zip(&x[..row]).map(|(a, y)| a * y).sum();
            x[row] = b[order[row]] - sum;
        }
        for row in (0..n).rev() {
            let u = &lu[row * n + row + 1..][..n - row - 1];
            let sum: f64 = u.iter().zip(&x[row + 1..]).map(|(a, y)| a * y).sum();
            x[row] = (x[row] - sum) / lu[row * n + row];
        }
    }
}

/// 行優先の n × n 行列 `a` をその場で分解し、`order`（0..n で始める）に行の並びを書く。ピボットが 0 か NaN なら偽
fn factor(a: &mut [f64], order: &mut [usize], n: usize) -> bool {
    for col in 0..n {
        let Some(best) =
            (col..n).max_by(|&i, &j| a[i * n + col].abs().total_cmp(&a[j * n + col].abs()))
        else {
            return false;
        };
        // 0 と NaN を弾く（NaN は total_cmp で最大になるので、ピボットに選ばれる）
        let pivot = a[best * n + col];
        if pivot == 0.0 || pivot.is_nan() {
            return false;
        }
        if best != col {
            for j in 0..n {
                a.swap(best * n + j, col * n + j);
            }
            order.swap(best, col);
        }
        for row in col + 1..n {
            let factor = a[row * n + col] / a[col * n + col];
            a[row * n + col] = factor;
            for j in col + 1..n {
                a[row * n + j] -= factor * a[col * n + j];
            }
        }
    }
    true
}
