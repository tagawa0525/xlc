//! 小さな実数の密行列の連立一次方程式（lattice の CMFD の群×群のブロック、UFG と SPH のリングの間の結合、B1）。

/// 部分ピボット選択付きの LU 分解。
#[derive(Debug, Clone)]
pub struct Lu {
    n: usize,
    /// 行優先の L（対角は 1、格納しない）と U
    lu: Vec<f64>,
    /// 列 col の消去で行 col と入れ替えた行（LAPACK の ipiv と同じく、順に適用する）
    swaps: Vec<usize>,
}

impl Lu {
    /// 行優先の n × n 行列を分解する。ピボットが 0 か NaN なら特異として `None`。
    pub fn new(mut a: Vec<f64>, n: usize) -> Option<Self> {
        let mut swaps = Vec::with_capacity(n);
        for col in 0..n {
            let best =
                (col..n).max_by(|&i, &j| a[i * n + col].abs().total_cmp(&a[j * n + col].abs()))?;
            // 0 と NaN を弾く（NaN は total_cmp で最大になるので、ピボットに選ばれる）
            let pivot = a[best * n + col];
            if pivot == 0.0 || pivot.is_nan() {
                return None;
            }
            if best != col {
                for j in 0..n {
                    a.swap(best * n + j, col * n + j);
                }
            }
            swaps.push(best);
            for row in col + 1..n {
                let factor = a[row * n + col] / a[col * n + col];
                a[row * n + col] = factor;
                for j in col + 1..n {
                    a[row * n + j] -= factor * a[col * n + j];
                }
            }
        }
        Some(Self { n, lu: a, swaps })
    }

    /// `b` を A x = b の解 x で置き換える。
    ///
    /// CMFD の反復ごとに全セルで呼ぶので、確保をせず、行ごとの内積で前進・後退代入する。
    pub fn solve(&self, b: &mut [f64]) {
        let n = self.n;
        for (col, &row) in self.swaps.iter().enumerate() {
            if row != col {
                b.swap(col, row);
            }
        }
        for row in 0..n {
            let l = &self.lu[row * n..][..row];
            let sum: f64 = l.iter().zip(&b[..row]).map(|(a, x)| a * x).sum();
            b[row] -= sum;
        }
        for row in (0..n).rev() {
            let u = &self.lu[row * n + row + 1..][..n - row - 1];
            let sum: f64 = u.iter().zip(&b[row + 1..]).map(|(a, x)| a * x).sum();
            b[row] = (b[row] - sum) / self.lu[row * n + row];
        }
    }
}

/// 行優先の n × n 行列 `a` で A x = b を一度だけ解く。特異か、解が有限でなければ `None`。
pub fn solve(a: Vec<f64>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let lu = Lu::new(a, b.len())?;
    lu.solve(&mut b);
    b.iter().all(|v| v.is_finite()).then_some(b)
}

#[cfg(test)]
mod tests;
