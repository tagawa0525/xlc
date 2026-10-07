//! 共鳴の計算で解く、小さな複素数の連立一次方程式。

use num_complex::Complex64;

/// 係数行列が特異だった。
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Singular;

/// 複素 n×n 行列の方程式 M X = B を部分ピボットつき Gauss-Jordan 法で解き、`b` を X で置き換える
/// （どちらも行優先の n×n。`m` は壊れる）。同じ大きさのピボットの候補があれば後の行を選ぶ。
pub(super) fn solve(m: &mut [Complex64], b: &mut [Complex64], n: usize) -> Result<(), Singular> {
    for col in 0..n {
        let pivot = (col + 1..n).fold(col, |best, i| {
            if m[i * n + col]
                .norm_sqr()
                .total_cmp(&m[best * n + col].norm_sqr())
                .is_ge()
            {
                i
            } else {
                best
            }
        });
        if m[pivot * n + col].norm_sqr() == 0.0 {
            return Err(Singular);
        }
        for j in 0..n {
            m.swap(col * n + j, pivot * n + j);
            b.swap(col * n + j, pivot * n + j);
        }
        let p = m[col * n + col].inv();
        for j in 0..n {
            m[col * n + j] *= p;
            b[col * n + j] *= p;
        }
        for row in 0..n {
            if row != col {
                let f = m[row * n + col];
                for j in 0..n {
                    m[row * n + j] -= f * m[col * n + j];
                    b[row * n + j] -= f * b[col * n + j];
                }
            }
        }
    }
    Ok(())
}
