//! 共鳴の計算で解く、小さな複素数の連立一次方程式。

use num_complex::Complex64;

/// 係数行列が特異だった。
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Singular;

/// 複素 n×n 行列の方程式 M X = B を部分ピボット付き Gauss-Jordan 法で解き、`b` を X で置き換える
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

#[cfg(test)]
mod tests {
    use super::*;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    #[test]
    fn solves_a_system_that_needs_pivoting() {
        // 先頭の対角が 0 なので行の入れ替えが要る。右辺は2列
        let m = [
            c(0.0, 0.0),
            c(2.0, 1.0),
            c(1.0, 0.0),
            c(1.0, -1.0),
            c(0.5, 0.0),
            c(0.0, 3.0),
            c(4.0, 0.0),
            c(-1.0, 2.0),
            c(2.0, 2.0),
        ];
        let b = [
            c(1.0, 0.0),
            c(0.0, 1.0),
            c(2.0, -1.0),
            c(3.0, 0.0),
            c(-1.0, 1.0),
            c(0.0, 0.5),
            c(1.0, 1.0),
            c(2.0, 0.0),
            c(0.0, -2.0),
        ];
        let mut work = m;
        let mut x = b;
        solve(&mut work, &mut x, 3).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let mx: Complex64 = (0..3).map(|k| m[i * 3 + k] * x[k * 3 + j]).sum();
                assert!((mx - b[i * 3 + j]).norm() < 1e-14, "({i}, {j}): {mx}");
            }
        }
    }

    #[test]
    fn rejects_a_zero_pivot() {
        let mut m = [c(1.0, 0.0), c(2.0, 0.0), c(2.0, 0.0), c(4.0, 0.0)];
        let mut b = [c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(1.0, 0.0)];
        assert_eq!(solve(&mut m, &mut b, 2), Err(Singular));
    }

    #[test]
    fn rejects_a_nan_pivot() {
        let mut m = [c(f64::NAN, 0.0), c(1.0, 0.0), c(1.0, 0.0), c(1.0, 0.0)];
        let mut b = [c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(1.0, 0.0)];
        assert_eq!(solve(&mut m, &mut b, 2), Err(Singular));
    }
}
