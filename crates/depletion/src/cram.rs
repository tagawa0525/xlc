//! CRAM（Chebyshev 有理近似）の不完全部分分数形（IPF）で exp(A t) n(0) を求める（ADR-0161）。
//!
//! 次数 k の近似 r(z) ≈ exp(z) を、Pusa (2016) の積の形
//! r(z) = α0 ∏_{l=1}^{k/2} (1 + 2 Re[α̃_l / (z − θ_l)]) で行列に当てる。各因子は y ← y + 2 Re[α̃_l (A t − θ_l I)⁻¹ y]
//! で、複素の疎な連立方程式を極ごとに 1 回解く。極ごとの行列は対角だけが違い、非零の位置は同じなので、記号的な
//! 分解は 1 回だけ行う。

mod coefficients;

use faer::Mat;
use faer::linalg::solvers::Solve;
use faer::sparse::linalg::solvers::{Lu, SymbolicLu};
use faer::sparse::{SparseColMatRef, SymbolicSparseColMat};
use num_complex::Complex64;

use crate::{DepletionError, DepletionMatrix};

/// CRAM の次数。次数 k は k/2 個の極で、1 回の `deplete` で疎な複素の連立方程式を k/2 回解く。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CramOrder {
    /// 次数 16（8 個の極）
    Order16,
    /// 次数 48（24 個の極）
    Order48,
}

/// IPF の極 θ、留数 α̃（どちらも（実部、虚部））と α0
struct Coefficients {
    theta: &'static [(f64, f64)],
    alpha: &'static [(f64, f64)],
    alpha0: f64,
}

impl CramOrder {
    fn coefficients(self) -> Coefficients {
        match self {
            Self::Order16 => Coefficients {
                theta: &coefficients::THETA_16,
                alpha: &coefficients::ALPHA_16,
                alpha0: coefficients::ALPHA0_16,
            },
            Self::Order48 => Coefficients {
                theta: &coefficients::THETA_48,
                alpha: &coefficients::ALPHA_48,
                alpha0: coefficients::ALPHA0_48,
            },
        }
    }
}

/// 時間 `t`（s）の後の数密度 exp(A t) n(0) を返す。
///
/// 数密度の単位は問わない（線形なので、返す値は `n0` と同じ単位）。近似の誤差で、ほぼ 0 の核種が小さな負の値に
/// なることがある（次の段の入力として、そのまま渡してよい）。
pub fn deplete(
    a: &DepletionMatrix,
    n0: &[f64],
    t: f64,
    order: CramOrder,
) -> Result<Vec<f64>, DepletionError> {
    let n = a.len();
    if n0.len() != n {
        return Err(DepletionError::LengthMismatch {
            expected: n,
            found: n0.len(),
        });
    }
    if let Some((nuclide, &value)) = n0.iter().enumerate().find(|(_, v)| !v.is_finite()) {
        return Err(DepletionError::NonFiniteDensity { nuclide, value });
    }
    if !(t.is_finite() && t >= 0.0) {
        return Err(DepletionError::InvalidTime(t));
    }
    // exp(0) = I。CRAM の r(0) は 1 から近似の誤差だけずれるので、当てない
    if t == 0.0 {
        return Ok(n0.to_vec());
    }

    let col_ptr = a.col_ptr();
    let row_idx = a.row_idx();
    let mut scaled = Vec::with_capacity(a.values().len());
    for (col, range) in col_ptr.windows(2).enumerate() {
        let (rows, values) = (
            &row_idx[range[0]..range[1]],
            &a.values()[range[0]..range[1]],
        );
        for (&row, &value) in rows.iter().zip(values) {
            let value = value * t;
            if !value.is_finite() {
                return Err(DepletionError::Overflow { row, col });
            }
            scaled.push(value);
        }
    }
    let symbolic =
        SymbolicSparseColMat::new_checked(n, n, col_ptr.to_vec(), None, row_idx.to_vec());
    let lu_symbolic =
        SymbolicLu::try_new(symbolic.as_ref()).map_err(DepletionError::SymbolicFactorization)?;

    let Coefficients {
        theta,
        alpha,
        alpha0,
    } = order.coefficients();
    let mut y = n0.to_vec();
    // 極ごとに変わるのは対角だけなので、非対角は 1 回だけ複素にし、対角を書き換える
    let mut values: Vec<Complex64> = scaled.iter().map(|&s| Complex64::new(s, 0.0)).collect();
    for (&(theta_re, theta_im), &(alpha_re, alpha_im)) in theta.iter().zip(alpha) {
        let theta = Complex64::new(theta_re, theta_im);
        let alpha = Complex64::new(alpha_re, alpha_im);
        for &k in a.diagonal() {
            values[k] = Complex64::new(scaled[k], 0.0) - theta;
        }
        let system = SparseColMatRef::new(symbolic.as_ref(), &values);
        let lu = Lu::try_new_with_symbolic(lu_symbolic.clone(), system)
            .map_err(DepletionError::NumericFactorization)?;
        let mut x = Mat::from_fn(n, 1, |i, _| Complex64::new(y[i], 0.0));
        lu.solve_in_place(x.as_mut());
        for (i, yi) in y.iter_mut().enumerate() {
            *yi += 2.0 * (alpha * x[(i, 0)]).re;
        }
    }
    for (i, yi) in y.iter_mut().enumerate() {
        *yi *= alpha0;
        if !yi.is_finite() {
            return Err(DepletionError::NonFiniteResult(i));
        }
    }
    Ok(y)
}

#[cfg(test)]
mod tests;
