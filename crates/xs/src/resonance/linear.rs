//! 共鳴の計算で解く、小さな複素数の連立一次方程式。

use num_complex::Complex64;

/// 係数行列が特異だった。
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Singular;

/// 複素 n×n 行列の方程式 M X = B を部分ピボット付き Gauss-Jordan 法で解き、`b` を X で置き換える
/// （どちらも行優先の n×n。`m` は壊れる）。同じ大きさのピボットの候補があれば後の行を選ぶ。
/// ピボットが 0 か NaN なら `Singular`。
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
        // 0 と NaN を弾く（NaN は total_cmp で最大になるので、ピボットに選ばれる）
        let size = m[pivot * n + col].norm_sqr();
        if size == 0.0 || size.is_nan() {
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

/// 極の共鳴を含む M = M' + c w wᵀ/d（w は実のベクトル、c は 0 でない定数）の、d → 0 での M⁻¹ の極限を
/// Sherman-Morrison の公式 M⁻¹ → M'⁻¹ − (M'⁻¹w)(wᵀM'⁻¹) / (wᵀM'⁻¹w) で求め、`inverse` に入れた
/// M'⁻¹ を置き換える（行優先の n×n）。M' は複素対称とする（wᵀM'⁻¹ = (M'⁻¹w)ᵀ）。
/// w ≠ 0 なら、M' の実部が正定値なので wᵀM'⁻¹w の実部は正になる。0 になったら `Singular`
/// （幅のない極の共鳴は、呼び出す前に除いておく）。
pub(super) fn pole_limit(inverse: &mut [Complex64], w: &[f64], n: usize) -> Result<(), Singular> {
    let u: Vec<Complex64> = (0..n)
        .map(|a| (0..n).fold(Complex64::ZERO, |acc, b| acc + inverse[a * n + b] * w[b]))
        .collect();
    let s = (0..n).fold(Complex64::ZERO, |acc, a| acc + u[a] * w[a]);
    if s == Complex64::ZERO {
        return Err(Singular);
    }
    let scale = s.inv();
    for a in 0..n {
        for b in 0..n {
            inverse[a * n + b] -= u[a] * u[b] * scale;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
