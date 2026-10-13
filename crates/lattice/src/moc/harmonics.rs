//! 実球面調和関数（散乱の非等方性を角度について展開する基底）。
//!
//! 半正規化（Schmidt）: Y_lm(Ω) = √((2 - δ_m0) (l - |m|)! / (l + |m|)!) P_l^|m|(μ) × (m ≥ 0 なら cos mφ、
//! m < 0 なら sin |m|φ)。Y_00 = 1 で、加法定理 Σ_m Y_lm(Ω) Y_lm(Ω') = P_l(Ω·Ω') を満たす。

/// 2次元（z について対称）の問題で残る (l, m) の組（l + m が偶数。奇数の組は μ について奇関数で、
/// 角度中性子束のモーメントが 0 になる）。l = 0..=`order`、各 l の中は m の昇順
pub(crate) fn planar_moments(order: usize) -> Vec<(usize, isize)> {
    (0..=order)
        .flat_map(|l| {
            let l_signed = l as isize;
            (-l_signed..=l_signed)
                .filter(move |m| (l_signed + m) % 2 == 0)
                .map(move |m| (l, m))
        })
        .collect()
}

/// 方向 Ω（極角の余弦 `mu`、方位角 `phi`）での Y_lm。|m| ≤ l であること
pub(crate) fn real_harmonic(l: usize, m: isize, mu: f64, phi: f64) -> f64 {
    let order = m.unsigned_abs();
    // (l - |m|)! / (l + |m|)!
    let ratio: f64 = (l - order + 1..=l + order)
        .map(|k| 1.0 / k as f64)
        .product();
    let norm = if m == 0 { 1.0 } else { (2.0 * ratio).sqrt() };
    let angle = order as f64 * phi;
    let trig = if m >= 0 { angle.cos() } else { angle.sin() };
    norm * associated_legendre(l, order, mu) * trig
}

/// ルジャンドル陪関数 P_l^m(x)（Condon–Shortley の位相 (-1)^m を含まない）。m ≤ l であること
fn associated_legendre(l: usize, m: usize, x: f64) -> f64 {
    // P_m^m = (2m - 1)!! (1 - x²)^(m/2)
    let sin = (1.0 - x * x).max(0.0).sqrt();
    let mut diagonal = 1.0;
    for k in 0..m {
        diagonal *= (2 * k + 1) as f64 * sin;
    }
    if l == m {
        return diagonal;
    }
    // (l - m) P_l^m = (2l - 1) x P_(l-1)^m - (l + m - 1) P_(l-2)^m
    let (mut previous, mut current) = (diagonal, x * (2 * m + 1) as f64 * diagonal);
    for n in m + 2..=l {
        let next =
            ((2 * n - 1) as f64 * x * current - (n + m - 1) as f64 * previous) / (n - m) as f64;
        previous = current;
        current = next;
    }
    current
}

#[cfg(test)]
mod tests;
