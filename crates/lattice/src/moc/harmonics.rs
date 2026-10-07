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
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// ルジャンドル多項式 P_l(x)
    fn legendre(l: usize, x: f64) -> f64 {
        let (mut previous, mut current) = (1.0, x);
        if l == 0 {
            return previous;
        }
        for n in 1..l {
            let n = n as f64;
            let next = ((2.0 * n + 1.0) * x * current - n * previous) / (n + 1.0);
            previous = current;
            current = next;
        }
        current
    }

    #[test]
    fn low_orders_are_the_direction_cosines() {
        let (mu, phi) = (0.3_f64, 1.1_f64);
        let sin = (1.0 - mu * mu).sqrt();
        assert_eq!(real_harmonic(0, 0, mu, phi), 1.0);
        assert!((real_harmonic(1, 0, mu, phi) - mu).abs() < 1e-15);
        assert!((real_harmonic(1, 1, mu, phi) - sin * phi.cos()).abs() < 1e-15);
        assert!((real_harmonic(1, -1, mu, phi) - sin * phi.sin()).abs() < 1e-15);
    }

    #[test]
    fn addition_theorem_holds() {
        let directions: [(f64, f64); 4] = [(0.3, 1.1), (-0.7, 4.0), (0.95, 2.5), (0.0, 0.2)];
        for (i, &(mu1, phi1)) in directions.iter().enumerate() {
            for &(mu2, phi2) in &directions[i..] {
                let s1 = (1.0 - mu1 * mu1).sqrt();
                let s2 = (1.0 - mu2 * mu2).sqrt();
                let cosine = mu1 * mu2 + s1 * s2 * (phi1 - phi2).cos();
                for l in 0..=5 {
                    let sum: f64 = (-(l as isize)..=l as isize)
                        .map(|m| real_harmonic(l, m, mu1, phi1) * real_harmonic(l, m, mu2, phi2))
                        .sum();
                    let expected = legendre(l, cosine);
                    assert!((sum - expected).abs() < 1e-13, "l {l}: {sum} vs {expected}");
                }
            }
        }
    }

    #[test]
    fn planar_moments_are_even_in_mu() {
        assert_eq!(planar_moments(0), vec![(0, 0)]);
        assert_eq!(
            planar_moments(3),
            vec![
                (0, 0),
                (1, -1),
                (1, 1),
                (2, -2),
                (2, 0),
                (2, 2),
                (3, -3),
                (3, -1),
                (3, 1),
                (3, 3)
            ]
        );
        let phi = 0.7 * PI;
        for l in 0..=4_usize {
            for m in -(l as isize)..=l as isize {
                let (up, down) = (
                    real_harmonic(l, m, 0.4, phi),
                    real_harmonic(l, m, -0.4, phi),
                );
                let even = planar_moments(l).contains(&(l, m));
                let expected = if even { up } else { -up };
                assert!((down - expected).abs() < 1e-15, "({l}, {m})");
            }
        }
    }
}
