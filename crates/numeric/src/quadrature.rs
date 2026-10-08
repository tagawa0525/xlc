//! Legendre 多項式と Gauss–Legendre の求積。

use std::f64::consts::PI;

/// Gauss–Legendre の点が求まらない。
#[derive(Debug, thiserror::Error)]
#[error("Gauss-Legendre nodes for {0} points did not converge")]
pub struct NodesError(pub usize);

/// Legendre 多項式 P_0(x)〜P_order(x)
pub fn legendre(x: f64, order: usize) -> Vec<f64> {
    let mut p = Vec::with_capacity(order + 1);
    p.push(1.0);
    if order >= 1 {
        p.push(x);
    }
    for l in 2..=order {
        let l_f = l as f64;
        p.push(((2.0 * l_f - 1.0) * x * p[l - 1] - (l_f - 1.0) * p[l - 2]) / l_f);
    }
    p
}

/// [−1, 1] の n 点の Gauss–Legendre の点と重み（Newton 法で Legendre 多項式の零点を求める）。
pub fn gauss_legendre(n: usize) -> Result<Vec<(f64, f64)>, NodesError> {
    let mut nodes = Vec::with_capacity(n);
    for i in 0..n {
        let mut x = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut converged = false;
        for _ in 0..100 {
            let (p, dp) = legendre_with_derivative(n, x);
            let step = p / dp;
            x -= step;
            if step.abs() <= 1e-15 {
                converged = true;
                break;
            }
        }
        if !converged {
            return Err(NodesError(n));
        }
        let (_, dp) = legendre_with_derivative(n, x);
        nodes.push((x, 2.0 / ((1.0 - x * x) * dp * dp)));
    }
    Ok(nodes)
}

/// Legendre 多項式 P_n(x) とその導関数
fn legendre_with_derivative(n: usize, x: f64) -> (f64, f64) {
    let (mut p0, mut p1) = (1.0, x);
    for k in 2..=n {
        let k = k as f64;
        (p0, p1) = (p1, ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k);
    }
    if n == 0 {
        return (1.0, 0.0);
    }
    (p1, n as f64 * (x * p1 - p0) / (x * x - 1.0))
}

/// 2 点の Gauss–Legendre の点と重み（[−1, 1]、昇順）
pub const GAUSS_LEGENDRE_2: [(f64, f64); 2] = [
    (-0.577_350_269_189_625_8, 1.0),
    (0.577_350_269_189_625_8, 1.0),
];

/// 4 点の Gauss–Legendre の点と重み（[−1, 1]、昇順）
pub const GAUSS_LEGENDRE_4: [(f64, f64); 4] = [
    (-0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
    (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
];

/// 8 点の Gauss–Legendre の点と重み（[−1, 1]、昇順）
pub const GAUSS_LEGENDRE_8: [(f64, f64); 8] = [
    (-0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
    (-0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (-0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (-0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn legendre_matches_closed_forms() {
        let x = 0.3_f64;
        let p = legendre(x, 3);
        assert_eq!(p[0], 1.0);
        assert_eq!(p[1], x);
        assert!((p[2] - 0.5 * (3.0 * x * x - 1.0)).abs() < 1e-15);
        assert!((p[3] - 0.5 * (5.0 * x * x * x - 3.0 * x)).abs() < 1e-15);
    }

    #[test]
    fn tables_match_the_computed_nodes() {
        for table in [
            &GAUSS_LEGENDRE_2[..],
            &GAUSS_LEGENDRE_4[..],
            &GAUSS_LEGENDRE_8[..],
        ] {
            let mut computed = gauss_legendre(table.len()).unwrap();
            computed.sort_by(|a, b| a.0.total_cmp(&b.0));
            for ((x, w), (cx, cw)) in table.iter().zip(&computed) {
                assert!((x - cx).abs() < 1e-15 && (w - cw).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn zero_points_is_an_error() {
        assert!(gauss_legendre(0).is_err());
    }

    #[test]
    fn nodes_integrate_polynomials_exactly() {
        // n 点は 2n − 1 次までの多項式を正確に積分する
        for n in [1, 3, 7, 20] {
            let nodes = gauss_legendre(n).unwrap();
            let degree = 2 * n - 2;
            let integral: f64 = nodes.iter().map(|&(x, w)| w * x.powi(degree as i32)).sum();
            assert!((integral - 2.0 / (degree as f64 + 1.0)).abs() < 1e-13);
        }
    }
}
