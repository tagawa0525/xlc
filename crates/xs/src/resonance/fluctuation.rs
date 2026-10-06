//! 非分離共鳴の幅のゆらぎ（平均 1 の χ² 分布）による期待値。
//!
//! 自由度 ν の χ² 分布の密度は k^k x^(k-1) e^(-kx) / Γ(k)（k = ν/2）。x = u² と置くと
//! 2 k^k u^(2k-1) e^(-k u²) / Γ(k) du になり、ν = 1 でも原点で有限になるので、u について適応積分する。

/// 自由度 `dof` の χ² 分布での f(x) の期待値。自由度 0 は幅がゆらがない（x = 1）ことを表す。
pub(super) fn expectation(dof: u32, f: impl Fn(f64) -> f64) -> f64 {
    let _ = (dof, f);
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::special::erfc;

    const SQRT_PI: f64 = 1.772_453_850_905_516;

    #[test]
    fn reproduces_moments_of_chi_square() {
        // 平均 1 の χ²_ν の k 次の積率は Π_{i<k} (1 + 2i/ν)
        for dof in 1..=4 {
            let nu = f64::from(dof);
            let mut expected = 1.0;
            for k in 0..=6 {
                let moment = expectation(dof, |x| x.powi(k));
                assert!(
                    (moment - expected).abs() <= 1e-9 * expected,
                    "ν={dof}, k={k}: {moment} vs {expected}"
                );
                expected *= 1.0 + 2.0 * f64::from(k) / nu;
            }
        }
    }

    #[test]
    fn matches_closed_form_near_the_pole() {
        // E[1/(x+c)] = k^k c^(k-1) e^(kc) Γ(1-k, kc)。
        // Γ(1/2, z) = √π erfc(√z)、Γ(-1/2, z) = 2 e^(-z)/√z - 2√π erfc(√z)
        for c in [1.0e-4_f64, 1.0e-2, 1.0, 10.0] {
            let z = 0.5 * c;
            let nu1 = 0.5_f64.sqrt() / c.sqrt() * z.exp() * SQRT_PI * erfc(z.sqrt());
            let z = 1.5 * c;
            let gamma = 2.0 * (-z).exp() / z.sqrt() - 2.0 * SQRT_PI * erfc(z.sqrt());
            let nu3 = 1.5_f64.powf(1.5) * c.sqrt() * z.exp() * gamma;
            for (dof, exact) in [(1, nu1), (3, nu3)] {
                let actual = expectation(dof, |x| 1.0 / (x + c));
                assert!(
                    (actual - exact).abs() <= 1e-8 * exact,
                    "ν={dof}, c={c}: {actual} vs {exact}"
                );
            }
        }
    }

    #[test]
    fn zero_degrees_of_freedom_do_not_fluctuate() {
        assert_eq!(expectation(0, |x| 3.0 * x + 1.0), 4.0);
    }
}
