//! ドップラー拡がりした共鳴の線形状と、ETOX 法の J, K 積分（NJOY UNRESR の ajku に相当）。
//!
//! x = 2(E - E_r)/Γ、θ = Γ/Δ（Δ はドップラー幅）として、線形状は
//! ψ(x, θ) = √π (θ/2) Re w(θx/2 + iθ/2)（θ → ∞ で 1/(1 + x²)）。
//! J(β, θ) = ∫₀^∞ ψ/(ψ + β) dx、K(β, θ) = ∫₀^∞ ψ(ψ + 2β)/(ψ + β)² dx。

/// 線形状 ψ(x, θ)。
pub(super) fn psi(x: f64, theta: f64) -> f64 {
    let _ = (x, theta);
    todo!()
}

/// 各 β について (J(β, θ), K(β, θ)) を求める。β は正であること。
pub(super) fn jk_integrals(betas: &[f64], theta: f64) -> Vec<(f64, f64)> {
    let _ = (betas, theta);
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn assert_rel(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() <= tol * expected.abs(),
            "{actual} vs {expected}"
        );
    }

    #[test]
    fn natural_line_shape_without_doppler() {
        // θ → ∞ で ψ = 1/(1+x²)。a² = (1+β)/β として
        // J = π/(2aβ)、K = (1/β)[π/(2a) + π(a²+1)/(4a³)]
        let theta = 1.0e5;
        let betas = [1.0e-4, 1.0e-2, 1.0, 1.0e2];
        for (&beta, (j, k)) in betas.iter().zip(jk_integrals(&betas, theta)) {
            let a = ((1.0 + beta) / beta).sqrt();
            assert_rel(j, PI / (2.0 * a * beta), 1e-8);
            assert_rel(
                k,
                (PI / (2.0 * a) + PI * (a * a + 1.0) / (4.0 * a * a * a)) / beta,
                1e-8,
            );
        }
        assert_rel(psi(0.7, theta), 1.0 / (1.0 + 0.49), 1e-8);
    }

    #[test]
    fn dilute_limit() {
        // β → ∞ で J → (1/β)∫ψ dx = π/(2β)、K → 2J
        for theta in [0.01, 0.3, 10.0] {
            let [(j, k)] = jk_integrals(&[1.0e8], theta)[..] else {
                panic!("one β")
            };
            assert_rel(j, PI / 2.0e8, 1e-7);
            assert_rel(k, PI / 1.0e8, 1e-7);
        }
    }

    /// J を x について Simpson 則で直接積分した値（x ≤ X、X より先は ψ ≈ 1/x² で足す）。
    fn j_by_simpson(beta: f64, theta: f64) -> f64 {
        let x_max = 1.0e4;
        let n = 2_000_000;
        let h = x_max / n as f64;
        let f = |x: f64| {
            let p = psi(x, theta);
            p / (p + beta)
        };
        let mut sum = f(0.0) + f(x_max);
        for i in 1..n {
            sum += if i % 2 == 1 { 4.0 } else { 2.0 } * f(h * i as f64);
        }
        sum * h / 3.0 + 1.0 / (beta * x_max)
    }

    #[test]
    fn doppler_dominated_line_matches_direct_integration() {
        // ドップラー幅が自然幅の 100 倍
        let theta = 0.01;
        let betas = [1.0e-3, 0.1];
        for (&beta, (j, _)) in betas.iter().zip(jk_integrals(&betas, theta)) {
            assert_rel(j, j_by_simpson(beta, theta), 1e-6);
        }
    }
}
