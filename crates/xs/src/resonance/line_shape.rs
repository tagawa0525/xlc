//! ドップラー拡がりした共鳴の線形状と、ETOX 法の J, K 積分（NJOY UNRESR の ajku に相当）。
//!
//! x = 2(E - E_r)/Γ、θ = Γ/Δ（Δ はドップラー幅）として、線形状は
//! ψ(x, θ) = √π (θ/2) Re w(θx/2 + iθ/2)（θ → ∞ で 1/(1 + x²)）。
//! J(β, θ) = ∫₀^∞ ψ/(ψ + β) dx、K(β, θ) = ∫₀^∞ ψ(ψ + 2β)/(ψ + β)² dx。

use crate::special::faddeeva;

/// J, K 積分の相対許容誤差（全成分の最大）
const TOLERANCE: f64 = 1.0e-10;

/// 適応積分の区間の数の上限
const MAX_INTERVALS: usize = 10_000;

/// 線形状 ψ(x, θ)。
pub(super) fn psi(x: f64, theta: f64) -> f64 {
    let y = 0.5 * theta;
    std::f64::consts::PI.sqrt() * y * faddeeva(y * x, y).0
}

/// 各 β について (J(β, θ), K(β, θ)) を求める。β は正であること。
///
/// x = tan φ と置いて [0, π/2] を適応 Gauss-Kronrod 積分する。線形状の裾は ψ ~ 1/x² なので
/// ψ/(ψ+β) dx = ψ/(ψ+β) (1 + x²) dφ は φ → π/2 でも有限（1/β に近づく）。
pub(super) fn jk_integrals(betas: &[f64], theta: f64) -> Vec<(f64, f64)> {
    assert!(
        betas.iter().all(|&b| b > 0.0),
        "β must be positive: {betas:?}"
    );
    let integrand = |phi: f64, out: &mut [f64]| {
        let x = phi.tan();
        let p = psi(x, theta);
        let jacobian = 1.0 + x * x;
        for (i, &beta) in betas.iter().enumerate() {
            let r = p / (p + beta);
            out[2 * i] = r * jacobian;
            out[2 * i + 1] = r * (p + 2.0 * beta) / (p + beta) * jacobian;
        }
    };
    let values = gauss_kronrod(integrand, 2 * betas.len(), 0.0, 0.5 * std::f64::consts::PI);
    values.chunks(2).map(|c| (c[0], c[1])).collect()
}

/// [a, b] でベクトル値の関数を適応 Gauss-Kronrod（7-15 点）積分する。誤差の見積もりが最も大きい
/// 区間を順に二分し、全成分で誤差の合計が積分値の `TOLERANCE` 倍以下になったら止める。
fn gauss_kronrod(f: impl Fn(f64, &mut [f64]), n: usize, a: f64, b: f64) -> Vec<f64> {
    let mut intervals = vec![kronrod_interval(&f, n, a, b)];
    loop {
        let mut total = vec![0.0; n];
        let mut error = vec![0.0; n];
        for iv in &intervals {
            for i in 0..n {
                total[i] += iv.value[i];
                error[i] += iv.error[i];
            }
        }
        let converged = (0..n).all(|i| error[i] <= TOLERANCE * total[i].abs());
        if converged {
            return total;
        }
        assert!(
            intervals.len() < MAX_INTERVALS,
            "Gauss-Kronrod integration did not converge"
        );
        // 相対誤差が最大の区間を分ける
        let worst = (0..intervals.len())
            .max_by(|&p, &q| {
                let score = |iv: &Interval| {
                    (0..n)
                        .map(|i| iv.error[i] / total[i].abs().max(f64::MIN_POSITIVE))
                        .fold(0.0, f64::max)
                };
                score(&intervals[p]).total_cmp(&score(&intervals[q]))
            })
            .expect("at least one interval");
        let iv = intervals.swap_remove(worst);
        let mid = 0.5 * (iv.a + iv.b);
        intervals.push(kronrod_interval(&f, n, iv.a, mid));
        intervals.push(kronrod_interval(&f, n, mid, iv.b));
    }
}

/// 1区間の Kronrod 15 点の値と、Gauss 7 点との差による誤差の見積もり。
struct Interval {
    a: f64,
    b: f64,
    value: Vec<f64>,
    error: Vec<f64>,
}

fn kronrod_interval(f: &impl Fn(f64, &mut [f64]), n: usize, a: f64, b: f64) -> Interval {
    // Kronrod の節点（正の側）と重み、Gauss 7 点の重み（節点は奇数番目と中心）
    const XK: [f64; 8] = [
        0.991_455_371_120_812_6,
        0.949_107_912_342_758_5,
        0.864_864_423_359_769_1,
        0.741_531_185_599_394_4,
        0.586_087_235_467_691_1,
        0.405_845_151_377_397_2,
        0.207_784_955_007_898_5,
        0.0,
    ];
    const WK: [f64; 8] = [
        0.022_935_322_010_529_22,
        0.063_092_092_629_978_55,
        0.104_790_010_322_250_2,
        0.140_653_259_715_525_9,
        0.169_004_726_639_267_9,
        0.190_350_578_064_785_4,
        0.204_432_940_075_298_9,
        0.209_482_141_084_727_8,
    ];
    const WG: [f64; 4] = [
        0.129_484_966_168_869_7,
        0.279_705_391_489_276_7,
        0.381_830_050_505_118_9,
        0.417_959_183_673_469_4,
    ];
    let (center, half) = (0.5 * (a + b), 0.5 * (b - a));
    let mut kronrod = vec![0.0; n];
    let mut gauss = vec![0.0; n];
    let mut buffer = vec![0.0; n];
    let mut add = |x: f64, wk: f64, wg: f64, kronrod: &mut [f64], gauss: &mut [f64]| {
        f(x, &mut buffer);
        for i in 0..n {
            kronrod[i] += wk * buffer[i];
            gauss[i] += wg * buffer[i];
        }
    };
    for (j, (&x, &wk)) in XK.iter().zip(&WK).enumerate() {
        // 奇数番目（j = 1, 3, 5）と中心（j = 7）は Gauss の節点でもある
        let wg = if j % 2 == 1 { WG[j / 2] } else { 0.0 };
        if x == 0.0 {
            add(center, wk, wg, &mut kronrod, &mut gauss);
        } else {
            add(center - half * x, wk, wg, &mut kronrod, &mut gauss);
            add(center + half * x, wk, wg, &mut kronrod, &mut gauss);
        }
    }
    Interval {
        a,
        b,
        value: kronrod.iter().map(|v| v * half).collect(),
        error: kronrod
            .iter()
            .zip(&gauss)
            .map(|(k, g)| ((k - g) * half).abs())
            .collect(),
    }
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
