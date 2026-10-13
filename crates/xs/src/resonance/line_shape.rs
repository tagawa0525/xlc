//! ドップラー拡がりした共鳴の線形状と、ETOX 法の J, K 積分（NJOY UNRESR の ajku に相当）。
//!
//! x = 2(E - E_r)/Γ、θ = Γ/Δ（Δ はドップラー幅）として、線形状は
//! ψ(x, θ) = √π (θ/2) Re w(θx/2 + iθ/2)（θ → ∞ で 1/(1 + x²)）。
//! J(β, θ) = ∫₀^∞ ψ/(ψ + β) dx、K(β, θ) = ∫₀^∞ ψ(ψ + 2β)/(ψ + β)² dx。

use std::sync::atomic::{AtomicUsize, Ordering};

use num_complex::Complex64;

use super::quadrature::gauss_kronrod;
use super::unresolved::UrrError;
use crate::numeric::NumericError;
use crate::special::faddeeva;

/// J, K 積分の相対許容誤差（全成分の最大）。外側の幅のゆらぎの積分より小さくする
const TOLERANCE: f64 = 1.0e-12;

/// 線形状 ψ(x, θ)。
pub(super) fn psi(x: f64, theta: f64) -> f64 {
    let y = 0.5 * theta;
    std::f64::consts::PI.sqrt() * y * faddeeva(Complex64::new(y * x, y)).re
}

/// 各 β について (J(β, θ), K(β, θ)) を求める。β は正であること。
///
/// x = tan φ と置いて [0, π/2] を適応 Gauss-Kronrod 積分する。線形状の裾は ψ ~ 1/x² なので
/// ψ/(ψ+β) dx = ψ/(ψ+β) (1 + x²) dφ は φ → π/2 でも有限（1/β に近づく）。
pub(super) fn jk_integrals(betas: &[f64], theta: f64) -> Result<Vec<(f64, f64)>, NumericError> {
    if let Some(beta) = betas.iter().find(|&&b| b <= 0.0 || b.is_nan()) {
        return Err(NumericError::Domain(format!(
            "J, K integrals for β = {beta}"
        )));
    }
    let integrand = |phi: f64, out: &mut [f64]| {
        let x = phi.tan();
        let p = psi(x, theta);
        let jacobian = 1.0 + x * x;
        for (i, &beta) in betas.iter().enumerate() {
            let r = p / (p + beta);
            out[2 * i] = r * jacobian;
            out[2 * i + 1] = r * (p + 2.0 * beta) / (p + beta) * jacobian;
        }
        Ok(())
    };
    let values = gauss_kronrod::<NumericError>(
        integrand,
        2 * betas.len(),
        0.0,
        0.5 * std::f64::consts::PI,
        TOLERANCE,
    )?;
    Ok(values.chunks(2).map(|c| (c[0], c[1])).collect())
}

/// J, K の表の範囲（log10 β, log10 θ）。燃料の核種の非分離共鳴で実際に使う範囲
/// （β ≳ 1e-2、θ ≈ 0.04〜500）を十分含む。構造材は低温で θ が大きくなる（Fe-58 の 550 K で θ ≈ 1.0e4）。
/// θ の上限は、その上で使う Lorentz の線形状との差（1/θ² の程度）が 1e-12 になる 10⁶ とする（ADR-0107）
const LOG_BETA: (f64, f64) = (-6.0, 10.0);
const LOG_THETA: (f64, f64) = (-3.0, 6.0);

/// 表の点の間隔（log10）
const STEP: f64 = 0.025;

/// ln J, ln K を (log10 β, log10 θ) の等間隔の格子で持つ表。
struct JkTable {
    n_beta: usize,
    ln_j: Vec<f64>,
    ln_k: Vec<f64>,
}

/// 表は最初に使うときに一度だけ作る（作れなかったときのエラーも保持する）。
fn table() -> Result<&'static JkTable, NumericError> {
    static TABLE: std::sync::OnceLock<Result<JkTable, NumericError>> = std::sync::OnceLock::new();
    let table = TABLE.get_or_init(|| {
        let points = |(lo, hi): (f64, f64)| ((hi - lo) / STEP).round() as usize + 1;
        let (n_beta, n_theta) = (points(LOG_BETA), points(LOG_THETA));
        let betas: Vec<f64> = (0..n_beta)
            .map(|i| 10f64.powf(LOG_BETA.0 + STEP * i as f64))
            .collect();
        // θ の行は互いに独立に積分するので、並列に作っても各点の値はスレッド数によらない（行どうしで和を取らない）。
        // 行のエラーは θ の小さい順に最初のものを返す
        let rows = in_parallel(n_theta, |t| {
            jk_integrals(&betas, 10f64.powf(LOG_THETA.0 + STEP * t as f64))
        })
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?;
        let (ln_j, ln_k) = rows
            .iter()
            .flatten()
            .map(|&(j, k)| (j.ln(), k.ln()))
            .unzip();
        Ok(JkTable { n_beta, ln_j, ln_k })
    });
    table.as_ref().map_err(Clone::clone)
}

/// `f(i)`（i = 0..n）を、使える CPU の数のスレッドで並列に求め、i の順に並べて返す。
/// 行の重さは θ で違うので、各スレッドは次の i を共有のカウンタから取る。
/// rayon を使わないのは、`OnceLock` の初期化の中で rayon のプールを待つと止まりうるため（#0088）。
fn in_parallel<T: Send>(n: usize, f: impl Fn(usize) -> T + Sync) -> Vec<T> {
    let threads = std::thread::available_parallelism()
        .map_or(1, std::num::NonZeroUsize::get)
        .min(n);
    let next = AtomicUsize::new(0);
    let mut done: Vec<(usize, T)> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                scope.spawn(|| {
                    let mut done = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= n {
                            break done;
                        }
                        done.push((i, f(i)));
                    }
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
            .collect()
    });
    done.sort_unstable_by_key(|&(i, _)| i);
    done.into_iter().map(|(_, value)| value).collect()
}

/// 表から (J(β, θ), K(β, θ)) を補間する（ln J, ln K を両方向に 4 点の Lagrange 補間）。
/// β が表より大きければ希釈の極限 J = π/(2β)、K = π/β を返す（誤差は ψ の最大値/β ≤ 1e-10）。
/// θ が表より大きければ Lorentz の線形状の J, K を返す（`natural_jk`）。
pub(super) fn jk_interpolated(beta: f64, theta: f64) -> Result<(f64, f64), UrrError> {
    let (lb, lt) = (beta.log10(), theta.log10());
    if lb > LOG_BETA.1 {
        let j = 0.5 * std::f64::consts::PI / beta;
        return Ok((j, 2.0 * j));
    }
    if lt > LOG_THETA.1 {
        return Ok(natural_jk(beta));
    }
    if !(LOG_BETA.0..=LOG_BETA.1).contains(&lb) || !(LOG_THETA.0..=LOG_THETA.1).contains(&lt) {
        return Err(UrrError::JkTable { beta, theta });
    }
    let tab = table()?;
    let n_theta = tab.ln_j.len() / tab.n_beta;
    let (ib, wb) = lagrange_weights((lb - LOG_BETA.0) / STEP, tab.n_beta);
    let (it, wt) = lagrange_weights((lt - LOG_THETA.0) / STEP, n_theta);
    let (mut ln_j, mut ln_k) = (0.0, 0.0);
    for (a, &w_t) in wt.iter().enumerate() {
        for (b, &w_b) in wb.iter().enumerate() {
            let index = (it + a) * tab.n_beta + ib + b;
            ln_j += w_t * w_b * tab.ln_j[index];
            ln_k += w_t * w_b * tab.ln_k[index];
        }
    }
    Ok((ln_j.exp(), ln_k.exp()))
}

/// ドップラー拡がりのない Lorentz の線形状 ψ = 1/(1 + x²) の (J, K)。a² = (1 + β)/β として
/// J = π/(2aβ)、K = 2J − ∫₀^∞ ψ²/(ψ + β)² dx = 2J − π/(4a³β²)。θ が表より大きいとき（自然幅がドップラー幅より
/// 十分大きい）に使う。Voigt の線形状との差は 1/θ² の程度（ADR-0107）。β > 0 ならどこでも成り立つ
fn natural_jk(beta: f64) -> (f64, f64) {
    let pi = std::f64::consts::PI;
    let a = ((1.0 + beta) / beta).sqrt();
    let j = pi / (2.0 * a * beta);
    (j, 2.0 * j - pi / (4.0 * a.powi(3) * beta * beta))
}

/// 格子座標 u（0 から n-1）での 4 点の Lagrange 補間の、最初の点の番号と重み。
fn lagrange_weights(u: f64, n: usize) -> (usize, [f64; 4]) {
    let first = (u.floor() as usize).saturating_sub(1).min(n - 4);
    let t = u - first as f64;
    (
        first,
        [
            -(t - 1.0) * (t - 2.0) * (t - 3.0) / 6.0,
            t * (t - 2.0) * (t - 3.0) / 2.0,
            -t * (t - 1.0) * (t - 3.0) / 2.0,
            t * (t - 1.0) * (t - 2.0) / 6.0,
        ],
    )
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
    fn jk_integrals_need_positive_beta() {
        for beta in [0.0, -1.0, f64::NAN] {
            assert!(matches!(
                jk_integrals(&[1.0, beta], 1.0),
                Err(NumericError::Domain(_))
            ));
        }
    }

    #[test]
    fn natural_line_shape_without_doppler() {
        // θ → ∞ で ψ = 1/(1+x²)。a² = (1+β)/β として
        // J = π/(2aβ)、K = (1/β)[π/(2a) + π(a²+1)/(4a³)]
        let theta = 1.0e5;
        let betas = [1.0e-4, 1.0e-2, 1.0, 1.0e2];
        for (&beta, (j, k)) in betas.iter().zip(jk_integrals(&betas, theta).unwrap()) {
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
            let [(j, k)] = jk_integrals(&[1.0e8], theta).unwrap()[..] else {
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
        for (&beta, (j, _)) in betas.iter().zip(jk_integrals(&betas, theta).unwrap()) {
            assert_rel(j, j_by_simpson(beta, theta), 1e-6);
        }
    }

    #[test]
    fn table_matches_direct_integration() {
        // 表の範囲から決まった疑似乱数で選んだ点で、補間と直接の積分の相対差が 1e-6 以下
        let mut seed: u64 = 12345;
        let mut uniform = || {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        for _ in 0..400 {
            let beta = 10f64.powf(LOG_BETA.0 + (LOG_BETA.1 - LOG_BETA.0) * uniform());
            let theta = 10f64.powf(LOG_THETA.0 + (LOG_THETA.1 - LOG_THETA.0) * uniform());
            let (j, k) = jk_interpolated(beta, theta).unwrap();
            let [(jd, kd)] = jk_integrals(&[beta], theta).unwrap()[..] else {
                panic!("one β")
            };
            assert_rel(j, jd, 1e-6);
            assert_rel(k, kd, 1e-6);
        }
    }

    #[test]
    fn table_uses_dilute_limit_above_its_range() {
        let (j, k) = jk_interpolated(1.0e12, 0.5).unwrap();
        assert_rel(j, PI / 2.0e12, 1e-12);
        assert_rel(k, PI / 1.0e12, 1e-12);
    }

    #[test]
    fn table_uses_natural_line_shape_above_its_theta_range() {
        // θ = Γ/Δ が表より大きいとき（低温の構造材の非分離共鳴。Fe-58 の 550 K で θ ≈ 1.0e4）は、ドップラー
        // 拡がりのない Lorentz の線形状の J, K。表の端で値が跳ぶと、幅のゆらぎの適応積分（許容誤差 1e-10）が
        // 収束しにくいので、表の端の内と外で直接の積分との差が 1e-10 より小さいこと
        let top = 10f64.powf(LOG_THETA.1);
        for theta in [1.0025e4, 0.999 * top, 1.001 * top, 1.0e8] {
            for beta in [1.0e-6, 1.0e-4, 1.0e-2, 1.0, 1.0e2, 1.0e10] {
                let (j, k) = jk_interpolated(beta, theta).unwrap();
                let direct = jk_integrals(&[beta], theta).unwrap()[0];
                assert_rel(j, direct.0, 1e-10);
                assert_rel(k, direct.1, 1e-10);
            }
        }
        // 表の β より小さくても、Lorentz の線形状の J, K は使える
        let (j, k) = jk_interpolated(1.0e-8, 1.0e8).unwrap();
        let direct = jk_integrals(&[1.0e-8], 1.0e8).unwrap()[0];
        assert_rel(j, direct.0, 1e-10);
        assert_rel(k, direct.1, 1e-10);
    }

    #[test]
    fn table_rejects_points_outside_its_range() {
        assert!(jk_interpolated(1.0e-8, 1.0).is_err());
        assert!(jk_interpolated(1.0, 1.0e-5).is_err());
    }
}
