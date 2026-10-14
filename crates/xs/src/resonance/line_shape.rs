//! ドップラー拡がりした共鳴の線形状と、ETOX 法の J, K 積分（NJOY UNRESR の ajku に相当）。
//!
//! x = 2(E - E_r)/Γ、θ = Γ/Δ（Δ はドップラー幅）として、線形状は
//! ψ(x, θ) = √π (θ/2) Re w(θx/2 + iθ/2)（θ → ∞ で 1/(1 + x²)）。
//! J(β, θ) = ∫₀^∞ ψ/(ψ + β) dx、K(β, θ) = ∫₀^∞ ψ(ψ + 2β)/(ψ + β)² dx。

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

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

/// これより小さい θ では、積分の変数を線形状の芯の幅 2/θ で伸ばす（`jk_integrals`）
const SCALED_BELOW: f64 = 1.0e-3;

/// 各 β について (J(β, θ), K(β, θ)) を求める。β は正であること。
///
/// x = s tan φ と置いて [0, π/2] を適応 Gauss-Kronrod 積分する（s は θ ≥ 1e-3 で 1、それより小さい θ で 2/θ）。線形状の裾は
/// ψ ~ 1/x² なので、ψ/(ψ+β) dx = ψ/(ψ+β) s (1 + tan²φ) dφ は φ → π/2 でも有限。
pub(super) fn jk_integrals(betas: &[f64], theta: f64) -> Result<Vec<(f64, f64)>, NumericError> {
    if let Some(beta) = betas.iter().find(|&&b| b <= 0.0 || b.is_nan()) {
        return Err(NumericError::Domain(format!(
            "J, K integrals for β = {beta}"
        )));
    }
    // θ が小さいと線形状の芯はドップラー幅の x ≈ 2/θ まで広がり、x = tan φ では φ = π/2 の近くの幅 θ ほどに潰れて
    // 適応積分が収束しない。そのときは x = (2/θ) tan φ と芯の幅で伸ばす。θ ≥ 1e-3 は x = tan φ で収束する（ADR-0170）
    let scale = if theta < SCALED_BELOW {
        2.0 / theta
    } else {
        1.0
    };
    let integrand = |phi: f64, out: &mut [f64]| {
        let t = phi.tan();
        let x = scale * t;
        let p = psi(x, theta);
        let jacobian = scale * (1.0 + t * t);
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

/// J, K の表の範囲（log10 β, log10 θ）。燃料の主な核種の非分離共鳴（β ≳ 1e-2、θ ≈ 0.04〜500）と、自然幅がドップラー幅
/// より狭い核分裂生成物など（Nd-148 の 900 K で θ ≈ 9.8e-4、Sm-154 の 1200 K で 1.4e-4。#0119）を含む。構造材は低温で θ が
/// 大きくなる（Fe-58 の 550 K で θ ≈ 1.0e4）。θ の上限は、その上で使う Lorentz の線形状との差（1/θ² の程度）が 1e-12 になる
/// 10⁶（ADR-0107）、下限は、その下で使うガウス形の線形状との差（θ の程度）が 1e-10 になる 10⁻¹⁰
const LOG_BETA: (f64, f64) = (-6.0, 10.0);
const LOG_THETA: (f64, f64) = (-10.0, 6.0);

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
        // θ の行は互いに独立に積分するので、並列に作っても各点の値はスレッド数によらない（行どうしで和を取らない）
        let rows = in_parallel(n_theta, |t| {
            jk_integrals(&betas, 10f64.powf(LOG_THETA.0 + STEP * t as f64))
        })?;
        let (ln_j, ln_k) = rows
            .iter()
            .flatten()
            .map(|&(j, k)| (j.ln(), k.ln()))
            .unzip();
        Ok(JkTable { n_beta, ln_j, ln_k })
    });
    table.as_ref().map_err(Clone::clone)
}

/// `f(i)`（i = 0..n）を、物理コアの数（ADR-0067。このプロセスに許された CPU の数で抑える）のスレッドで並列に求め、
/// i の順に並べて返す。
/// 行の重さは θ で違うので、各スレッドは次の i を共有のカウンタから取る。
/// rayon を使わないのは、`OnceLock` の初期化の中で rayon のプールを待つと止まりうるため（#0088）。
///
/// エラーが出たら、各スレッドは次に i を取る前にそれを見て止まる（見たあとに取った i は、スレッドごとに多くとも 1 つ
/// 余分に求める）。返すのは i の最も小さいエラーで、カウンタは増える一方なので、エラーの i より小さい i はすべて
/// 取られて最後まで求められており、返すエラーはスレッドの進み方によらない。
fn in_parallel<T: Send, E: Send>(
    n: usize,
    f: impl Fn(usize) -> Result<T, E> + Sync,
) -> Result<Vec<T>, E> {
    let allowed = std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get);
    let threads = num_cpus::get_physical().min(allowed).min(n);
    let next = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    let mut done: Vec<(usize, Result<T, E>)> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                scope.spawn(|| {
                    let mut done = Vec::new();
                    while !failed.load(Ordering::Relaxed) {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= n {
                            break;
                        }
                        let value = f(i);
                        if value.is_err() {
                            failed.store(true, Ordering::Relaxed);
                        }
                        done.push((i, value));
                    }
                    done
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
/// θ が表より大きければ Lorentz の線形状の J, K を返す（`natural_jk`）。θ が表より小さければ（β は表の範囲）、ガウス形の
/// 線形状の J, K を返す（`gaussian_jk`）。
pub(super) fn jk_interpolated(beta: f64, theta: f64) -> Result<(f64, f64), UrrError> {
    if theta.is_nan() || theta <= 0.0 {
        return Err(UrrError::JkTable { beta, theta });
    }
    let (lb, lt) = (beta.log10(), theta.log10());
    if lb > LOG_BETA.1 {
        let j = 0.5 * std::f64::consts::PI / beta;
        return Ok((j, 2.0 * j));
    }
    if lt > LOG_THETA.1 {
        return Ok(natural_jk(beta));
    }
    if lt < LOG_THETA.0 && lb >= LOG_BETA.0 {
        return Ok(gaussian_jk(beta, theta));
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

/// ドップラー幅に比べて自然幅が無視できる線形状 ψ = (√π θ/2) e^{-u²}（u = θx/2）の (J, K)。θ が表より小さいときに使う。
/// b = 2β/(√π θ)、g = e^{-u²} として J = (2/θ) ∫₀^∞ g/(g + b) du、K = (2/θ) ∫₀^∞ g(g + 2b)/(g + b)² du を 1/b の冪に展開し、
/// ∫₀^∞ gⁿ du = √π/(2√n) で項ごとに積分する:
/// J = (√π/θ) Σₙ₌₁ (−1)ⁿ⁺¹ b⁻ⁿ/√n、K = (√π/θ) Σₙ₌₁ (−1)ⁿ⁺¹ (n + 1) b⁻ⁿ/√n。
/// 表の下限（θ = 10⁻¹⁰）より下で β ≥ 10⁻⁶ なら b > 10⁴ なので、級数は数項で倍精度に収束する。
/// Voigt の線形状のローレンツの裾による差は相対で θ の程度（表の下限で 1e-10 以下。`quadrature_approaches_the_gaussian_limit_for_narrow_lines`）
fn gaussian_jk(beta: f64, theta: f64) -> (f64, f64) {
    // √π/θ と b⁻ⁿ を別に作ると、θ が小さいと b が倍精度を超えて 0·∞ になるので、先頭の項 √π/(θ b) = π/(2β) でくくる:
    // J = π/(2β) Σₙ₌₁ (−1)ⁿ⁺¹ b¹⁻ⁿ/√n
    let inverse_b = std::f64::consts::PI.sqrt() * theta / (2.0 * beta);
    let (mut j, mut k) = (0.0, 0.0);
    let mut power = 1.0;
    for n in 1..=SERIES_TERMS {
        let nf = f64::from(n);
        let sign = if n % 2 == 1 { 1.0 } else { -1.0 };
        j += sign * power / nf.sqrt();
        k += sign * (nf + 1.0) * power / nf.sqrt();
        power *= inverse_b;
    }
    let leading = 0.5 * std::f64::consts::PI / beta;
    (leading * j, leading * k)
}

/// `gaussian_jk` の級数の項の数。b > 10⁴ では、6 項目は先頭の項に比べて相対 4e-21（J）、1.4e-20（K）より小さく、交代級数なので打ち切りの誤差はそれ以下
const SERIES_TERMS: u32 = 5;

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
mod tests;
