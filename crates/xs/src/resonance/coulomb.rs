//! 荷電粒子のチャネルの Coulomb の透過因子 P_l(η, ρ) = ρ / (F_l² + G_l²)（ENDF-102 の式 D.84）。
//!
//! H = G + iF の対数微分を Y = H'/H とすると、Wronskian F'G − FG' = 1 から ρY = S + iP（S はシフト因子）。
//! 転回点 ρ_t = η + √(η² + l(l+1)) より外では、Y を Steed の連分数 CF2（DLMF 33.8.2）で直接求める。
//! 障壁の下では CF2 の収束が遅いので、ρ₀ = max(ρ_t, `RHO_MIN`) で CF2 から Y を求め、Coulomb の方程式
//! u'' = (2η/ρ + l(l+1)/ρ² − 1) u から導いた Riccati の式を x = ln ρ で ρ まで内向きに積分する:
//!
//!   dS/dx = S − S² + P² + l(l+1) + 2ηρ − ρ²,   d(ln P)/dx = 1 − 2S
//!
//! P は障壁の下で e^{−2πη} の桁まで小さくなるので、P でなく ln P を積分し、最後に指数を取る（アンダーフローは 0 になる）。
//! 積分は Taylor 級数の方法（係数を漸化式で求め、末尾の係数から刻みを決める）で行う。方法の比較は ADR-0273。

use num_complex::Complex64;

/// CF2 を始める ρ の下限。転回点の小さい（η と l の小さい）チャネルで CF2 の収束を速くする
/// （l = 0、η = 0.01 で、ρ = 0.02 では 2515 項、ρ = 2 では 52 項）。
const RHO_MIN: f64 = 2.0;
/// Taylor 級数の次数
const ORDER: usize = 20;
/// ln P がこれより小さければ P は倍精度で 0（最小の正の非正規化数は e^{−744.4}）
const LN_UNDERFLOW: f64 = -746.0;
/// CF2 の項の数の上限（l ≤ 4、η ≤ 1e8 で 2624 項）
const MAX_TERMS: usize = 100_000;

/// Steed の CF2 が `MAX_TERMS` 項で収束しなかった（ρ₀ ≥ max(ρ_t, `RHO_MIN`) では起きない見込みで、NaN の入力などで起きる）。
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
#[error(
    "Steed's CF2 for the Coulomb functions with l={l} eta={eta} rho={rho} did not converge in {MAX_TERMS} terms"
)]
pub struct NotConverged {
    pub l: u32,
    pub eta: f64,
    pub rho: f64,
}

/// 軌道角運動量 l、Coulomb のパラメータ η ≥ 0、ρ > 0 の透過因子 P_l(η, ρ)。
pub(super) fn penetrability(l: u32, eta: f64, rho: f64) -> Result<f64, NotConverged> {
    let ll = f64::from(l * (l + 1));
    let turning = eta + (eta * eta + ll).sqrt();
    let start = rho.max(turning).max(RHO_MIN);
    let y = cf2(l, eta, start)?;
    if start == rho {
        return Ok(rho * y.im);
    }
    let mut s = start * y.re;
    let mut p = (start * y.im).ln();
    let mut x = start.ln();
    let (end, x_turning) = (rho.ln(), turning.ln());
    while x > end {
        let r = x.exp();
        let (sc, pc) = taylor(s, p, ll, eta, r);
        let h = -step(&sc).min(step(&pc)).min(x - end);
        s = horner(&sc, h);
        p = horner(&pc, h);
        x = if x + h <= end { end } else { x + h };
        // 障壁の下（ρ < ρ_t）で S < 0 なら、内向きに S < 0 のまま（S = 0 では dS/dx > 0）で、d(ln P)/dx = 1 − 2S > 0 から
        // ln P は ρ まで減り続ける。ここで倍精度の 0 になれば、ρ でも 0
        if x < x_turning && s < 0.0 && p < LN_UNDERFLOW {
            return Ok(0.0);
        }
    }
    Ok(p.exp())
}

/// H'/H = i(1 − η/ρ) + (i/ρ) ab / (2(ρ − η + i) + (a+1)(b+1) / (2(ρ − η + 2i) + ...))、a = 1 + l + iη、b = −l + iη
/// （DLMF 33.8.2）。分母の連分数を修正 Lentz 法で求める（先頭の項を分子に置くと、η の大きいときに初期値で
/// 桁があふれる）。
fn cf2(l: u32, eta: f64, rho: f64) -> Result<Complex64, NotConverged> {
    let a = Complex64::new(1.0 + f64::from(l), eta);
    let b = Complex64::new(-f64::from(l), eta);
    let tiny = Complex64::new(f64::MIN_POSITIVE, 0.0);
    let mut g = Complex64::new(2.0 * (rho - eta), 2.0);
    let (mut c, mut d) = (g, Complex64::ZERO);
    for n in 2..MAX_TERMS {
        let k = n as f64;
        let an = (a + k - 1.0) * (b + k - 1.0);
        let bn = Complex64::new(2.0 * (rho - eta), 2.0 * k);
        d = bn + an * d;
        if d == Complex64::ZERO {
            d = tiny;
        }
        c = bn + an / c;
        if c == Complex64::ZERO {
            c = tiny;
        }
        d = d.inv();
        let delta = c * d;
        g *= delta;
        if (delta - 1.0).norm() < 2.0 * f64::EPSILON {
            return Ok(Complex64::new(0.0, 1.0 - eta / rho) + Complex64::I * a * b / (rho * g));
        }
    }
    Err(NotConverged { l, eta, rho })
}

/// x = ln r での S と ln P の Taylor 級数の係数（0 から `ORDER` 次）。P² の係数も同じ漸化式で求める
/// （d(P²)/dx = 2P²(1 − 2S)）。ρ = r eᵗ、ρ² = r² e²ᵗ の係数は r/n!、r² 2ⁿ/n!。
fn taylor(s0: f64, p0: f64, ll: f64, eta: f64, r: f64) -> ([f64; ORDER + 1], [f64; ORDER + 1]) {
    let mut s = [0.0; ORDER + 1];
    let mut p = [0.0; ORDER + 1];
    let mut q = [0.0; ORDER + 1];
    s[0] = s0;
    p[0] = p0;
    q[0] = (2.0 * p0).exp();
    // ρ と ρ² の n 次の係数
    let (mut rho_n, mut rho2_n) = (r, r * r);
    for n in 0..ORDER {
        let m = n as f64 + 1.0;
        let ss: f64 = (0..=n).map(|j| s[j] * s[n - j]).sum();
        let sq: f64 = (0..=n).map(|j| s[j] * q[n - j]).sum();
        let constant = if n == 0 { ll } else { 0.0 };
        s[n + 1] = (s[n] - ss + q[n] + constant + 2.0 * eta * rho_n - rho2_n) / m;
        q[n + 1] = (2.0 * q[n] - 4.0 * sq) / m;
        p[n + 1] = (if n == 0 { 1.0 } else { 0.0 } - 2.0 * s[n]) / m;
        rho_n /= m;
        rho2_n *= 2.0 / m;
    }
    (s, p)
}

/// 末尾の 2 つの係数の大きさが、値の大きさ（1 以上）の丸めの誤差になる刻み（Jorba と Zou の刻みの決め方）。
fn step(c: &[f64; ORDER + 1]) -> f64 {
    let scale = c[0].abs().max(1.0) * f64::EPSILON;
    let radius = |k: usize| {
        let a = c[k].abs();
        if a == 0.0 {
            f64::INFINITY
        } else {
            ((scale.ln() - a.ln()) / k as f64).exp()
        }
    };
    0.9 * radius(ORDER).min(radius(ORDER - 1))
}

fn horner(c: &[f64; ORDER + 1], h: f64) -> f64 {
    c.iter().rev().fold(0.0, |acc, &v| acc * h + v)
}

#[cfg(test)]
mod tests;
