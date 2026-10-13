//! 自由ガス模型による熱中性子の弾性散乱。
//!
//! 標的核を温度 T の理想気体（Maxwell 分布で運動する自由な原子）とみなし、実験室系の
//! 散乱を等方として、出射エネルギーの分布 S(E→E') を解析的に求める。
//!
//! S は自由原子の断面積を 1 としたときの値で、規格化しない。∫S dE' = N(E) は自由ガスの実効断面積と
//! 自由原子の断面積の比で、低エネルギーほど 1 より大きい。散乱行列では、温度に拡げた断面積
//! （一定の自由原子の断面積なら σ N(E)）に S/N(E) を掛ける。自由原子の断面積が一定なら σ S に等しく、
//! 詳細釣り合いが厳密に成り立つ。

use crate::special::{erf, erfcx};
use xlc_format::group::GroupStructure;
use xlc_numeric::quadrature::{GAUSS_LEGENDRE_4, GAUSS_LEGENDRE_8, legendre};

/// 熱領域の弾性散乱を自由ガス模型で扱う条件。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FreeGas {
    /// kT (eV)
    pub kt: f64,
    /// この入射エネルギー (eV) 未満で自由ガス模型を使う
    pub cutoff: f64,
}

/// 自由原子の散乱断面積を 1 としたときの、エネルギー移行の確率密度 S(E→E') (1/eV)。
/// `awr` は標的の質量比。
pub fn free_gas_kernel(e: f64, e_out: f64, awr: f64, kt: f64) -> f64 {
    if e <= 0.0 || e_out <= 0.0 {
        return 0.0;
    }
    let sqrt_a = awr.sqrt();
    let theta = (awr + 1.0) / (2.0 * sqrt_a);
    let zeta = (awr - 1.0) / (2.0 * sqrt_a);
    let (s, s_out) = ((e / kt).sqrt(), (e_out / kt).sqrt());
    // e^(-y²) erfcx(x) = e^(-y²+x²) erfc(x)。指数を先にまとめて桁あふれを避ける
    let scaled = |y: f64, x: f64| (-y * y).exp() * erfcx(x);
    let bracket = if e_out <= e {
        // erf(θ√ε' - ζ√ε) + erf(θ√ε' + ζ√ε) + e^(ε-ε') [erfc(θ√ε + ζ√ε') - erfc(θ√ε - ζ√ε')]
        let (x1, x2) = (theta * s + zeta * s_out, theta * s - zeta * s_out);
        let (u, v) = (zeta * s - theta * s_out, zeta * s + theta * s_out);
        if u > 0.0 {
            // 遠い下方散乱: O(1) の項が打ち消し合うので、同じ指数因子をもつ項どうしをまとめる
            (-u * u).exp() * (erfcx(u) - erfcx(x2)) + (-v * v).exp() * (erfcx(x1) - erfcx(v))
        } else {
            erf(-u) + erf(v) + scaled(v, x1) - scaled(u, x2)
        }
    } else {
        // erf(θ√ε' - ζ√ε) - erf(θ√ε' + ζ√ε) + e^(ε-ε') [erf(θ√ε - ζ√ε') + erf(θ√ε + ζ√ε')]
        let (p, q) = (theta * s_out - zeta * s, theta * s_out + zeta * s);
        let (r, t) = (theta * s - zeta * s_out, theta * s + zeta * s_out);
        if r < 0.0 {
            // 遠い上方散乱: 同様に、同じ指数因子をもつ項どうしをまとめる
            (-p * p).exp() * (erfcx(-r) - erfcx(p)) + (-q * q).exp() * (erfcx(q) - erfcx(t))
        } else {
            scaled(q, q) - scaled(p, p) + ((e - e_out) / kt).exp() * (erf(r) + erf(t))
        }
    };
    // 打ち消し合いの丸めで負になった裾は 0 にする
    (theta * theta / (2.0 * e) * bracket).max(0.0)
}

/// 入射エネルギー e での ∫_g' S(E→E') dE'（群ごと、和は N(E)）。
/// 群構造の外に散乱する分は端の群に入れる（ほかの反応と同じく、中性子を失わない）。
pub fn free_gas_group_integrals(groups: &GroupStructure, e: f64, awr: f64, kt: f64) -> Vec<f64> {
    free_gas_group_moments(groups, e, awr, kt, 0)
        .into_iter()
        .next()
        .unwrap_or_default()
}

/// 入射エネルギー e での ∫_g' S_l(E→E') dE'（[l][群]、l = 0..=order）。S_l は自由ガスの二重微分断面積の
/// ルジャンドルのモーメント（[`free_gas_kernel_moments`]）。
pub fn free_gas_group_moments(
    groups: &GroupStructure,
    e: f64,
    awr: f64,
    kt: f64,
    order: usize,
) -> Vec<Vec<f64>> {
    // 核の幅の目安。重い核では E の近くに集中し、軽い核では kT 程度に広がる
    let width = 2.0 * (e * kt / awr).sqrt() + kt / awr;
    let reach = e + 30.0 * width + 40.0 * kt;
    // 下方散乱は静止標的なら αE まで届き、熱運動でさらに広がる
    let alpha = ((awr - 1.0) / (awr + 1.0)).powi(2);
    let lower = (alpha * e - 30.0 * width).max(0.0);
    let kernel = |x: f64| free_gas_kernel_moments(e, x, awr, kt, order);

    // 区間の境界: 群境界、E' = E（核の折れ点）、αE、E の周りの核の幅ごとの点
    let mut points: Vec<f64> = groups
        .bounds()
        .iter()
        .copied()
        .chain((-60..=60).map(|k| e + 0.5 * f64::from(k) * width))
        .chain([lower, reach, e, alpha * e])
        .filter(|&x| x >= lower && x <= reach)
        .collect();
    points.sort_by(f64::total_cmp);
    points.dedup();
    integrate_points(groups, &points, order + 1, &kernel)
}

/// 自由ガスの二重微分断面積のルジャンドルのモーメント S_l(E→E') = ∫ S(E→E', μ) P_l(μ) dμ（自由原子の断面積を
/// 1 としたとき、1/eV）。P0 は解析式（[`free_gas_kernel`]）、P1 以上は
///   S_l = σb A/(4E) ∫ e^(-β/2) S(α, β) P_l(μ(α)) dα、S(α, β) = e^(-(α²+β²)/(4α))/√(4πα)、σb = ((A+1)/A)²、
///   μ(α) = (E + E' − α A kT)/(2√(EE'))
/// を [`gaussian_law_moments`]（T_eff = T）で積分する。
pub fn free_gas_kernel_moments(e: f64, e_out: f64, awr: f64, kt: f64, order: usize) -> Vec<f64> {
    let mut moments = vec![0.0; order + 1];
    moments[0] = free_gas_kernel(e, e_out, awr, kt);
    if order == 0 || e <= 0.0 || e_out <= 0.0 {
        return moments;
    }
    let scale = awr * kt;
    let (root, root_out) = (e.sqrt(), e_out.sqrt());
    let mu = |alpha: f64| (e + e_out - alpha * scale) / (2.0 * root * root_out);
    let mut integrals = vec![0.0; order + 1];
    gaussian_law_moments(
        (root - root_out).powi(2) / scale,
        (root + root_out).powi(2) / scale,
        (e_out - e) / kt,
        1.0,
        &mu,
        &mut integrals,
    );
    let sigma_bound = ((awr + 1.0) / awr).powi(2);
    let factor = sigma_bound * awr / (4.0 * e);
    for (m, i) in moments.iter_mut().zip(integrals).skip(1) {
        *m = factor * i;
    }
    moments
}

/// 短時間衝突近似の形の散乱則
///   S(α, β) = e^(−(α − |β|)² r/(4α) − |β|/2)/√(4πα/r)、r = T/T_eff
/// （r = 1 で自由ガスの S）について、∫ e^(−β/2) S(α, β) P_l(μ(α)) dα（α は [alpha_lo, alpha_hi]、l = 0..=moments.len()−1）を
/// `moments` に足す。α, β は物質の温度 T の値（β は符号つき）、`mu` は α での実験室系の散乱角の余弦。
///
/// e^(−β/2) と S を別に計算すると、遠い下方散乱（−β が約 1400 を超える）で e^(−β/2) が桁あふれするので、指数を
/// まとめた e^(−β/2) S = e^(−(α − |β|)² r/(4α) − (|β| + β)/2)/√(4πα/r) を使う。u = √α で積分する（被積分関数は
/// u − |β|/u の幅 √(2/r) ほどのガウス形になり、S dα = 2√r e^(...)/√(4π) du）。
pub(super) fn gaussian_law_moments(
    alpha_lo: f64,
    alpha_hi: f64,
    beta: f64,
    ratio: f64,
    mu: &impl Fn(f64) -> f64,
    moments: &mut [f64],
) {
    let b = beta.abs();
    // 指数 (u − |β|/u)² r/4 が 50 を超える u は寄与しない
    let c = (200.0 / ratio).sqrt();
    let root = (c * c + 4.0 * b).sqrt();
    let u_lo = alpha_lo.max(0.0).sqrt().max(0.5 * (root - c));
    let u_hi = alpha_hi.sqrt().min(0.5 * (root + c));
    if u_hi <= u_lo {
        return;
    }
    let step = 0.1 / ratio.sqrt().max(1.0);
    let pieces = ((u_hi - u_lo) / step).ceil().max(1.0) as usize;
    let width = (u_hi - u_lo) / pieces as f64;
    let order = moments.len() - 1;
    let norm = 2.0 * (ratio / (4.0 * std::f64::consts::PI)).sqrt();
    let shift = 0.5 * (b + beta);
    let half = 0.5 * width;
    for k in 0..pieces {
        let a = u_lo + k as f64 * width;
        for &(x, w) in &GAUSS_LEGENDRE_8 {
            let u = a + half * (x + 1.0);
            let alpha = u * u;
            let v = u - b / u;
            let s = w * half * norm * (-0.25 * v * v * ratio - shift).exp();
            for (m, p) in moments.iter_mut().zip(legendre(mu(alpha), order)) {
                *m += s * p;
            }
        }
    }
}

/// 区間の境界 `points` の間ごとに、ベクトル値の関数 f（長さ n）を積分して群に足す [n][群]。
/// 群構造の外の区間は端の群に入れる。
pub(super) fn integrate_points(
    groups: &GroupStructure,
    points: &[f64],
    n: usize,
    f: &impl Fn(f64) -> Vec<f64>,
) -> Vec<Vec<f64>> {
    let highest = groups.bounds()[0];
    let mut integrals = vec![vec![0.0; groups.len()]; n];
    for w in points.windows(2) {
        let (a, b) = (w[0], w[1]);
        let mid = 0.5 * (a + b);
        let g = groups
            .find(mid)
            .unwrap_or(if mid > highest { 0 } else { groups.len() - 1 });
        let values = if a > 0.0 {
            integrate_log_multi(f, a, b, n)
        } else {
            integrate_linear_multi(f, a, b, n)
        };
        for (row, v) in integrals.iter_mut().zip(values) {
            row[g] += v;
        }
    }
    integrals
}

/// [a, b]（a > 0）をエネルギー比が一定以下になるよう対数的に分割して、ベクトル値の関数を積分する。
fn integrate_log_multi(f: &impl Fn(f64) -> Vec<f64>, a: f64, b: f64, n: usize) -> Vec<f64> {
    let pieces = ((b / a).ln() / MAX_RATIO.ln()).ceil().max(1.0) as usize;
    let ratio = (b / a).powf(1.0 / pieces as f64);
    let mut lo = a;
    let mut sum = vec![0.0; n];
    for i in 0..pieces {
        let hi = if i + 1 == pieces { b } else { lo * ratio };
        gauss_legendre_multi(f, lo, hi, &mut sum);
        lo = hi;
    }
    sum
}

/// 0 を含む区間 [a, b] を等分して、ベクトル値の関数を積分する。
fn integrate_linear_multi(f: &impl Fn(f64) -> Vec<f64>, a: f64, b: f64, n: usize) -> Vec<f64> {
    const PIECES: usize = 16;
    let width = (b - a) / PIECES as f64;
    let mut sum = vec![0.0; n];
    for i in 0..PIECES {
        gauss_legendre_multi(
            f,
            a + i as f64 * width,
            a + (i + 1) as f64 * width,
            &mut sum,
        );
    }
    sum
}

fn gauss_legendre_multi(f: &impl Fn(f64) -> Vec<f64>, lo: f64, hi: f64, sum: &mut [f64]) {
    let (mid, half) = (0.5 * (lo + hi), 0.5 * (hi - lo));
    for &(x, wk) in &GAUSS_LEGENDRE_4 {
        for (s, v) in sum.iter_mut().zip(f(mid + half * x)) {
            *s += wk * half * v;
        }
    }
}

/// 1つの積分区間の上限と下限のエネルギー比の上限
const MAX_RATIO: f64 = 1.05;

#[cfg(test)]
mod tests;
