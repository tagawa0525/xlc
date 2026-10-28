//! 非分離共鳴の幅のゆらぎ（平均 1 の χ² 分布）による期待値。
//!
//! 自由度 ν の χ² 分布の密度は k^k x^(k-1) e^(-kx) / Γ(k)（k = ν/2）。x = u² と置くと
//! 2 k^k u^(2k-1) e^(-k u²) / Γ(k) du になり、ν = 1 でも原点で有限になるので、u について適応 Gauss-Kronrod 積分する。

use super::quadrature::gauss_kronrod;
use crate::endf::WidthDof;
use crate::numeric::NumericError;
use crate::special::{erfcx, exp_e1};

/// u の積分の上限を k u² = TAIL とする。密度は e^(-40) ≈ 4e-18 倍になり、倍精度で無視できる
const TAIL: f64 = 40.0;

/// これより大きい c では、閉じた式の代わりに 1/c の級数を使う。閉じた式は
/// E[x²/(x+c)] = 1 - c E[x/(x+c)] の引き算で相対誤差が eps c² 程度になり（c = 100 で 2e-12）、
/// 級数の項 m_n/c^n は c ≥ 100 なら 20 項ほどで倍精度の大きさを下回るため
const SERIES_ABOVE: f64 = 100.0;

/// 自由度 `dof` の χ² 分布での [E[x/(x + c)], E[x²/(x + c)]]（c > 0）。
///
/// x/(x+c) = 1 - c/(x+c) と x²/(x+c) = x - c x/(x+c) により E[1/(x+c)] に帰着させる。
/// E[1/(x+c)] = k^k c^(k-1) e^(kc) Γ(1-k, kc)（k = ν/2）は、ν = 1, 3 では erfc、ν = 2, 4 では E1 で書ける。
pub(super) fn width_ratios(dof: WidthDof, c: f64) -> Result<[f64; 2], NumericError> {
    let k = f64::from(dof.number()) / 2.0;
    let sqrt_pi = std::f64::consts::PI.sqrt();
    let inverse = match dof {
        WidthDof::Fixed => {
            let r = 1.0 / (1.0 + c);
            return Ok([r, r]);
        }
        _ if c > SERIES_ABOVE => return width_ratio_series(dof, c),
        WidthDof::One => (std::f64::consts::PI * k / c).sqrt() * erfcx((k * c).sqrt()),
        WidthDof::Two => exp_e1(c)?,
        WidthDof::Three => 2.0 * k - 2.0 * sqrt_pi * k.powf(1.5) * c.sqrt() * erfcx((k * c).sqrt()),
        WidthDof::Four => k - k * k * c * exp_e1(k * c)?,
    };
    let r = 1.0 - c * inverse;
    Ok([r, 1.0 - c * r])
}

/// 大きな c での [E[x/(x+c)], E[x²/(x+c)]] の漸近級数。t_n = m_n/c^n（m_n = Π_{i<n} (1 + 2i/ν) は
/// 平均 1 の χ²_ν の積率）として、E[x/(x+c)] = Σ_{n≥1} (-1)^(n+1) t_n、E[x²/(x+c)] = Σ_{n≥1} (-1)^(n+1) c t_(n+1)。
fn width_ratio_series(dof: WidthDof, c: f64) -> Result<[f64; 2], NumericError> {
    const MAX_TERMS: u32 = 60;
    let nu = f64::from(dof.number());
    let (mut first, mut second) = (0.0, 0.0);
    let mut t = 1.0 / c;
    let mut sign = 1.0;
    for n in 1..=MAX_TERMS {
        let next = t * (1.0 + 2.0 * f64::from(n) / nu) / c;
        first += sign * t;
        second += sign * c * next;
        if t <= f64::EPSILON * first.abs() && c * next <= f64::EPSILON * second.abs() {
            return Ok([first, second]);
        }
        t = next;
        sign = -sign;
    }
    Err(NumericError::NotConverged(format!(
        "width ratio series for c = {c}"
    )))
}

/// 自由度 `dof` の χ² 分布での f(x) の期待値（相対許容誤差 `tolerance`）。
/// 自由度 0 は幅がゆらがない（x = 1）ことを表す。
///
/// 期待値を入れ子にするときは、外側の許容誤差を内側より十分大きくする。内側の結果は、分割の仕方が
/// 切り替わるところで許容誤差の大きさだけ飛ぶので、外側が同じ精度を求めると分割が止まらない。
/// f のエラーはそのまま返す。
pub(super) fn expectation<E: From<NumericError>>(
    dof: WidthDof,
    tolerance: f64,
    f: impl Fn(f64) -> Result<f64, E>,
) -> Result<f64, E> {
    let values = expectation_vec::<E>(dof, tolerance, 1, |x, out| {
        out[0] = f(x)?;
        Ok(())
    })?;
    Ok(values[0])
}

/// `expectation` のベクトル版。f は x と長さ `n` の出力を受け取り、値を書き込む。
pub(super) fn expectation_vec<E: From<NumericError>>(
    dof: WidthDof,
    tolerance: f64,
    n: usize,
    f: impl Fn(f64, &mut [f64]) -> Result<(), E>,
) -> Result<Vec<f64>, E> {
    let k = f64::from(dof.number()) / 2.0;
    // Γ(k)（k = 1/2, 1, 3/2, 2）
    let gamma = match dof {
        WidthDof::Fixed => {
            let mut out = vec![0.0; n];
            f(1.0, &mut out)?;
            return Ok(out);
        }
        WidthDof::One => std::f64::consts::PI.sqrt(),
        WidthDof::Two | WidthDof::Four => 1.0,
        WidthDof::Three => 0.5 * std::f64::consts::PI.sqrt(),
    };
    let scale = 2.0 * k.powf(k) / gamma;
    let integrand = |u: f64, out: &mut [f64]| {
        let x = u * u;
        f(x, out)?;
        let density = scale * u.powf(2.0 * k - 1.0) * (-k * x).exp();
        out.iter_mut().for_each(|v| *v *= density);
        Ok(())
    };
    gauss_kronrod(integrand, n, 0.0, (TAIL / k).sqrt(), tolerance)
}

#[cfg(test)]
mod tests;
