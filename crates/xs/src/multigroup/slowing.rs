//! 超微細群（UFG）の減速計算による、σ0 ごとの群定数の重み。
//!
//! NJOY GROUPR の flux calculator と同じく、吸収体と背景の減速材の無限の均質な混合物を考える。背景の散乱は
//! 重み関数 w(E) どおりの滑らかなスペクトルを与え、吸収体の弾性散乱（重心系で等方、標的は静止）による減速は
//! 積分形の減速方程式で解く。レサジー u = ln(E_high/E) と F(u) = E φ(E)、W(u) = E w(E) で
//!
//! F(u) (σt + σ0) = σ0 W(u) + 1/(1 − α) ∫_{u−Δ}^{u} σs(u') F(u') e^{−(u−u')} du'、α = ((A − 1)/(A + 1))²、Δ = ln(1/α)
//!
//! を高いエネルギーから順に解く。σs F は格子点の間で線形とし、線形 × 指数関数を区間ごとに解析的に積分する。
//! 上端 E_high より上（最初の窓）は NR 近似とし、吸収体の散乱源はポテンシャル散乱 σp が重み関数のスペクトルで
//! 散乱する σp W(u) とする（NJOY GROUPR の flux calculator と同じ。共鳴の散乱 σs(u) F(u) を使うと、範囲のすぐ上の
//! 散乱の共鳴が、範囲の上端の中性子束を NJOY より大きく見積もる）。上端の値は F = W (σp + σ0) / (σt + σ0)。

use super::{SlowingDown, Weight};
use crate::endf::Tab1;

/// 1つの窓（散乱で届くレサジーの幅 Δ）に入れる格子点の最少の数。格子はポイントワイズの断面積の格子に、これより
/// 粗いところだけ点を足す（`finer_grids_do_not_change_the_flux` で収束を確かめた値）。lattice のピンごとの
/// 減速計算も同じ値を使っている（lattice の `ufg` を参照）
pub(crate) const WINDOW_POINTS: usize = 16;

/// 減速計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum SlowingDownError {
    #[error("energy range {0}..{1} eV must be positive and increasing")]
    Range(f64, f64),
    #[error(
        "atomic weight ratio {0} must be greater than 1 (hydrogen needs no slowing-down weight)"
    )]
    Mass(f64),
    #[error("background cross section {0} b must be positive")]
    Sigma0(f64),
    #[error("potential scattering {0} b must be positive")]
    Potential(f64),
    #[error("cross sections are not defined at {0} eV")]
    Undefined(f64),
    #[error("a scattering window needs at least 2 grid points, got {0}")]
    WindowPoints(usize),
}

/// 範囲 `slowing.range` = [e_low, e_high] (eV) での中性子束の形 f(E) = φ(E) / w(E)（エネルギーの昇順の lin-lin の表）。
pub fn slowing_down_shape(
    total: &Tab1,
    elastic: &Tab1,
    slowing: &SlowingDown,
    weight: &Weight,
    sigma0: f64,
) -> Result<Tab1, SlowingDownError> {
    slowing_down_shape_with(total, elastic, slowing, weight, sigma0, WINDOW_POINTS)
}

/// 窓あたりの最少の格子点の数 `window_points` を指定する `slowing_down_shape`。
pub fn slowing_down_shape_with(
    total: &Tab1,
    elastic: &Tab1,
    slowing: &SlowingDown,
    weight: &Weight,
    sigma0: f64,
    window_points: usize,
) -> Result<Tab1, SlowingDownError> {
    let SlowingDown {
        awr,
        potential,
        range: (e_low, e_high),
    } = *slowing;
    if !(e_low > 0.0 && e_high > e_low && e_high.is_finite()) {
        return Err(SlowingDownError::Range(e_low, e_high));
    }
    if !(awr > 1.0 && awr.is_finite()) {
        return Err(SlowingDownError::Mass(awr));
    }
    if !(sigma0 > 0.0 && sigma0.is_finite()) {
        return Err(SlowingDownError::Sigma0(sigma0));
    }
    if !(potential > 0.0 && potential.is_finite()) {
        return Err(SlowingDownError::Potential(potential));
    }
    if window_points < 2 {
        return Err(SlowingDownError::WindowPoints(window_points));
    }
    let alpha = ((awr - 1.0) / (awr + 1.0)).powi(2);
    let width = -alpha.ln();
    let one_minus_alpha = 1.0 - alpha;
    let e_top = e_high / alpha;

    // 格子（エネルギーの降順）: ポイントワイズの断面積の格子に、窓あたり window_points 点より粗いところだけ点を足す
    let mut energies: Vec<f64> = total
        .x
        .iter()
        .chain(&elastic.x)
        .copied()
        .filter(|&e| e > e_low && e < e_top)
        .chain([e_low, e_high, e_top])
        .collect();
    energies.sort_by(|a, b| b.total_cmp(a));
    energies.dedup();
    let step = width / window_points as f64;
    // レサジーと、その点のエネルギー。格子の点のエネルギーは計算し直さずに持つ
    // （計算し直すと丸めで範囲の端からずれ、範囲の端で表を引けなくなる）
    let mut u: Vec<f64> = Vec::with_capacity(energies.len());
    let mut points: Vec<f64> = Vec::with_capacity(energies.len());
    for e in energies {
        let ue = (e_high / e).ln();
        // 別のエネルギーでもレサジーが丸めで同じなら、同じ点とする（幅 0 の区間を作らない）。
        // 範囲の端なら、端のエネルギーのほうを残す（形の表がちょうど端で始まり終わるように）
        if u.last().is_some_and(|&last| ue <= last) {
            if (e == e_low || e == e_high)
                && let Some(p) = points.last_mut()
            {
                *p = e;
            }
            continue;
        }
        if let Some(&last) = u.last() {
            let gap: f64 = ue - last;
            let extra = (gap / step).ceil() as usize;
            for k in 1..extra {
                let fill = last + gap * k as f64 / extra as f64;
                u.push(fill);
                points.push(e_high * (-fill).exp());
            }
        }
        u.push(ue);
        points.push(e);
    }

    let at = |t: &Tab1, e: f64| t.eval(e).ok_or(SlowingDownError::Undefined(e));
    let mut f = Vec::with_capacity(u.len());
    // g = σs F（格子点の値。最初の窓は σp W）
    let mut g = Vec::with_capacity(u.len());
    for (i, (&ui, &e)) in u.iter().zip(&points).enumerate() {
        let big_w = e * weight.value(e);
        if ui < 0.0 {
            // 範囲より上: 散乱源だけを使い、中性子束は表に入れない（下でエネルギーが上端以下の点だけを残す）
            f.push(f64::NAN);
            g.push(potential * big_w);
            continue;
        }
        let (st, ss) = (at(total, e)?, at(elastic, e)?);
        let (value, source) = if ui <= 0.0 {
            // 最初の窓: NR 近似。散乱源はポテンシャル散乱が重み関数のスペクトルで散乱するもの
            (
                big_w * (potential + sigma0) / (st + sigma0),
                potential * big_w,
            )
        } else {
            // 窓 [ui − Δ, ui] の積分のうち、既知の部分と、最後の区間の g_{i−1} の項
            let lower = ui - width;
            let k = u.partition_point(|&x| x < lower);
            let mut known = 0.0;
            if k > 0 && u[k] > lower {
                let t = (lower - u[k - 1]) / (u[k] - u[k - 1]);
                let g_lower = g[k - 1] + t * (g[k] - g[k - 1]);
                let (ca, cb) = segment(lower, u[k], ui);
                known += ca * g_lower + cb * g[k];
            }
            for j in k..i - 1 {
                let (ca, cb) = segment(u[j], u[j + 1], ui);
                known += ca * g[j] + cb * g[j + 1];
            }
            let (ca, cb) = segment(u[i - 1], ui, ui);
            known += ca * g[i - 1];
            let value = (sigma0 * big_w + known / one_minus_alpha)
                / (st + sigma0 - ss * cb / one_minus_alpha);
            (value, ss * value)
        };
        f.push(value);
        g.push(source);
    }

    // [e_low, e_high] の点を、形 f = F / W にしてエネルギーの昇順に並べる
    let (x, y): (Vec<f64>, Vec<f64>) = points
        .iter()
        .zip(&f)
        .rev()
        .filter(|&(&e, _)| e <= e_high)
        .map(|(&e, &value)| (e, value / (e * weight.value(e))))
        .unzip();
    Ok(Tab1::linear(x, y))
}

/// 区間 [a, b] で線形な g の ∫_a^b g(u') e^{−(U − u')} du' を、両端の値 g(a), g(b) の係数 (ca, cb) で表す。
fn segment(a: f64, b: f64, big_u: f64) -> (f64, f64) {
    let s = b - a;
    // 1 − e^{−s} と、(s − (1 − e^{−s})) / s
    let em = -(-s).exp_m1();
    let tail = (s - em) / s;
    let decay = (-(big_u - b)).exp();
    (decay * (em - tail), decay * tail)
}

#[cfg(test)]
mod tests;
