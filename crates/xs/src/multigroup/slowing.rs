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
mod tests {
    use super::*;
    use crate::endf::Tab1;
    use crate::multigroup::SlowingDown;

    const RANGE: (f64, f64) = (1.0, 1.0e3);

    fn constant(value: f64) -> Tab1 {
        Tab1::linear(vec![1.0e-5, 2.0e7], vec![value, value])
    }

    /// 範囲 RANGE の減速計算（ポテンシャル散乱 `potential` b）
    fn slowing(awr: f64, potential: f64) -> SlowingDown {
        SlowingDown {
            awr,
            potential,
            range: RANGE,
        }
    }

    #[test]
    fn pure_scatterer_keeps_the_weight_spectrum() {
        let (total, elastic) = (constant(10.0), constant(10.0));
        let shape = slowing_down_shape(
            &total,
            &elastic,
            &slowing(236.0, 10.0),
            &Weight::lwr(),
            50.0,
        )
        .unwrap();
        for y in &shape.y {
            assert!((y - 1.0).abs() < 1e-12, "{y}");
        }
    }

    #[test]
    fn source_above_the_range_is_potential_scattering() {
        // 範囲の上の吸収体の散乱は、NR 近似どおりポテンシャル散乱 σp が重み関数のスペクトルで散乱するとする
        // （NJOY GROUPR の flux calculator と同じ）。範囲の上に散乱の共鳴があっても、範囲の中の純粋な
        // 散乱体（σs = σp）の中性子束は重み関数のまま。共鳴は上端のすぐ上（最初の窓の中）に置く
        let above = 1.005 * RANGE.1;
        let x = vec![1.0e-5, above, above, 2.0e7];
        let step = Tab1::linear(x, vec![10.0, 10.0, 30.0, 30.0]);
        let shape =
            slowing_down_shape(&step, &step, &slowing(236.0, 10.0), &Weight::lwr(), 50.0).unwrap();
        for y in &shape.y {
            assert!((y - 1.0).abs() < 1e-12, "{y}");
        }
    }

    #[test]
    fn shape_ends_exactly_at_the_range() {
        // 端のエネルギーを計算し直すと丸めでずれ、範囲の端で表を引けなくなる
        let (total, elastic) = (constant(10.0), constant(8.0));
        let range = (0.1, 1.83e5);
        let slowing = SlowingDown {
            awr: 56.44629,
            potential: 8.0,
            range,
        };
        let shape = slowing_down_shape(&total, &elastic, &slowing, &Weight::lwr(), 50.0).unwrap();
        assert_eq!(shape.x.first(), Some(&range.0));
        assert_eq!(shape.x.last(), Some(&range.1));
    }

    #[test]
    fn energies_with_equal_lethargy_give_finite_flux() {
        // 別の値なのに、レサジー ln(E_high / E) が丸めで同じになる2つのエネルギーを探す
        let lethargy = |e: f64| (RANGE.1 / e).ln();
        let mut e = 1.0_f64;
        while lethargy(e.next_up()) != lethargy(e) {
            e = e.next_up();
        }
        let x = vec![1.0e-5, e, e.next_up(), 2.0e7];
        let total = Tab1::linear(x.clone(), vec![10.0, 10.0, 12.0, 12.0]);
        let elastic = Tab1::linear(x, vec![8.0, 8.0, 9.0, 9.0]);
        let shape =
            slowing_down_shape(&total, &elastic, &slowing(236.0, 9.0), &Weight::lwr(), 50.0)
                .unwrap();
        assert!(shape.y.iter().all(|y| y.is_finite()), "{:?}", shape.y);
    }

    #[test]
    fn rejects_fewer_than_two_points_per_window() {
        let (total, elastic) = (constant(10.0), constant(8.0));
        for points in [0, 1] {
            let result = slowing_down_shape_with(
                &total,
                &elastic,
                &slowing(236.0, 8.0),
                &Weight::lwr(),
                50.0,
                points,
            );
            assert!(
                matches!(result, Err(SlowingDownError::WindowPoints(_))),
                "{points}"
            );
        }
    }

    #[test]
    fn constant_absorber_approaches_the_asymptotic_flux() {
        // σa = 2, σs = 8, σ0 = 50: 上端から十分下では F = σ0 / (σa + σ0)
        let (total, elastic) = (constant(10.0), constant(8.0));
        let shape =
            slowing_down_shape(&total, &elastic, &slowing(236.0, 8.0), &Weight::lwr(), 50.0)
                .unwrap();
        let low = shape.eval(2.0).unwrap();
        assert!((low - 50.0 / 52.0).abs() < 1e-10, "{low}");
    }

    /// 10 eV に共鳴（幅 0.5 eV、ピークで捕獲 2000 b・散乱 1000 b）、ほかは散乱 10 b
    fn resonance() -> (Tab1, Tab1) {
        let x: Vec<f64> = (0..=20_000)
            .map(|i| 5.0 + 10.0 * f64::from(i) / 20_000.0)
            .collect();
        let peak = |e: f64| 1.0 / (1.0 + ((e - 10.0) / 0.25).powi(2));
        let elastic: Vec<f64> = x.iter().map(|&e| 10.0 + 1000.0 * peak(e)).collect();
        let total = x
            .iter()
            .zip(&elastic)
            .map(|(&e, s)| s + 2000.0 * peak(e))
            .collect();
        let wrap = |y: Vec<f64>| {
            // 範囲の外は一定値で延ばす
            let mut xs = vec![1.0e-5];
            xs.extend(&x);
            xs.push(2.0e7);
            let mut ys = vec![y[0]];
            ys.extend(&y);
            ys.push(y[y.len() - 1]);
            Tab1::linear(xs, ys)
        };
        (wrap(total), wrap(elastic))
    }

    #[test]
    fn heavy_absorber_scatters_without_losing_energy() {
        // 質量が大きく、散乱でのエネルギーの損失（窓の幅 Δ ≈ 4/A = 4e-4）が共鳴の幅（5%）よりずっと
        // 小さいと、共鳴の中でも F ≈ σ0 / (σa + σ0)（散乱断面積によらない。NR 近似の 1/(σt + σ0) とは違う）
        let (total, elastic) = resonance();
        let shape = slowing_down_shape(
            &total,
            &elastic,
            &slowing(1.0e4, 10.0),
            &Weight::lwr(),
            50.0,
        )
        .unwrap();
        let at_peak = shape.eval(10.0).unwrap();
        assert!(
            (at_peak - 50.0 / 2050.0).abs() < 1e-3 * at_peak,
            "{at_peak}"
        );
    }

    #[test]
    fn finer_grids_do_not_change_the_flux() {
        let (total, elastic) = resonance();
        let coarse = slowing_down_shape(
            &total,
            &elastic,
            &slowing(236.0, 10.0),
            &Weight::lwr(),
            50.0,
        )
        .unwrap();
        let fine = slowing_down_shape_with(
            &total,
            &elastic,
            &slowing(236.0, 10.0),
            &Weight::lwr(),
            50.0,
            4 * WINDOW_POINTS,
        )
        .unwrap();
        for e in [9.0, 9.9, 10.0, 10.1, 11.0, 3.0] {
            let (a, b) = (coarse.eval(e).unwrap(), fine.eval(e).unwrap());
            assert!((a - b).abs() < 1e-4 * b, "{e}: {a} vs {b}");
        }
    }
}
