//! ドップラー拡がり（NJOY BROADR と同じ SIGMA1 法の核関数）。
//!
//! 0 K の lin-lin のポイントワイズの断面積を、標的核の熱運動（自由ガスの Maxwell 分布）で平均する。
//! 速度に比例する変数 x = √(αE)（α = A/(kT)）を使うと、拡げた断面積は
//!
//!   σ̄(y) = 1/(y²√π) ∫₀^∞ x² σ(x) [exp(-(x-y)²) - exp(-(x+y)²)] dx
//!
//! になる。表の最初の点より下は 1/v で x = 0 まで、最後の点より上は一定値で延ばす。

use crate::endf::{Mf2, Tab1};

/// 拡げる上限エネルギーの既定の最大値 (eV)。NJOY BROADR と同じ。
pub const MAX_LIMIT: f64 = 6.5e6;

/// 0 K の表を、温度 `temperature` (K) に拡げる。
///
/// `tables` は lin-lin の表で、同じエネルギー格子の上でまとめて拡げる。
/// `e_max` 未満のエネルギーは拡げた値を、格子の中点での線形補間の誤差が `tol` 以下になるまで
/// 細かくした格子に置き、`e_max` 以上は元の表の点をそのまま使う。
/// 最初の点が格子全体の最初の点より上にある表（しきい反応）は、しきい値より下を 0 のままにする。
pub fn broaden(tables: &[Tab1], awr: f64, temperature: f64, e_max: f64, tol: f64) -> Vec<Tab1> {
    let _ = (tables, awr, temperature, e_max, tol);
    todo!()
}

/// 拡げる上限エネルギー（NJOY BROADR の既定値と同じ規則）。
///
/// 分離共鳴があればその上限、非分離共鳴だけならその下限とする（どちらも `MAX_LIMIT` まで）。
/// 共鳴パラメータがなければ `MAX_LIMIT` とし、1 eV より上から始まるしきい反応があれば、
/// そのしきい値の手前までにする。
pub fn broadening_limit<'a>(mf2: Option<&Mf2>, tables: impl Iterator<Item = &'a Tab1>) -> f64 {
    let _ = (mf2, tables.count());
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::endf::{Interpolation, RangeKind, ResonanceRange};
    use crate::physics::BK;
    use crate::pointwise::linearize;
    use crate::special::erf;

    const SQRT_PI: f64 = 1.772_453_850_905_516;

    fn alpha(awr: f64, temperature: f64) -> f64 {
        awr / (BK * temperature)
    }

    /// 一定の断面積を拡げた値: σ [(1 + 1/(2y²)) erf(y) + exp(-y²)/(y√π)]。
    fn broadened_constant(sigma: f64, y: f64) -> f64 {
        sigma * ((1.0 + 0.5 / (y * y)) * erf(y) + (-y * y).exp() / (y * SQRT_PI))
    }

    fn below<'a>(t: &'a Tab1, e_max: f64) -> impl Iterator<Item = (f64, f64)> + 'a {
        t.x.iter()
            .zip(&t.y)
            .map(|(&e, &s)| (e, s))
            .filter(move |&(e, _)| e < e_max)
    }

    #[test]
    fn constant_cross_section_matches_analytic_broadening() {
        let (awr, temperature, e_max) = (1.0, 300.0, 1.0e3);
        let a = alpha(awr, temperature);
        // 最初の点を十分低くして、1/v で延ばす部分の影響をなくす
        let table = Tab1::linear(vec![1.0e-11, 1.0e5], vec![10.0, 10.0]);
        let [out] = broaden(&[table], awr, temperature, e_max, 1.0e-3)
            .try_into()
            .unwrap();
        assert!(out.x.len() > 10);
        for (e, s) in below(&out, e_max) {
            let exact = broadened_constant(10.0, (a * e).sqrt());
            assert!(
                (s - exact).abs() <= 1.0e-9 * exact,
                "at {e}: {s} vs {exact}"
            );
        }
    }

    #[test]
    fn grid_resolves_broadened_cross_section() {
        let (awr, temperature, e_max, tol) = (1.0, 300.0, 1.0e3, 1.0e-3);
        let a = alpha(awr, temperature);
        let table = Tab1::linear(vec![1.0e-11, 1.0e5], vec![10.0, 10.0]);
        let [out] = broaden(&[table], awr, temperature, e_max, tol)
            .try_into()
            .unwrap();
        for w in out.x.windows(2).filter(|w| w[1] <= e_max) {
            let mid = 0.5 * (w[0] + w[1]);
            let exact = broadened_constant(10.0, (a * mid).sqrt());
            let approx = out.eval(mid).unwrap();
            assert!((approx - exact).abs() <= tol * exact, "at {mid}");
        }
    }

    #[test]
    fn one_over_v_is_unchanged() {
        let (awr, temperature, e_max) = (236.0, 600.0, 1.0e3);
        let exact = |e: f64| 1.0 / e.sqrt();
        let source = Tab1 {
            breaks: vec![2],
            laws: vec![Interpolation::LogLog],
            x: vec![1.0e-5, 1.0e5],
            y: vec![exact(1.0e-5), exact(1.0e5)],
        };
        let table = linearize(&source, 1.0e-5);
        let [out] = broaden(&[table], awr, temperature, e_max, 1.0e-3)
            .try_into()
            .unwrap();
        for (e, s) in below(&out, e_max) {
            assert!((s - exact(e)).abs() <= 2.0e-5 * exact(e), "at {e}");
        }
    }

    /// 被積分関数を、元の表の区間ごとに Simpson 則で積分した拡げた断面積。
    fn simpson_reference(table: &Tab1, a: f64, e: f64) -> f64 {
        let y = (a * e).sqrt();
        let integrand = |x: f64| {
            let sigma = table.eval(x * x / a).unwrap_or(0.0);
            x * x * sigma * ((-(x - y) * (x - y)).exp() - (-(x + y) * (x + y)).exp())
        };
        let (lo, hi) = ((y - 6.0).max(0.0), y + 6.0);
        let mut sum = 0.0;
        for w in table.x.windows(2) {
            let (x0, x1) = ((a * w[0]).sqrt().max(lo), (a * w[1]).sqrt().min(hi));
            if x1 <= x0 {
                continue;
            }
            const N: usize = 20_000;
            let h = (x1 - x0) / N as f64;
            let mut s = integrand(x0) + integrand(x1);
            for i in 1..N {
                let weight = if i % 2 == 1 { 4.0 } else { 2.0 };
                s += weight * integrand(x0 + h * i as f64);
            }
            sum += s * h / 3.0;
        }
        sum / (y * y * SQRT_PI)
    }

    #[test]
    fn narrow_resonance_matches_direct_integration() {
        // ドップラー幅（約 0.09 eV）より狭い三角形の共鳴
        let (awr, temperature, e_max) = (236.0, 600.0, 1.0e3);
        let table = Tab1::linear(
            vec![1.0e-5, 9.99, 10.0, 10.01, 1.0e5],
            vec![10.0, 10.0, 1000.0, 10.0, 10.0],
        );
        let [out] = broaden(
            std::slice::from_ref(&table),
            awr,
            temperature,
            e_max,
            1.0e-3,
        )
        .try_into()
        .unwrap();
        let a = alpha(awr, temperature);
        let near: Vec<_> = below(&out, e_max)
            .filter(|&(e, _)| (9.5..10.5).contains(&e))
            .collect();
        assert!(near.len() > 10);
        for (e, s) in near {
            let exact = simpson_reference(&table, a, e);
            assert!(
                (s - exact).abs() <= 1.0e-8 * exact,
                "at {e}: {s} vs {exact}"
            );
        }
    }

    #[test]
    fn points_at_and_above_limit_are_kept() {
        let (awr, temperature, e_max) = (236.0, 600.0, 1.0e3);
        let table = Tab1::linear(
            vec![1.0e-5, 1.0e3, 1.0e3, 2.0e3, 1.0e5],
            vec![10.0, 10.0, 20.0, 30.0, 40.0],
        );
        let [out] = broaden(
            std::slice::from_ref(&table),
            awr,
            temperature,
            e_max,
            1.0e-3,
        )
        .try_into()
        .unwrap();
        let k = out.x.partition_point(|&e| e < e_max);
        // e_max では拡げた値（段差がならされて中間の値）と、元の表の右側の値の2点を持つ
        assert_eq!(out.x[k], e_max);
        assert!((out.y[k] - 15.0).abs() < 0.5, "{}", out.y[k]);
        assert_eq!(&out.x[k + 1..], &[1.0e3, 2.0e3, 1.0e5]);
        assert_eq!(&out.y[k + 1..], &[20.0, 30.0, 40.0]);
    }

    #[test]
    fn threshold_reaction_stays_zero_below_threshold() {
        let (awr, temperature, e_max) = (236.0, 600.0, 1.0e3);
        let elastic = Tab1::linear(vec![1.0e-5, 1.0e5], vec![10.0, 10.0]);
        let inelastic = Tab1::linear(vec![100.0, 1.0e5], vec![0.0, 5.0]);
        let [elastic, inelastic] = broaden(&[elastic, inelastic], awr, temperature, e_max, 1.0e-3)
            .try_into()
            .unwrap();
        assert_eq!((inelastic.x[0], inelastic.y[0]), (100.0, 0.0));
        // しきい値より上は拡げた値で、弾性散乱と同じ格子に乗る
        let above: Vec<f64> = below(&inelastic, e_max).map(|(e, _)| e).collect();
        assert!(above.len() > 2);
        assert!(above.iter().all(|e| elastic.x.contains(e)));
        assert!(below(&inelastic, e_max).skip(1).all(|(_, s)| s > 0.0));
    }

    fn range(el: f64, eh: f64, kind: RangeKind) -> ResonanceRange {
        ResonanceRange { el, eh, kind }
    }

    fn mf2(ranges: Vec<ResonanceRange>) -> Mf2 {
        Mf2 {
            za: 92238.0,
            awr: 236.0,
            ranges,
        }
    }

    fn mlbw() -> RangeKind {
        RangeKind::Mlbw(crate::endf::Mlbw {
            spi: 0.0,
            ap: 0.9,
            naps: 0,
            l_values: vec![],
        })
    }

    #[test]
    fn limit_is_top_of_resolved_range() {
        let tables = [Tab1::linear(vec![1.0e-5, 2.0e7], vec![1.0, 1.0])];
        let resolved = mf2(vec![
            range(1.0e-5, 2.0e4, mlbw()),
            range(2.0e4, 1.5e5, RangeKind::Unresolved { lssf: 1 }),
        ]);
        assert_eq!(broadening_limit(Some(&resolved), tables.iter()), 2.0e4);
        let wide = mf2(vec![range(1.0e-5, 2.0e7, mlbw())]);
        assert_eq!(broadening_limit(Some(&wide), tables.iter()), MAX_LIMIT);
    }

    #[test]
    fn limit_is_bottom_of_unresolved_range_without_resolved_range() {
        let tables = [Tab1::linear(vec![1.0e-5, 2.0e7], vec![1.0, 1.0])];
        let unresolved = mf2(vec![range(3.0e4, 1.0e5, RangeKind::Unresolved { lssf: 1 })]);
        assert_eq!(broadening_limit(Some(&unresolved), tables.iter()), 3.0e4);
    }

    #[test]
    fn limit_without_resonances_stops_before_thresholds() {
        let low = Tab1::linear(vec![1.0e-5, 2.0e7], vec![1.0, 1.0]);
        let near_zero = Tab1::linear(vec![0.5, 2.0e7], vec![0.0, 1.0]);
        let high = Tab1::linear(vec![2.2e6, 2.0e7], vec![0.0, 1.0]);
        let radius = mf2(vec![range(
            1.0e-5,
            2.0e7,
            RangeKind::ScatteringRadius { ap: 0.5 },
        )]);
        let tables = [low.clone(), near_zero.clone()];
        assert_eq!(broadening_limit(Some(&radius), tables.iter()), MAX_LIMIT);
        assert_eq!(broadening_limit(None, tables.iter()), MAX_LIMIT);
        let tables = [low, near_zero, high];
        assert_eq!(broadening_limit(None, tables.iter()), 0.99999 * 2.2e6);
    }
}
