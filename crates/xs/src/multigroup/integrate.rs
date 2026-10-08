//! 群ごとの積分。

use super::{GroupStructure, Weight};
use xlc_numeric::quadrature::GAUSS_LEGENDRE_4;

/// 各群について ∫_g f(E) w(E) dE を求める。
///
/// f は `breakpoints` の間で滑らかであること（lin-lin の表なら、その x を渡す）。
/// 区間ごとに、エネルギー比が一定以下になるよう対数的に分割して Gauss-Legendre 積分する。
pub fn group_integrals(
    groups: &GroupStructure,
    weight: &Weight,
    breakpoints: &[f64],
    f: impl Fn(f64) -> f64,
) -> Vec<f64> {
    group_integrals_multi(groups, weight, breakpoints, 1, |e, out| out[0] = f(e))
        .into_iter()
        .map(|row| row[0])
        .collect()
}

/// `group_integrals` の多成分版。f は E と、長さ `n` の出力を受け取り、値を書き込む。
/// 戻り値は [入射の群][成分]。
pub fn group_integrals_multi(
    groups: &GroupStructure,
    weight: &Weight,
    breakpoints: &[f64],
    n: usize,
    f: impl Fn(f64, &mut [f64]),
) -> Vec<Vec<f64>> {
    let bounds = groups.bounds();
    let (lowest, highest) = (bounds[groups.len()], bounds[0]);
    let mut points: Vec<f64> = breakpoints
        .iter()
        .chain(bounds)
        .chain(&weight.breakpoints())
        .copied()
        .filter(|e| (lowest..=highest).contains(e))
        .collect();
    points.sort_by(f64::total_cmp);
    points.dedup();

    let mut integrals = vec![vec![0.0; n]; groups.len()];
    let mut values = vec![0.0; n];
    for w in points.windows(2) {
        let (a, b) = (w[0], w[1]);
        // 点は群構造の範囲にあり、群の境界をすべて含むので、区間はちょうど1つの群に入る
        // （境界は降順。中点以上の境界の数 - 1 が群の番号。`GroupStructure::find` と同じ判定）
        let middle = 0.5 * (a + b);
        let g = bounds
            .partition_point(|&bound| bound >= middle)
            .saturating_sub(1)
            .min(groups.len() - 1);
        let pieces = ((b / a).ln() / MAX_RATIO.ln()).ceil().max(1.0) as usize;
        let ratio = (b / a).powf(1.0 / pieces as f64);
        let mut lo = a;
        for i in 0..pieces {
            let hi = if i + 1 == pieces { b } else { lo * ratio };
            let (mid, half) = (0.5 * (lo + hi), 0.5 * (hi - lo));
            for &(x, wk) in &GAUSS_LEGENDRE_4 {
                let e = mid + half * x;
                values.fill(0.0);
                f(e, &mut values);
                let factor = wk * weight.value(e) * half;
                for (sum, v) in integrals[g].iter_mut().zip(&values) {
                    *sum += factor * v;
                }
            }
            lo = hi;
        }
    }
    integrals
}

/// 1つの積分区間の上限と下限のエネルギー比の上限
const MAX_RATIO: f64 = 1.05;

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 1e-10 * expected.abs(),
            "actual {actual}, expected {expected}"
        );
    }

    #[test]
    fn weight_only_in_one_over_e_region() {
        let groups = GroupStructure::new(vec![1.0e5, 1.0e3, 1.0]).unwrap();
        let r = group_integrals(&groups, &Weight::lwr(), &[], |_| 1.0);
        assert_close(r[0], (1.0e5_f64 / 1.0e3).ln());
        assert_close(r[1], (1.0e3_f64).ln());
    }

    #[test]
    fn weight_only_in_maxwellian_region() {
        let w = Weight::lwr();
        let (a, b, tb) = (0.01, 0.05, 0.0253);
        let groups = GroupStructure::new(vec![b, a]).unwrap();
        let r = group_integrals(&groups, &w, &[], |_| 1.0);
        // ∫ ab E exp(-E/tb) dE = ab tb² [(1 + a/tb) e^{-a/tb} - (1 + b/tb) e^{-b/tb}]
        let ab = w.value(a) / (a * (-a / tb).exp());
        let primitive = |e: f64| (1.0 + e / tb) * (-e / tb).exp();
        assert_close(r[0], ab * tb * tb * (primitive(a) - primitive(b)));
    }

    #[test]
    fn linear_function_times_one_over_e() {
        let groups = GroupStructure::new(vec![100.0, 10.0]).unwrap();
        let r = group_integrals(&groups, &Weight::lwr(), &[], |e| e);
        assert_close(r[0], 90.0);
    }

    #[test]
    fn step_at_breakpoint() {
        let groups = GroupStructure::new(vec![100.0, 10.0]).unwrap();
        let r = group_integrals(&groups, &Weight::lwr(), &[30.0], |e| {
            if e < 30.0 { 1.0 } else { 2.0 }
        });
        assert_close(r[0], (3.0_f64).ln() + 2.0 * (100.0_f64 / 30.0).ln());
    }

    #[test]
    fn multi_component_matches_scalar() {
        let groups = GroupStructure::new(vec![100.0, 10.0, 1.0]).unwrap();
        let w = Weight::lwr();
        let r = group_integrals_multi(&groups, &w, &[30.0], 2, |e, out| {
            out[0] = e;
            out[1] = 1.0;
        });
        let e = group_integrals(&groups, &w, &[30.0], |e| e);
        let one = group_integrals(&groups, &w, &[30.0], |_| 1.0);
        for g in 0..2 {
            assert_close(r[g][0], e[g]);
            assert_close(r[g][1], one[g]);
        }
    }

    #[test]
    fn group_spanning_weight_breakpoint() {
        let w = Weight::lwr();
        let groups = GroupStructure::new(vec![1.0, 0.05]).unwrap();
        let r = group_integrals(&groups, &w, &[], |_| 1.0);
        let tb = 0.0253;
        let ab = w.value(0.05) / (0.05 * (-0.05_f64 / tb).exp());
        let primitive = |e: f64| (1.0 + e / tb) * (-e / tb).exp();
        let expected = ab * tb * tb * (primitive(0.05) - primitive(0.1)) + (1.0_f64 / 0.1).ln();
        assert_close(r[0], expected);
    }
}
