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
mod tests;
