//! 出射エネルギー分布の群への振り分け。

use super::GroupStructure;
use crate::endf::{Interpolation, OutgoingDistributions, Tab1};
use crate::pointwise::linearize;

/// 出射エネルギーの確率密度の表を、群ごとに積分する。
///
/// 群構造の外に出る分は端の群に入れる（NJOY GROUPR と同じ）。たとえば 10 MeV を超える
/// 核分裂中性子は最上位の群に入る。数えずに捨てると、χ を正規化したときにその分が
/// すべての群に比例配分され、熱群にまで配られてしまう。
pub fn outgoing_group_integrals(groups: &GroupStructure, density: &Tab1) -> Vec<f64> {
    let linear = linearize(density, DENSITY_TOLERANCE);
    let bounds = groups.bounds();
    let mut integrals = vec![0.0; groups.len()];
    for k in 0..linear.x.len().saturating_sub(1) {
        let (x0, x1) = (linear.x[k], linear.x[k + 1]);
        if x1 <= x0 {
            continue;
        }
        let (y0, y1) = (linear.y[k], linear.y[k + 1]);
        let y = |e: f64| y0 + (y1 - y0) * (e - x0) / (x1 - x0);
        // 区間を群境界で切り、それぞれを台形公式で積分する（区間内は線形なので厳密）
        let mut cuts: Vec<f64> = bounds
            .iter()
            .copied()
            .filter(|&b| b > x0 && b < x1)
            .collect();
        cuts.push(x0);
        cuts.push(x1);
        cuts.sort_by(f64::total_cmp);
        for w in cuts.windows(2) {
            let (a, b) = (w[0], w[1]);
            let mid = 0.5 * (a + b);
            let g = groups
                .find(mid)
                .unwrap_or(if mid > bounds[0] { 0 } else { groups.len() - 1 });
            integrals[g] += 0.5 * (b - a) * (y(a) + y(b));
        }
    }
    integrals
}

/// 出射分布を線形化するときの許容誤差
const DENSITY_TOLERANCE: f64 = 1.0e-4;

/// 入射エネルギーごとの出射分布を群に振り分けたもの。入射エネルギーの間は補間する。
#[derive(Debug, Clone)]
pub struct GroupedOutgoing {
    breaks: Vec<usize>,
    laws: Vec<Interpolation>,
    incident: Vec<f64>,
    probabilities: Vec<Vec<f64>>,
}

impl GroupedOutgoing {
    pub fn new(groups: &GroupStructure, distributions: &OutgoingDistributions) -> Self {
        Self {
            breaks: distributions.breaks.clone(),
            laws: distributions.laws.clone(),
            incident: distributions.incident.clone(),
            probabilities: distributions
                .spectra
                .iter()
                .map(|density| outgoing_group_integrals(groups, density))
                .collect(),
        }
    }

    /// 入射エネルギー
    pub fn incident(&self) -> &[f64] {
        &self.incident
    }

    /// 入射エネルギー e での、出射の群ごとの確率に `scale` を掛けて `out` に足す。
    /// 入射エネルギーの範囲外では、端の分布を使う。
    pub fn accumulate(&self, e: f64, scale: f64, out: &mut [f64]) {
        let last = self.incident.len() - 1;
        let add = |p: &[f64], factor: f64, out: &mut [f64]| {
            for (o, v) in out.iter_mut().zip(p) {
                *o += factor * v;
            }
        };
        if e <= self.incident[0] {
            return add(&self.probabilities[0], scale, out);
        }
        if e >= self.incident[last] {
            return add(&self.probabilities[last], scale, out);
        }
        let k = self.incident.partition_point(|&x| x <= e) - 1;
        let law = self.laws[self.breaks.partition_point(|&b| b < k + 2)];
        match law {
            Interpolation::Histogram => add(&self.probabilities[k], scale, out),
            _ => {
                let (x0, x1) = (self.incident[k], self.incident[k + 1]);
                let t = (e - x0) / (x1 - x0);
                add(&self.probabilities[k], scale * (1.0 - t), out);
                add(&self.probabilities[k + 1], scale * t, out);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integrates_piecewise_linear_density_by_group() {
        let groups = GroupStructure::new(vec![4.0, 2.0, 1.0]);
        let density = Tab1::linear(vec![1.0, 3.0], vec![0.5, 0.5]);
        assert_eq!(outgoing_group_integrals(&groups, &density), vec![0.5, 0.5]);
        // 群構造の外に出る分は端の群に入れる（NJOY と同じ。中性子を失わない）
        let density = Tab1::linear(vec![0.0, 3.0, 5.0], vec![0.2, 0.2, 0.2]);
        let r = outgoing_group_integrals(&groups, &density);
        assert!(
            (r[0] - 0.6).abs() < 1e-15 && (r[1] - 0.4).abs() < 1e-15,
            "{r:?}"
        );
    }

    #[test]
    fn linear_density_split_inside_a_segment() {
        let groups = GroupStructure::new(vec![2.0, 1.0, 0.5]);
        let density = Tab1::linear(vec![0.0, 2.0], vec![0.0, 1.0]);
        let r = outgoing_group_integrals(&groups, &density);
        // ∫ E/2 dE。0.5 eV 未満の分（0.0625）は最下位の群に入る
        assert!(
            (r[0] - 0.75).abs() < 1e-15 && (r[1] - 0.25).abs() < 1e-15,
            "{r:?}"
        );
    }

    #[test]
    fn interpolates_between_incident_energies() {
        let groups = GroupStructure::new(vec![4.0, 2.0, 1.0]);
        let distributions = OutgoingDistributions {
            breaks: vec![2],
            laws: vec![Interpolation::LinLin],
            incident: vec![10.0, 20.0],
            spectra: vec![
                Tab1::linear(vec![2.0, 4.0], vec![0.5, 0.5]),
                Tab1::linear(vec![1.0, 2.0], vec![1.0, 1.0]),
            ],
        };
        let grouped = GroupedOutgoing::new(&groups, &distributions);
        let mut out = vec![0.0; 2];
        grouped.accumulate(12.5, 2.0, &mut out);
        assert_eq!(out, vec![1.5, 0.5]);
        // 範囲外は端の分布
        let mut out = vec![0.0; 2];
        grouped.accumulate(30.0, 1.0, &mut out);
        assert_eq!(out, vec![0.0, 1.0]);
    }
}
