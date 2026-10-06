//! 出射エネルギー分布の群への振り分け。

use super::GroupStructure;
use super::integrate::GAUSS_LEGENDRE;
use crate::endf::{EnergyLaw, Interpolation, OutgoingDistributions, Tab1};
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

/// 重心系の出射エネルギーの確率密度を、実験室系の群ごとに積分する。
///
/// 重心系でエネルギー E'cm の中性子は、重心系で等方なら、実験室系で
/// [(√E'cm - √Ec)², (√E'cm + √Ec)²] に一様に分布する。Ec = E/(A+1)² は重心の運動による
/// 中性子1個あたりのエネルギー。
fn laboratory_group_integrals(groups: &GroupStructure, density: &Tab1, ec: f64) -> Vec<f64> {
    const PIECES: usize = 16;
    let linear = linearize(density, DENSITY_TOLERANCE);
    let mut integrals = vec![0.0; groups.len()];
    for k in 0..linear.x.len().saturating_sub(1) {
        let (x0, x1) = (linear.x[k], linear.x[k + 1]);
        if x1 <= x0 {
            continue;
        }
        let (y0, y1) = (linear.y[k], linear.y[k + 1]);
        let width = (x1 - x0) / PIECES as f64;
        for i in 0..PIECES {
            let (mid, half) = (x0 + (i as f64 + 0.5) * width, 0.5 * width);
            for &(t, wk) in &GAUSS_LEGENDRE {
                let e = mid + half * t;
                let density = y0 + (y1 - y0) * (e - x0) / (x1 - x0);
                let (lo, hi) = (
                    (e.sqrt() - ec.sqrt()).powi(2),
                    (e.sqrt() + ec.sqrt()).powi(2),
                );
                uniform_fractions(groups, lo, hi, density * wk * half, &mut integrals);
            }
        }
    }
    integrals
}

/// [lo, hi] の一様分布を群に振り分け、`scale` を掛けて `out` に足す。
/// 群構造の外に出る分は端の群に入れる。
pub(super) fn uniform_fractions(
    groups: &GroupStructure,
    lo: f64,
    hi: f64,
    scale: f64,
    out: &mut [f64],
) {
    let bounds = groups.bounds();
    let edge = |e: f64| {
        groups
            .find(e)
            .unwrap_or(if e > bounds[0] { 0 } else { groups.len() - 1 })
    };
    if hi <= lo {
        out[edge(lo)] += scale;
        return;
    }
    let mut cuts: Vec<f64> = bounds
        .iter()
        .copied()
        .filter(|&b| b > lo && b < hi)
        .collect();
    cuts.push(lo);
    cuts.push(hi);
    cuts.sort_by(f64::total_cmp);
    for w in cuts.windows(2) {
        out[edge(0.5 * (w[0] + w[1]))] += scale * (w[1] - w[0]) / (hi - lo);
    }
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

    /// 重心系で与えられた分布を、重心系で等方として実験室系に変換して群に振り分ける。
    /// `awr` は標的の質量比。
    pub fn center_of_mass(
        groups: &GroupStructure,
        distributions: &OutgoingDistributions,
        awr: f64,
    ) -> Self {
        Self {
            breaks: distributions.breaks.clone(),
            laws: distributions.laws.clone(),
            incident: distributions.incident.clone(),
            probabilities: distributions
                .incident
                .iter()
                .zip(&distributions.spectra)
                .map(|(&e, density)| {
                    laboratory_group_integrals(groups, density, e / (awr + 1.0).powi(2))
                })
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

/// MF5 の部分分布を群に振り分けたもの。
#[derive(Debug, Clone)]
pub enum GroupedLaw {
    Table(GroupedOutgoing),
    Analytic(AnalyticSpectrum),
}

/// 解析式の出射分布（LF=7, 9）。
#[derive(Debug, Clone)]
pub struct AnalyticSpectrum {}

impl GroupedLaw {
    pub fn new(groups: &GroupStructure, law: &EnergyLaw) -> Self {
        todo!("{groups:?} {law:?}")
    }

    /// 入射エネルギーについての積分で区間を分ける点
    pub fn breakpoints(&self) -> Vec<f64> {
        todo!()
    }

    /// 入射エネルギー e での、出射の群ごとの確率に `scale` を掛けて `out` に足す。
    pub fn accumulate(&self, e: f64, scale: f64, out: &mut [f64]) {
        todo!("{e} {scale} {out:?}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::special::{gamma_p_2, gamma_p_3_2};

    fn assert_close(actual: &[f64], expected: &[f64]) {
        for (a, e) in actual.iter().zip(expected) {
            assert!((a - e).abs() <= 1e-14, "{actual:?} vs {expected:?}");
        }
    }

    #[test]
    fn evaporation_spectrum_is_cut_at_e_minus_u() {
        // E' ≤ E - U = 2 MeV。群ごとの割合は P(2, x) の差（x = E'/θ）
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e5, 1.0e-5]);
        let law = EnergyLaw::Evaporation {
            restriction: 1.0e6,
            theta: Tab1::linear(vec![1.0e6, 2.0e7], vec![1.0e5, 1.0e5]),
        };
        let grouped = GroupedLaw::new(&groups, &law);
        let mut out = vec![0.0; 3];
        grouped.accumulate(3.0e6, 2.0, &mut out);
        let p = |x: f64| gamma_p_2(x) / gamma_p_2(20.0);
        assert_close(
            &out,
            &[
                2.0 * (1.0 - p(10.0)),
                2.0 * (p(10.0) - p(1.0)),
                2.0 * p(1.0),
            ],
        );
        // E ≤ U では中性子を出さない
        let mut out = vec![0.0; 3];
        grouped.accumulate(1.0e6, 1.0, &mut out);
        assert_eq!(out, vec![0.0; 3]);
        // E - U が群境界を横切る点で、群ごとの割合の傾きが変わる
        let breakpoints = grouped.breakpoints();
        for e in [1.0e6, 1.1e6, 2.0e6] {
            assert!(breakpoints.contains(&e), "{e}: {breakpoints:?}");
        }
    }

    #[test]
    fn maxwellian_interpolates_theta() {
        // θ(2 MeV) = 200 keV。E - U が群構造の上限を超える分は最上位の群に入る
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e5, 1.0e-5]);
        let law = EnergyLaw::Maxwellian {
            restriction: -3.0e7,
            theta: Tab1::linear(vec![1.0e6, 3.0e6], vec![1.0e5, 3.0e5]),
        };
        let grouped = GroupedLaw::new(&groups, &law);
        let mut out = vec![0.0; 3];
        grouped.accumulate(2.0e6, 1.0, &mut out);
        let p = |x: f64| gamma_p_3_2(x) / gamma_p_3_2(160.0);
        assert_close(&out, &[1.0 - p(5.0), p(5.0) - p(0.5), p(0.5)]);
        assert!(grouped.breakpoints().contains(&3.0e6));
    }

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
    fn center_of_mass_energy_spreads_in_laboratory() {
        // 重心系で E'cm にほぼ単色の中性子は、実験室系で
        // [(√E'cm - √Ec)², (√E'cm + √Ec)²] に一様に分布する（Ec = E/(A+1)²）
        let (e, awr) = (1.6e7, 15.0);
        let ec = e / (awr + 1.0_f64).powi(2);
        assert_eq!(ec, 62500.0);
        let distributions = OutgoingDistributions {
            breaks: vec![1],
            laws: vec![Interpolation::LinLin],
            incident: vec![e],
            spectra: vec![Tab1 {
                breaks: vec![2],
                laws: vec![Interpolation::Histogram],
                x: vec![1.0e6 - 1.0, 1.0e6 + 1.0],
                y: vec![0.5, 0.5],
            }],
        };
        let groups = GroupStructure::new(vec![2.0e6, 1.0e6, 1.0e5]);
        let grouped = GroupedOutgoing::center_of_mass(&groups, &distributions, awr);
        let mut out = vec![0.0; 2];
        grouped.accumulate(e, 1.0, &mut out);
        // 実験室系の範囲は [562500, 1562500]
        assert!(
            (out[0] - 0.5625).abs() < 1e-4 && (out[1] - 0.4375).abs() < 1e-4,
            "{out:?}"
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
