//! 出射エネルギー分布の群への振り分け。

use super::GroupStructure;
use super::integrate::GAUSS_LEGENDRE;
use crate::endf::{
    ClampedTab1, EnergyLaw, IncidentLaw, Interpolation, OutgoingDistributions, Tab1,
};
use crate::pointwise::linearize;
use crate::special::{gamma_p_2, gamma_p_3_2};

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
///
/// ふつうの補間則では、入射エネルギーごとに群に振り分けた値を補間する。単位基底の補間則では、
/// 隣り合う分布を出射エネルギーの範囲で正規化して補間した分布を、入射エネルギーごとに作ってから
/// 群に振り分ける（範囲が入射エネルギーとともに動くので、前もって群に振り分けておけない）。
#[derive(Debug, Clone)]
pub struct GroupedOutgoing {
    groups: GroupStructure,
    breaks: Vec<usize>,
    laws: Vec<IncidentLaw>,
    incident: Vec<f64>,
    /// [入射エネルギー] の出射分布（線形化したもの）
    spectra: Vec<Tab1>,
    /// [入射エネルギー][群] の出射の確率
    probabilities: Vec<Vec<f64>>,
    /// 重心系の分布なら標的の質量比（実験室系への変換に使う）
    center_of_mass: Option<f64>,
}

impl GroupedOutgoing {
    pub fn new(groups: &GroupStructure, distributions: &OutgoingDistributions) -> Self {
        Self::grouped(groups, distributions, None)
    }

    /// 重心系で与えられた分布を、重心系で等方として実験室系に変換して群に振り分ける。
    /// `awr` は標的の質量比。
    pub fn center_of_mass(
        groups: &GroupStructure,
        distributions: &OutgoingDistributions,
        awr: f64,
    ) -> Self {
        Self::grouped(groups, distributions, Some(awr))
    }

    fn grouped(
        groups: &GroupStructure,
        distributions: &OutgoingDistributions,
        center_of_mass: Option<f64>,
    ) -> Self {
        let spectra: Vec<Tab1> = distributions
            .spectra
            .iter()
            .map(|density| linearize(density, DENSITY_TOLERANCE))
            .collect();
        let mut grouped = Self {
            groups: groups.clone(),
            breaks: distributions.breaks.clone(),
            laws: distributions.laws.clone(),
            incident: distributions.incident.clone(),
            spectra,
            probabilities: Vec::new(),
            center_of_mass,
        };
        grouped.probabilities = grouped
            .incident
            .iter()
            .zip(&grouped.spectra)
            .map(|(&e, density)| grouped.group_integrals(e, density))
            .collect();
        grouped
    }

    /// 入射エネルギー e での出射分布（線形化したもの）を群ごとに積分する。
    fn group_integrals(&self, e: f64, density: &Tab1) -> Vec<f64> {
        match self.center_of_mass {
            None => outgoing_group_integrals(&self.groups, density),
            Some(awr) => laboratory_group_integrals(&self.groups, density, e / (awr + 1.0).powi(2)),
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
        let add = |p: &[f64], out: &mut [f64]| {
            for (o, v) in out.iter_mut().zip(p) {
                *o += scale * v;
            }
        };
        if e <= self.incident[0] {
            return add(&self.probabilities[0], out);
        }
        if e >= self.incident[last] {
            return add(&self.probabilities[last], out);
        }
        let k = self.incident.partition_point(|&x| x <= e) - 1;
        let IncidentLaw { law, unit_base } = self.laws[self.breaks.partition_point(|&b| b < k + 2)];
        let (x0, x1) = (self.incident[k], self.incident[k + 1]);
        if law == Interpolation::Histogram {
            add(&self.probabilities[k], out);
        } else if unit_base {
            let density =
                unit_base_density(law, e, (x0, &self.spectra[k]), (x1, &self.spectra[k + 1]));
            add(&self.group_integrals(e, &density), out);
        } else {
            let p: Vec<f64> = self.probabilities[k]
                .iter()
                .zip(&self.probabilities[k + 1])
                .map(|(&p0, &p1)| law.interpolate(e, (x0, p0), (x1, p1)))
                .collect();
            add(&p, out);
        }
    }
}

/// 入射エネルギー x0 と x1 の線形化した出射分布 `lo`, `hi` から、単位基底の補間則で入射エネルギー e の
/// 分布を作る。
///
/// 各分布の出射エネルギーの範囲 [a, b] を 0〜1 の座標 u に移し、範囲の両端と、u での密度 ×
/// 範囲の幅（u についての密度）を、補間則 `law` で入射エネルギーについて補間する。両方の分布の点を
/// u で合わせた格子の上では、どちらも u について線形なので、線形の補間則なら補間した密度も線形になる。
/// 範囲の幅が正であることは、読み込みのとき（`check_unit_base`）に確かめている。
fn unit_base_density(
    law: Interpolation,
    e: f64,
    (x0, lo): (f64, &Tab1),
    (x1, hi): (f64, &Tab1),
) -> Tab1 {
    let range = |t: &Tab1| (t.x[0], t.x[t.x.len() - 1]);
    let ((a0, b0), (a1, b1)) = (range(lo), range(hi));
    let a = law.interpolate(e, (x0, a0), (x1, a1));
    let b = law.interpolate(e, (x0, b0), (x1, b1));
    let mut u: Vec<f64> =
        lo.x.iter()
            .map(|x| (x - a0) / (b0 - a0))
            .chain(hi.x.iter().map(|x| (x - a1) / (b1 - a1)))
            .map(|u| u.clamp(0.0, 1.0))
            .collect();
    u.sort_by(f64::total_cmp);
    u.dedup();
    let y: Vec<f64> = u
        .iter()
        .map(|&u| {
            let p0 = linear_at(lo, a0 + u * (b0 - a0)) * (b0 - a0);
            let p1 = linear_at(hi, a1 + u * (b1 - a1)) * (b1 - a1);
            law.interpolate(e, (x0, p0), (x1, p1)) / (b - a)
        })
        .collect();
    Tab1::linear(u.iter().map(|u| a + u * (b - a)).collect(), y)
}

/// 線形化した表の、範囲 [x の最初, x の最後] の中の値（2点以上の表）。
fn linear_at(t: &Tab1, x: f64) -> f64 {
    let k =
        t.x.partition_point(|&xi| xi <= x)
            .saturating_sub(1)
            .min(t.x.len() - 2);
    let (xa, xb, ya, yb) = (t.x[k], t.x[k + 1], t.y[k], t.y[k + 1]);
    if xb == xa {
        yb
    } else {
        ya + (yb - ya) * (x - xa) / (xb - xa)
    }
}

/// MF5 の部分分布を群に振り分けたもの。
#[derive(Debug, Clone)]
pub enum GroupedLaw {
    Table(GroupedOutgoing),
    Analytic(AnalyticSpectrum),
}

/// 解析式の出射分布（LF=7, 9）。入射エネルギーごとに、群ごとの割合を不完全ガンマ関数の差で
/// 厳密に求める（NJOY GROUPR の anased と同じ）。
#[derive(Debug, Clone)]
pub struct AnalyticSpectrum {
    /// x = E'/θ 未満の出射の割合（正則化した下側不完全ガンマ関数）
    fraction: fn(f64) -> f64,
    restriction: f64,
    theta: ClampedTab1,
    bounds: Vec<f64>,
}

impl AnalyticSpectrum {
    fn accumulate(&self, e: f64, scale: f64, out: &mut [f64]) {
        let top = e - self.restriction;
        if top <= 0.0 {
            return;
        }
        // 入射エネルギーの範囲外では、端の θ を使う（表形式の分布と同じ）
        let theta = self.theta.eval(e);
        let total = (self.fraction)(top / theta);
        // 最下位の群から積み上げる。群構造の外に出る分は端の群に入れる
        let mut below = 0.0;
        for g in (0..out.len()).rev() {
            let upper = if g == 0 { top } else { self.bounds[g].min(top) };
            let cumulative = (self.fraction)(upper / theta);
            out[g] += scale * (cumulative - below) / total;
            below = cumulative;
            if upper >= top {
                break;
            }
        }
    }
}

impl GroupedLaw {
    pub fn new(groups: &GroupStructure, law: &EnergyLaw) -> Self {
        let analytic = |fraction: fn(f64) -> f64, restriction: f64, theta: &ClampedTab1| {
            Self::Analytic(AnalyticSpectrum {
                fraction,
                restriction,
                theta: theta.clone(),
                bounds: groups.bounds().to_vec(),
            })
        };
        match law {
            EnergyLaw::Tabulated(distributions) => {
                Self::Table(GroupedOutgoing::new(groups, distributions))
            }
            EnergyLaw::Maxwellian { restriction, theta } => {
                analytic(gamma_p_3_2, *restriction, theta)
            }
            EnergyLaw::Evaporation { restriction, theta } => {
                analytic(gamma_p_2, *restriction, theta)
            }
        }
    }

    /// 入射エネルギーについての積分で区間を分ける点
    pub fn breakpoints(&self) -> Vec<f64> {
        match self {
            Self::Table(grouped) => grouped.incident().to_vec(),
            // θ の表の点と、E - U が 0 や群境界を横切る点（群ごとの割合の傾きが変わる）
            Self::Analytic(spectrum) => spectrum
                .theta
                .table()
                .x
                .iter()
                .copied()
                .chain(std::iter::once(spectrum.restriction))
                .chain(spectrum.bounds.iter().map(|b| b + spectrum.restriction))
                .filter(|&e| e > 0.0)
                .collect(),
        }
    }

    /// 入射エネルギー e での、出射の群ごとの確率に `scale` を掛けて `out` に足す。
    pub fn accumulate(&self, e: f64, scale: f64, out: &mut [f64]) {
        match self {
            Self::Table(grouped) => grouped.accumulate(e, scale, out),
            Self::Analytic(spectrum) => spectrum.accumulate(e, scale, out),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::endf::IncidentLaw;

    fn assert_close(actual: &[f64], expected: &[f64]) {
        for (a, e) in actual.iter().zip(expected) {
            assert!((a - e).abs() <= 1e-14, "{actual:?} vs {expected:?}");
        }
    }

    #[test]
    fn evaporation_spectrum_is_cut_at_e_minus_u() {
        // E' ≤ E - U = 2 MeV。群ごとの割合は P(2, x) の差（x = E'/θ）
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e5, 1.0e-5]).unwrap();
        let law = EnergyLaw::Evaporation {
            restriction: 1.0e6,
            theta: ClampedTab1::new(Tab1::linear(vec![1.0e6, 2.0e7], vec![1.0e5, 1.0e5])).unwrap(),
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
        let groups = GroupStructure::new(vec![2.0e7, 1.0e6, 1.0e5, 1.0e-5]).unwrap();
        let law = EnergyLaw::Maxwellian {
            restriction: -3.0e7,
            theta: ClampedTab1::new(Tab1::linear(vec![1.0e6, 3.0e6], vec![1.0e5, 3.0e5])).unwrap(),
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
        let groups = GroupStructure::new(vec![4.0, 2.0, 1.0]).unwrap();
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
        let groups = GroupStructure::new(vec![2.0, 1.0, 0.5]).unwrap();
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
            laws: vec![IncidentLaw {
                law: Interpolation::LinLin,
                unit_base: false,
            }],
            incident: vec![e],
            spectra: vec![Tab1 {
                breaks: vec![2],
                laws: vec![Interpolation::Histogram],
                x: vec![1.0e6 - 1.0, 1.0e6 + 1.0],
                y: vec![0.5, 0.5],
            }],
        };
        let groups = GroupStructure::new(vec![2.0e6, 1.0e6, 1.0e5]).unwrap();
        let grouped = GroupedOutgoing::center_of_mass(&groups, &distributions, awr);
        let mut out = vec![0.0; 2];
        grouped.accumulate(e, 1.0, &mut out);
        // 実験室系の範囲は [562500, 1562500]
        assert!(
            (out[0] - 0.5625).abs() < 1e-4 && (out[1] - 0.4375).abs() < 1e-4,
            "{out:?}"
        );
    }

    /// 入射エネルギー 10 と 20 で、[0, 1] と [0, 3] の一様分布
    fn widening(law: IncidentLaw) -> OutgoingDistributions {
        OutgoingDistributions {
            breaks: vec![2],
            laws: vec![law],
            incident: vec![10.0, 20.0],
            spectra: vec![
                Tab1::linear(vec![0.0, 1.0], vec![1.0, 1.0]),
                Tab1::linear(vec![0.0, 3.0], vec![1.0 / 3.0, 1.0 / 3.0]),
            ],
        }
    }

    #[test]
    fn unit_base_keeps_the_shape_between_incident_energies() {
        // 単位基底では、中間の入射エネルギーで範囲が [0, 2] に伸びた一様分布になる
        let groups = GroupStructure::new(vec![3.0, 2.0, 1.0, 1.0e-5]).unwrap();
        let unit_base = IncidentLaw {
            law: Interpolation::LinLin,
            unit_base: true,
        };
        let grouped = GroupedOutgoing::new(&groups, &widening(unit_base));
        let mut out = vec![0.0; 3];
        grouped.accumulate(15.0, 1.0, &mut out);
        let expected = [0.0, 0.5, 0.5];
        for (o, e) in out.iter().zip(expected) {
            assert!((o - e).abs() < 1e-6, "{out:?}");
        }
        // ふつうの線形補間なら、2つの分布の平均（[0, 1] に 2/3、[1, 3] に 1/3）
        let cartesian = IncidentLaw {
            law: Interpolation::LinLin,
            unit_base: false,
        };
        let grouped = GroupedOutgoing::new(&groups, &widening(cartesian));
        let mut out = vec![0.0; 3];
        grouped.accumulate(15.0, 1.0, &mut out);
        let expected = [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0];
        for (o, e) in out.iter().zip(expected) {
            assert!((o - e).abs() < 1e-6, "{out:?}");
        }
    }

    #[test]
    fn unit_base_keeps_histogram_steps() {
        // 出射分布がヒストグラム（段差がある）でも、確率の合計と形を保つ。
        // [0, 0.5] に密度 2 と [0, 1.5] に密度 2/3 の中間は、[0, 1] に密度 1
        let step = |width: f64| Tab1 {
            breaks: vec![3],
            laws: vec![Interpolation::Histogram],
            x: vec![0.0, 0.5 * width, width],
            y: vec![2.0 / width, 0.0, 0.0],
        };
        let distributions = OutgoingDistributions {
            breaks: vec![2],
            laws: vec![IncidentLaw {
                law: Interpolation::LinLin,
                unit_base: true,
            }],
            incident: vec![10.0, 20.0],
            spectra: vec![step(1.0), step(3.0)],
        };
        let groups = GroupStructure::new(vec![3.0, 2.0, 1.0, 1.0e-5]).unwrap();
        let grouped = GroupedOutgoing::new(&groups, &distributions);
        let mut out = vec![0.0; 3];
        grouped.accumulate(15.0, 1.0, &mut out);
        let expected = [0.0, 0.0, 1.0];
        for (o, e) in out.iter().zip(expected) {
            assert!((o - e).abs() < 1e-9, "{out:?}");
        }
    }

    #[test]
    fn unit_base_in_the_center_of_mass_matches_the_laboratory_for_a_heavy_target() {
        // 標的が十分重ければ、重心系と実験室系はほぼ同じ
        let groups = GroupStructure::new(vec![3.0, 2.0, 1.0, 1.0e-5]).unwrap();
        let law = IncidentLaw {
            law: Interpolation::LinLin,
            unit_base: true,
        };
        let lab = GroupedOutgoing::new(&groups, &widening(law));
        let cm = GroupedOutgoing::center_of_mass(&groups, &widening(law), 1.0e9);
        let (mut a, mut b) = (vec![0.0; 3], vec![0.0; 3]);
        lab.accumulate(15.0, 1.0, &mut a);
        cm.accumulate(15.0, 1.0, &mut b);
        for (x, y) in a.iter().zip(&b) {
            assert!((x - y).abs() < 1e-4, "{a:?} vs {b:?}");
        }
    }

    #[test]
    fn incident_interpolation_follows_its_law() {
        // 線形・対数（LinLog）なら、入射エネルギーの対数の中点で2つの分布の平均
        let groups = GroupStructure::new(vec![3.0, 2.0, 1.0, 1.0e-5]).unwrap();
        let law = IncidentLaw {
            law: Interpolation::LinLog,
            unit_base: false,
        };
        let grouped = GroupedOutgoing::new(&groups, &widening(law));
        let mut out = vec![0.0; 3];
        grouped.accumulate((10.0_f64 * 20.0).sqrt(), 1.0, &mut out);
        let expected = [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0];
        for (o, e) in out.iter().zip(expected) {
            assert!((o - e).abs() < 1e-12, "{out:?}");
        }
    }

    #[test]
    fn interpolates_between_incident_energies() {
        let groups = GroupStructure::new(vec![4.0, 2.0, 1.0]).unwrap();
        let distributions = OutgoingDistributions {
            breaks: vec![2],
            laws: vec![IncidentLaw {
                law: Interpolation::LinLin,
                unit_base: false,
            }],
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
