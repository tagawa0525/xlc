//! 出射エネルギー分布の群への振り分け。

use crate::endf::{
    ClampedTab1, EnergyLaw, IncidentLaw, Interpolation, OutgoingDistributions, Tab1,
};
use crate::pointwise::linearize;
use crate::special::{erf, exp_e1, gamma_p_2, gamma_p_3_2, lower_gamma_series};
use xlc_format::group::GroupStructure;
use xlc_numeric::quadrature::GAUSS_LEGENDRE_4;

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
            for &(t, wk) in &GAUSS_LEGENDRE_4 {
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
/// 範囲の幅（u についての密度）を、補間則 `law` で入射エネルギーについて補間する。両方の分布の区切りを
/// u で合わせた区間ごとに、区間の中点を含む線分で両端の値を求めて点を出す。段差（同じ x の2点）は
/// 区間の境目に左右の値として残る。範囲の幅が正であることは、読み込みのとき（`check_unit_base`）に
/// 確かめている。
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
    let mut cuts: Vec<f64> =
        lo.x.iter()
            .map(|x| (x - a0) / (b0 - a0))
            .chain(hi.x.iter().map(|x| (x - a1) / (b1 - a1)))
            .map(|u| u.clamp(0.0, 1.0))
            .collect();
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    let (mut x, mut y) = (Vec::new(), Vec::new());
    for w in cuts.windows(2) {
        let (ua, ub) = (w[0], w[1]);
        let middle = 0.5 * (ua + ub);
        // u についての密度を、区間の中点を含む線分で両端について求める
        let density = |t: &Tab1, start: f64, width: f64, u: f64| {
            segment_value(t, start + middle * width, start + u * width) * width
        };
        for u in [ua, ub] {
            let p0 = density(lo, a0, b0 - a0, u);
            let p1 = density(hi, a1, b1 - a1, u);
            x.push(a + u * (b - a));
            y.push(law.interpolate(e, (x0, p0), (x1, p1)) / (b - a));
        }
    }
    Tab1::linear(x, y)
}

/// 線形化した表（2点以上）の、点 `inside` を含む線分を延ばした x での値。
fn segment_value(t: &Tab1, inside: f64, x: f64) -> f64 {
    let k =
        t.x.partition_point(|&xi| xi <= inside)
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

/// θ(E) で決まる出射分布（LF=5, 7, 9, 11, 12）。入射エネルギーごとに、群ごとの割合を x = E'/θ までの累積の差で
/// 厳密に求める（LF=7, 9 は不完全ガンマ関数、LF=11 は誤差関数、LF=12 は閉じた式。累積の差で分けるのは NJOY GROUPR の anased と
/// 同じ）。LF=11 の θ は a(E)。LF=11 の √(E'/a) < 0.01 の累積は、anased が x^(3/2) の初項だけで近似するのに対し、ここでは
/// 正の級数で求める（`watt_small_cumulative`）。初項は級数の打ち切りで、相対の差は (0.6 - 0.4 a b/4) x ほど（U-233 の a b/4 = 0.62 で、
/// √x = 0.01 の 3.5e-5）。近似の差なので、NJOY の誤りとして warn は出さない。
#[derive(Debug, Clone)]
pub struct AnalyticSpectrum {
    /// x = E'/θ 未満の出射の割合
    fraction: Cumulative,
    /// E' ≤ E - U の U。`None` は切らない（LF=12）
    restriction: Option<f64>,
    theta: Theta,
    bounds: Vec<f64>,
}

/// 分布の温度 θ。LF=5 はエネルギーによらない値（`EnergyLaw::GeneralEvaporation`）。
#[derive(Debug, Clone)]
enum Theta {
    Constant(f64),
    /// 入射エネルギーの範囲外では端の θ を使う（表形式の分布と同じ）
    Table(ClampedTab1),
}

impl Theta {
    fn eval(&self, e: f64) -> f64 {
        match self {
            Self::Constant(theta) => *theta,
            Self::Table(theta) => theta.eval(e),
        }
    }
}

/// x = E'/θ 未満の出射の（正規化していない）累積。
#[derive(Debug, Clone)]
enum Cumulative {
    /// 正則化した下側不完全ガンマ関数（LF=7, 9）
    Gamma(fn(f64) -> f64),
    /// 表の g(x)（LF=5）を線形化した区分線形の密度の積分。`integral[k]` は x[k] までの積分
    Table { shape: Tab1, integral: Vec<f64> },
    /// Watt の分布（LF=11）の累積（`watt_cumulative`）。θ は a(E)、`b` は b(E)
    Watt { b: ClampedTab1 },
    /// Madland-Nix のスペクトル（LF=12）の 0〜∞ で 1 の累積。θ は T_M、`light` と `heavy` は E_F (eV)
    MadlandNix { light: f64, heavy: f64 },
}

impl Cumulative {
    fn table(shape: &Tab1) -> Self {
        let shape = linearize(shape, DENSITY_TOLERANCE);
        let mut integral = vec![0.0; shape.x.len()];
        for k in 1..shape.x.len() {
            integral[k] = integral[k - 1]
                + 0.5 * (shape.x[k] - shape.x[k - 1]) * (shape.y[k] + shape.y[k - 1]);
        }
        Self::Table { shape, integral }
    }

    /// 入射エネルギー e（θ は `theta`）で決まり、x によらない累積の形の量（`at` に渡す）。LF=11 は a b/4、LF=12 は θ、
    /// ほかは使わない。群ごとに呼ぶ `at` の外で 1 度だけ求める（LF=11 の b(E) の補間を群ごとに繰り返さない）
    fn parameter(&self, e: f64, theta: f64) -> f64 {
        match self {
            Self::Watt { b } => 0.25 * theta * b.eval(e),
            _ => theta,
        }
    }

    /// x = E'/θ までの累積。`parameter` は x を出した入射エネルギーの `Cumulative::parameter`
    fn at(&self, x: f64, parameter: f64) -> f64 {
        match self {
            Self::Gamma(p) => p(x),
            Self::Watt { .. } => watt_cumulative(x, parameter),
            Self::MadlandNix { light, heavy } => {
                0.5 * (madland_nix_fraction(x, light / parameter)
                    + madland_nix_fraction(x, heavy / parameter))
            }
            Self::Table { shape, integral } => {
                let (xs, ys) = (&shape.x, &shape.y);
                let k = xs.partition_point(|&xk| xk <= x);
                if k == 0 {
                    return 0.0;
                }
                if k == xs.len() {
                    return integral[k - 1];
                }
                // xs[k - 1] ≤ x < xs[k]（段差の重複点は partition_point が越える）
                let (x0, x1, y0, y1) = (xs[k - 1], xs[k], ys[k - 1], ys[k]);
                let y = y0 + (y1 - y0) * (x - x0) / (x1 - x0);
                integral[k - 1] + 0.5 * (x - x0) * (y0 + y)
            }
        }
    }
}

/// x = E'/a がこれより小さいところは、閉じた式でなく級数で累積を求める（`watt_small_cumulative`）。閉じた式は 1 程度の項の
/// 差で、累積（x^(3/2) で 0 に近づく）の相対誤差が ε x^(-3/2) ほどになる。級数は項がすべて正で、x が大きいと項の数が増える
/// （`special` の不完全ガンマ関数と同じ切り替え）。x = 1 の前後の両者の相対差は 6e-15 以下（c = a b/4 = 0.05〜5）。数値積分との
/// 比較は c = 1e-6〜5（テスト）。読み込みが c ≤ 5 に限る（endf の mf5 の `WATT_MAX_C`）
const WATT_SERIES_LIMIT: f64 = 1.0;

/// ENDF-102 の 5.1.1.5 の exp(-E'/a) sinh(√(b E')) の、E'/a < x の（正規化しない）累積を a e^c で割ったもの。`c` は a b/4。
///
/// E' = a u² と置くと、r = √c として 2u e^(-u²) sinh(2ru) = u (e^(-(u-r)²) - e^(-(u+r)²)) e^c なので、
/// W(x) = (√π/2) r (erf(√x - r) + erf(√x + r)) - (e^(-(√x - r)²) - e^(-(√x + r)²))/2（ENDF-102 の正規化の定数 I を
/// a e^c で割ったものと同じ形）。小さい x（`WATT_SERIES_LIMIT`）は
/// `watt_small_cumulative`。
fn watt_cumulative(x: f64, c: f64) -> f64 {
    if x < WATT_SERIES_LIMIT {
        return watt_small_cumulative(x, c);
    }
    let (r, s) = (c.sqrt(), x.sqrt());
    0.5 * std::f64::consts::PI.sqrt() * r * (erf(s - r) + erf(s + r))
        - 0.5 * ((-(s - r) * (s - r)).exp() - (-(s + r) * (s + r)).exp())
}

/// x < `WATT_SERIES_LIMIT` での `watt_cumulative`。sinh(2r√s) を級数に展開して項ごとに積分すると
/// W(x) = e^(-c) Σ_k (2r)^(2k+1)/(2k+1)! γ(k + 3/2, x) で、γ(p, x) = x^p e^(-x)/p Σ_n xⁿ/((p+1)…(p+n))
/// （`lower_gamma_series`）。どちらの和も項が正なので、桁落ちしない。
fn watt_small_cumulative(x: f64, c: f64) -> f64 {
    // k = 0 の (2r)^(2k+1) x^(k+3/2)/((2k+1)! (k+3/2))
    let mut term = 2.0 * c.sqrt() * x * x.sqrt() / 1.5;
    let mut sum = 0.0;
    for k in 0.. {
        let p = f64::from(k) + 1.5;
        let part = term * lower_gamma_series(p, x);
        sum += part;
        if part <= f64::EPSILON * sum {
            break;
        }
        let k = f64::from(k);
        term *= 4.0 * c * x / ((2.0 * k + 2.0) * (2.0 * k + 3.0)) * p / (p + 1.0);
    }
    (-(c + x)).exp() * sum
}

/// E'/T_M が E_F/T_M のこの倍より小さいところは、閉じた式でなく級数で累積を求める（`madland_nix_small_fraction`）。
/// 閉じた式は 1 程度の項の差で、累積（x^(3/2) で 0 に近づく）に 1e-16 ほどの絶対誤差が残る。この倍の x では、E_F/T_M = 0.5 で
/// 累積が 2e-7 ほどになる。この x の前後の相対誤差は、閉じた式で 3e-9 以下、級数で 1.1e-11 以下（ε = 0.05〜5 で、数値積分と
/// 比べた）。
const MADLAND_NIX_SMALL: f64 = 1.0e-4;

/// ENDF-102 の 5.1.1.6 の g(E', E_F) の、E'/T_M < x の割合（0〜∞ で 1）。`fragment` は E_F/T_M。
///
/// NJOY GROUPR の anased の LF=12 の 3 つの領域（ENDF-102 の 5.4.10）の式を、1 つの原始関数 H(√E') にまとめたもの:
/// A = (√E' + √E_F)²/T_M、B = (√E' - √E_F)²/T_M、σ = sgn(√E' - √E_F) として、領域 I と II の式は H(b) - H(a) で、
/// B の項の符号が σ で替わる。領域 III の定数 -2 は、σ e^(-B)(1 + B) の項を σ (e^(-B)(1 + B) - 1) に替えて H を
/// E' = E_F で連続にしたもの。H(∞) = 1/2、H(0) = -1/2 なので、割合は H + 1/2。小さい x（`MADLAND_NIX_SMALL`）は
/// `madland_nix_small_fraction`。式の正しさは、ENDF-102 の g の数値積分との比較（テスト）で確かめる。
fn madland_nix_fraction(x: f64, fragment: f64) -> f64 {
    if x == f64::INFINITY {
        return 1.0;
    }
    if x < MADLAND_NIX_SMALL * fragment {
        return madland_nix_small_fraction(x, fragment);
    }
    let r = 1.0 / fragment.sqrt();
    let signed = x.sqrt() - fragment.sqrt();
    let sigma = signed.signum();
    let a = (x.sqrt() + fragment.sqrt()).powi(2);
    let b = signed * signed;
    let gamma_3_2 = |u: f64| 0.5 * std::f64::consts::PI.sqrt() * gamma_p_3_2(u);
    let (gamma_a, gamma_b) = (gamma_3_2(a), gamma_3_2(b));
    // γ(5/2, u) = 3/2 γ(3/2, u) - u^(3/2) e^(-u)
    let gamma_5_2 = |u: f64, gamma: f64| 1.5 * gamma - u.powf(1.5) * (-u).exp();
    let decay = |u: f64| (-u).exp() * (1.0 + u);
    // u^n E1(u) は u → 0 で 0（n ≥ 2）
    let e1_b = if b == 0.0 { 0.0 } else { madland_nix_e1(b) };
    let upper = (2.0 / 15.0 * r * a.powf(2.5) - a * a / 6.0) * madland_nix_e1(a)
        + (r / 3.0 * a - 2.0 / 3.0 * a.sqrt()) * gamma_a
        - r / 5.0 * gamma_5_2(a, gamma_a)
        - 0.5 * decay(a);
    let lower = (2.0 / 15.0 * r * b.powf(2.5) + sigma * b * b / 6.0) * e1_b
        + (r / 3.0 * b + sigma * 2.0 / 3.0 * b.sqrt()) * gamma_b
        - r / 5.0 * gamma_5_2(b, gamma_b)
        + 0.5 * sigma * (decay(b) - 1.0);
    0.5 + upper - lower
}

/// x < `MADLAND_NIX_SMALL` E_F/T_M での `madland_nix_fraction`。項が打ち消し合わない級数で求める（T_M を単位とし、
/// ε = E_F/T_M、w = √ε、X = √x）。φ(u) = u^(3/2) E1(u) + γ(3/2, u) の導関数は (3/2) √u E1(u) なので、√u = w + t と置くと
/// g(s²) = (1/(2w)) ∫_{-s}^{s} q(w + t) dt、q(v) = 2v² E1(v²)。t について展開して E' = s² で積分すると、累積は
/// (2/w) Σ_k q^(2k)(w) X^(2k+3)/((2k+1)! (2k+3)) で、k = 0〜2 の項を使う:
/// q(w) = 2ε E1(ε)、q''(w) = 4 E1(ε) + e^(-ε)(8ε - 12)、q''''(w) = e^(-ε)(8/ε + 56 - 128ε + 32ε²)。
fn madland_nix_small_fraction(x: f64, fragment: f64) -> f64 {
    let e1 = madland_nix_e1(fragment);
    let decay = (-fragment).exp();
    let q0 = 2.0 * fragment * e1;
    let q2 = 4.0 * e1 + decay * (8.0 * fragment - 12.0);
    let q4 = decay * (8.0 / fragment + 56.0 - 128.0 * fragment + 32.0 * fragment * fragment);
    2.0 / fragment.sqrt() * x.sqrt().powi(3) * (q0 / 3.0 + q2 * x / 30.0 + q4 * x * x / 840.0)
}

/// E1(u)（u > 0）。u は正の有限の値で、そこで exp_e1 は収束する。万一の失敗は NaN にして、χ を有限の値に限るライブラリの
/// 検査（xlc_format の NuclideData の check）で止める（`accumulate` はエラーを返さない）
fn madland_nix_e1(u: f64) -> f64 {
    exp_e1(u).map_or(f64::NAN, |scaled| scaled * (-u).exp())
}

impl AnalyticSpectrum {
    fn accumulate(&self, e: f64, scale: f64, out: &mut [f64]) {
        let top = match self.restriction {
            Some(restriction) => e - restriction,
            None => f64::INFINITY,
        };
        if top <= 0.0 {
            return;
        }
        let theta = self.theta.eval(e);
        let parameter = self.fraction.parameter(e, theta);
        let total = self.fraction.at(top / theta, parameter);
        // 最下位の群から積み上げる。群構造の外に出る分は端の群に入れる
        let mut below = 0.0;
        for g in (0..out.len()).rev() {
            let upper = if g == 0 { top } else { self.bounds[g].min(top) };
            let cumulative = self.fraction.at(upper / theta, parameter);
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
        let analytic = |fraction: Cumulative, restriction: Option<f64>, theta: Theta| {
            Self::Analytic(AnalyticSpectrum {
                fraction,
                restriction,
                theta,
                bounds: groups.bounds().to_vec(),
            })
        };
        match law {
            EnergyLaw::Tabulated(distributions) => {
                Self::Table(GroupedOutgoing::new(groups, distributions))
            }
            EnergyLaw::GeneralEvaporation {
                restriction,
                theta,
                shape,
            } => analytic(
                Cumulative::table(shape),
                Some(*restriction),
                Theta::Constant(*theta),
            ),
            EnergyLaw::Maxwellian { restriction, theta } => analytic(
                Cumulative::Gamma(gamma_p_3_2),
                Some(*restriction),
                Theta::Table(theta.clone()),
            ),
            EnergyLaw::Evaporation { restriction, theta } => analytic(
                Cumulative::Gamma(gamma_p_2),
                Some(*restriction),
                Theta::Table(theta.clone()),
            ),
            EnergyLaw::Watt { restriction, a, b } => analytic(
                Cumulative::Watt { b: b.clone() },
                Some(*restriction),
                Theta::Table(a.clone()),
            ),
            EnergyLaw::MadlandNix {
                light,
                heavy,
                temperature,
            } => analytic(
                Cumulative::MadlandNix {
                    light: *light,
                    heavy: *heavy,
                },
                None,
                Theta::Table(temperature.clone()),
            ),
        }
    }

    /// 入射エネルギーについての積分で区間を分ける点
    pub fn breakpoints(&self) -> Vec<f64> {
        match self {
            Self::Table(grouped) => grouped.incident().to_vec(),
            // θ の表の点と、E - U が 0 や群境界を横切る点（群ごとの割合の傾きが変わる）。表の g(x)（θ は一定）では、
            // E - U が g の表の点の θ 倍を横切る点でも傾きが変わる。LF=11 は b の表の点でも傾きが変わる。U で切らない
            // 分布（LF=12）は θ の表の点だけ
            Self::Analytic(spectrum) => {
                let mut points = Vec::new();
                if let Some(restriction) = spectrum.restriction {
                    points.push(restriction);
                    points.extend(spectrum.bounds.iter().map(|b| b + restriction));
                    if let (Cumulative::Table { shape, .. }, Theta::Constant(theta)) =
                        (&spectrum.fraction, &spectrum.theta)
                    {
                        points.extend(shape.x.iter().map(|x| restriction + theta * x));
                    }
                }
                if let Theta::Table(theta) = &spectrum.theta {
                    points.extend(&theta.table().x);
                }
                if let Cumulative::Watt { b } = &spectrum.fraction {
                    points.extend(&b.table().x);
                }
                points.retain(|&e| e > 0.0);
                points
            }
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
mod tests;
