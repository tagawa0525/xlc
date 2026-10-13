//! ドップラー拡がり（NJOY BROADR と同じ SIGMA1 法の核関数）。
//!
//! 0 K の lin-lin のポイントワイズの断面積を、標的核の熱運動（自由ガスの Maxwell 分布）で平均する。
//! 速度に比例する変数 x = √(αE)（α = A/(kT)）を使うと、拡げた断面積は
//!
//!   σ̄(y) = 1/(y²√π) ∫₀^∞ x² σ(x) [exp(-(x-y)²) - exp(-(x+y)²)] dx
//!
//! になる。表の最初の点より下は 1/v で x = 0 まで、最後の点より上は一定値で延ばす。
//!
//! NJOY は区間ごとの積分を誤差関数の差で解析的に求める。ここでは区間ごとに Gauss-Legendre 求積する。
//! 被積分関数は区間の中で滑らかなので、区間の幅に応じた少ない点数で誤差を 1e-11 以下に抑えられ、
//! その上限を誤差公式で示せる（`quadrature`）。NJOY の方法を移した計算とは 1e-7 程度で一致する。
//!
//! 格子の作り方は NJOY broadn と同じ規則による（`grid`、`nodes`、`split`）。

use crate::endf::{Mf2, RangeKind, Tab1};
use crate::physics::BK;
use crate::pointwise::{left_limit, right_limit};
use xlc_numeric::quadrature::{GAUSS_LEGENDRE_2, GAUSS_LEGENDRE_4, GAUSS_LEGENDRE_8};

/// 拡げる上限エネルギーの既定の最大値 (eV)。NJOY BROADR と同じ。
pub const MAX_LIMIT: f64 = 6.5e6;

/// 共鳴パラメータがない核種で、これより上から始まる反応をしきい反応とみなす (eV)。
const LOW_THRESHOLD: f64 = 1.0;

/// しきい反応があるとき、上限をしきい値のこの倍にする（NJOY と同じ）。
const BELOW_THRESHOLD: f64 = 0.99999;

/// |x - y| がこれを超える範囲は積分しない（核関数が exp(-25) ≈ 1.4e-11 以下になる）。
/// NJOY は 4 で打ち切るが、それでは低エネルギーで打ち切りの誤差が 2e-6 程度になる。
const WINDOW: f64 = 5.0;

// 格子の作り方の定数。どれも NJOY broadn と同じ値。

/// 収束判定の許容誤差 errmax、errint の、`tol`（errthn）に対する比（NJOY の既定値）
const ERRMAX_PER_TOL: f64 = 10.0;
const ERRINT_PER_TOL: f64 = 1.0 / 20_000.0;
/// これより下のエネルギー (eV) では、許容誤差を `LOW_ENERGY_SCALE` 倍にする
const LOW_ENERGY: f64 = 0.4999;
const LOW_ENERGY_SCALE: f64 = 0.2;
/// 節の間に挟む元の格子点の数の上限
const MAX_SKIPPED: usize = 10;
/// 節のエネルギーの、前の節に対する比の上限
const MAX_RATIO: f64 = 2.01;
/// 増減の向きを決めるとき、隣の点との差が値のこの割合未満なら平らとみなす
const FLAT: f64 = 1.0e-3;
/// 区間の両端の値の比がこれを超えると、中点の誤差によらず区間を分ける
const MAX_END_RATIO: f64 = 3.0;
/// 区間の幅が直前の区間のこの倍を超えると、収束とみなさない（誤った収束を防ぐ）
const MAX_GROWTH: f64 = 4.1;
/// 誤差の判定に足す絶対誤差
const ERRMIN: f64 = 1.0e-15;
/// 全反応の和に対する比がこれより小さい反応は、収束判定に使わない
const NEGLIGIBLE: f64 = 1.0e-6;

const SQRT_PI: f64 = 1.772_453_850_905_516;

/// 0 K の表を、温度 `temperature` (K) に拡げる。
///
/// `tables` は lin-lin の表で、同じエネルギー格子の上でまとめて拡げる。
/// `e_max` 未満のエネルギーは拡げた値を、NJOY broadn と同じ判定（`tol` を errthn とする）で
/// 細かくした格子に置き、`e_max` 以上は元の表の点をそのまま使う。
/// 最初の点が格子全体の最初の点より上にある表（しきい反応）は、しきい値より下を 0 のままにする。
pub fn broaden(tables: &[Tab1], awr: f64, temperature: f64, e_max: f64, tol: f64) -> Vec<Tab1> {
    let alpha = awr / (BK * temperature);
    let union = Union::new(tables);
    // 点が2つ未満なら拡げる区間がない
    let [first, .., last] = union.energies[..] else {
        return tables.to_vec();
    };
    let top = e_max.min(last);
    if top <= first {
        return tables.to_vec();
    }
    let integral = Integral {
        alpha,
        x: union.energies.iter().map(|e| (alpha * e).sqrt()).collect(),
        values: &union.values,
    };

    let points = grid(&union, &integral, top, tol);

    tables
        .iter()
        .enumerate()
        .map(|(i, t)| {
            let start = t.x[0];
            if start >= top {
                return t.clone();
            }
            let (mut x, mut y) = (Vec::new(), Vec::new());
            if start > first {
                // しきい反応: しきい値より下に漏れた分は捨てる
                x.push(start);
                y.push(t.y[0]);
            }
            for (e, v) in points.iter().filter(|(e, _)| *e > start || start == first) {
                x.push(*e);
                y.push(v[i]);
            }
            // e_max から上は元の表の点。e_max に段差があれば、右側の値から使う
            let mut k = t.x.partition_point(|&e| e < top);
            if t.x.get(k + 1) == Some(&top) {
                k += 1;
            }
            x.extend_from_slice(&t.x[k..]);
            y.extend_from_slice(&t.y[k..]);
            Tab1::linear(x, y)
        })
        .collect()
}

/// 拡げた値の格子（NJOY broadn と同じ作り方）。
///
/// 元の格子から選んだ節と上限 `top` で拡げた値を計算し、節の間は中点で収束を判定して分ける。
/// NJOY と違う点は次のとおり。
/// - エネルギーを有効数字 7 桁に丸めない。NJOY は丸めた節のエネルギーで積分の区間を探すため、
///   丸めで節より下になった点で1〜3区間ずれた値を出す（BROADR の出力との食い違いの原因）
/// - 分ける深さに上限（NJOY は 12）を設けない。NJOY の上限は固定長の配列のためのもので、
///   上限に達した区間は収束していなくてもそのまま残る。分割は、中点が両端の間に入らなくなれば止まる
fn grid(union: &Union, integral: &Integral, top: f64, tol: f64) -> Vec<(f64, Vec<f64>)> {
    let node = |k: usize| {
        let values = integral.at(union.energies[k]);
        let cold = &union.values[k];
        let matches_cold = values
            .iter()
            .zip(cold)
            .all(|(hot, cold)| (hot - cold).abs() <= tol * cold);
        Point {
            e: union.energies[k],
            values,
            node: Some(k),
            matches_cold,
        }
    };
    let selected = nodes(&union.energies, &union.values, top);
    let mut uppers = selected[1..].iter().map(|&k| node(k)).collect::<Vec<_>>();
    uppers.push(Point {
        e: top,
        values: integral.at(top),
        node: None,
        matches_cold: false,
    });

    let mut out: Vec<(f64, Vec<f64>)> = Vec::new();
    let mut lower = node(selected[0]);
    for upper in uppers {
        // stack の末尾が、lower と組む区間の上端
        let mut stack = vec![upper];
        while let Some(next) = stack.pop() {
            match split(&lower, &next, out.last(), integral, tol) {
                None => {
                    out.push((lower.e, lower.values));
                    lower = next;
                }
                // 分けるときは、上端を戻してから中点を積む
                Some(mid) => {
                    stack.push(next);
                    stack.push(mid);
                }
            }
        }
    }
    out.push((lower.e, lower.values));
    out
}

/// 格子の点。`node` は元の格子の節の番号（中点なら `None`）、
/// `matches_cold` は拡げた値が全反応で 0 K の値と `tol` 以内か。
struct Point {
    e: f64,
    values: Vec<f64>,
    node: Option<usize>,
    matches_cold: bool,
}

/// 区間 (lower, upper) を中点で分けるなら、その中点を返す（NJOY broadn の収束判定）。
/// `previous` は直前に格子に加えた点。
fn split(
    lower: &Point,
    upper: &Point,
    previous: Option<&(f64, Vec<f64>)>,
    integral: &Integral,
    tol: f64,
) -> Option<Point> {
    // 隣り合う元の格子点で、どちらも拡げても値が変わらなければ、そのまま線形補間でよい
    if let (Some(a), Some(b)) = (lower.node, upper.node)
        && b == a + 1
        && lower.matches_cold
        && upper.matches_cold
    {
        return None;
    }
    let em = 0.5 * (lower.e + upper.e);
    if em <= lower.e || em >= upper.e {
        return None;
    }
    let scale = if upper.e < LOW_ENERGY {
        LOW_ENERGY_SCALE
    } else {
        1.0
    };
    let (errthn, errmax) = (scale * tol, scale * ERRMAX_PER_TOL * tol);
    let errint = ERRINT_PER_TOL * tol;
    let dx = upper.e - lower.e;
    let sn = integral.at(em);
    let total: f64 = sn.iter().sum();
    let converged = (0..sn.len()).all(|i| {
        if total < ERRMIN || (sn[i] / total).abs() < NEGLIGIBLE {
            return true;
        }
        let (a, b) = (lower.values[i], upper.values[i]);
        if b > MAX_END_RATIO * a || b < a / MAX_END_RATIO {
            return false;
        }
        let dy = (sn[i] - 0.5 * (a + b)).abs();
        if dy <= errthn * sn[i].abs() + ERRMIN {
            return true;
        }
        dy <= errmax * sn[i].abs() + ERRMIN && dy * dx / 2.0 <= errint * em
    });
    // 直前の区間より極端に広い区間は、収束を取り違えている恐れがある
    let too_wide = previous.is_some_and(|(e, _)| dx > MAX_GROWTH * (lower.e - e));
    if converged && !too_wide {
        return None;
    }
    Some(Point {
        e: em,
        values: sn,
        node: None,
        matches_cold: false,
    })
}

/// 拡げた値を計算する元の格子の点（節）の番号。`top` 未満の点から NJOY broadn と同じ規則で選ぶ。
///
/// 最初の点と `top` 未満の最後の点のほか、次の点を節にする: いずれかの反応で増減の向きが
/// 前の節と変わる点、前の節から `MAX_SKIPPED` 点を超えた点、エネルギーが前の節の `MAX_RATIO` 倍を
/// 超えた点、有効数字3桁で表せるエネルギー（0.0253 eV を含む）。前の節と同じエネルギーの点は飛ばす。
fn nodes(energies: &[f64], values: &[Vec<f64>], top: f64) -> Vec<usize> {
    let below = energies.partition_point(|&e| e < top);
    if below == 0 {
        return Vec::new();
    }
    let last = below - 1;
    // 点 k から次の点への、反応ごとの増減の向き（平らは増加と同じ）
    let rising = |k: usize| -> Vec<bool> {
        values[k]
            .iter()
            .zip(&values[k + 1])
            .map(|(a, b)| b - a >= 0.0 || (b - a).abs() < FLAT * a.abs())
            .collect()
    };
    let mut nodes = vec![0];
    if last == 0 {
        return nodes;
    }
    let mut direction = rising(0);
    // 最後に選んだ節
    let mut previous = 0;
    for k in 1..=last {
        let e = energies[k];
        if e == energies[previous] {
            continue;
        }
        if k == last {
            nodes.push(k);
            break;
        }
        let here = rising(k);
        if k > previous + MAX_SKIPPED
            || e > MAX_RATIO * energies[previous]
            || is_round(e)
            || here != direction
        {
            nodes.push(k);
            previous = k;
            direction = here;
        }
    }
    nodes
}

/// 有効数字3桁で表せるエネルギーか。
fn is_round(e: f64) -> bool {
    let unit = 10f64.powi(e.log10().floor() as i32 - 2);
    let rounded = (e / unit).round() * unit;
    (e - rounded).abs() < 1.0e-9 * rounded
}

/// 複数の lin-lin の表を、共通のエネルギー格子の上の値の組にしたもの。
/// 不連続点はエネルギーの重複（左側の値、右側の値の順）で表し、表の定義域の外は 0 とする。
struct Union {
    energies: Vec<f64>,
    values: Vec<Vec<f64>>,
}

impl Union {
    fn new(tables: &[Tab1]) -> Self {
        let mut breakpoints: Vec<f64> = tables.iter().flat_map(|t| t.x.iter().copied()).collect();
        breakpoints.sort_by(f64::total_cmp);
        breakpoints.dedup();
        let last = breakpoints.len() - 1;
        let mut energies = Vec::with_capacity(breakpoints.len());
        let mut values = Vec::with_capacity(breakpoints.len());
        for (j, &e) in breakpoints.iter().enumerate() {
            let left: Vec<f64> = tables.iter().map(|t| left_limit(t, e)).collect();
            let right: Vec<f64> = tables.iter().map(|t| right_limit(t, e)).collect();
            let jump = left != right;
            if j > 0 {
                energies.push(e);
                values.push(left);
            }
            if j == 0 || (j < last && jump) {
                energies.push(e);
                values.push(right);
            }
        }
        Self { energies, values }
    }
}

/// 拡がりの積分。`x[k]` = √(α E_k) の各点で、表の値の組 `values[k]` を持つ。
struct Integral<'a> {
    alpha: f64,
    x: Vec<f64>,
    values: &'a [Vec<f64>],
}

impl Integral<'_> {
    /// エネルギー e での、拡げた断面積の組。
    fn at(&self, e: f64) -> Vec<f64> {
        let y = (self.alpha * e).sqrt();
        let x = &self.x;
        let n = self.values[0].len();
        let mut sum = vec![0.0; n];
        let mut add = |weight: f64, values: &[f64]| {
            for (s, v) in sum.iter_mut().zip(values) {
                *s += weight * v;
            }
        };
        let (lo, hi) = ((y - WINDOW).max(0.0), y + WINDOW);

        // 最初の点より下は 1/v: x² σ(x) = x x₀ σ₀
        let x0 = x[0];
        if lo < x0 {
            let mut weights = 0.0;
            quadrature(lo, hi.min(x0), |xi, w| {
                weights += w * xi * x0 * kernel(xi, y)
            });
            add(weights, &self.values[0]);
        }

        let start = x.partition_point(|&xk| xk <= lo).saturating_sub(1);
        for k in start..x.len() - 1 {
            let (x0, x1) = (x[k], x[k + 1]);
            if x0 >= hi {
                break;
            }
            let (a, b) = (x0.max(lo), x1.min(hi));
            if b <= a {
                continue;
            }
            // 区間の中で σ は E（= x²/α）について線形: σ = σ₀ + t (σ₁ - σ₀)、
            // t = (x - x₀)(x + x₀) / ((x₁ - x₀)(x₁ + x₀))。両端の重みに分けて積分する
            let denom = (x1 - x0) * (x1 + x0);
            let (mut w0, mut w1) = (0.0, 0.0);
            quadrature(a, b, |xi, w| {
                let f = w * xi * xi * kernel(xi, y);
                let t = (xi - x0) * (xi + x0) / denom;
                w0 += f * (1.0 - t);
                w1 += f * t;
            });
            add(w0, &self.values[k]);
            add(w1, &self.values[k + 1]);
        }

        // 最後の点より上は一定値
        let last = x.len() - 1;
        if hi > x[last] {
            let mut weights = 0.0;
            quadrature(lo.max(x[last]), hi, |xi, w| {
                weights += w * xi * xi * kernel(xi, y)
            });
            add(weights, &self.values[last]);
        }

        let scale = 1.0 / (y * y * SQRT_PI);
        sum.iter_mut().for_each(|s| *s *= scale);
        sum
    }
}

/// exp(-(x-y)²) - exp(-(x+y)²) を、y が小さいときにも桁落ちしない形で計算する。
fn kernel(x: f64, y: f64) -> f64 {
    let d = x - y;
    (-d * d).exp() * -(-4.0 * x * y).exp_m1()
}

/// 区間 [a, b] の Gauss-Legendre 求積。各点 x と重み w で `f(x, w)` を呼ぶ。
/// 被積分関数は x で幅 1 程度の構造（核関数）を持つので、区間の幅に応じて点数を選び、
/// 幅 1 を超える区間は幅 1 以下の小区間に分ける。
///
/// n 点則の誤差は w^(2n+1) (n!)⁴ / ((2n+1) ((2n)!)³) · f⁽²ⁿ⁾ で、核関数 exp(-z²) の 2n 階微分は
/// (2n)!/n! 倍以下なので、幅 w の区間の相対誤差は 2点則・w = 0.005 で 2e-12、4点則・w = 0.2 で 2e-12、
/// 8点則・w = 1 で 9e-15 になる（x² σ(x) の変化は核関数よりずっと緩やかなので無視できる）。
fn quadrature(a: f64, b: f64, mut f: impl FnMut(f64, f64)) {
    let width = b - a;
    let (rule, panels): (&[(f64, f64)], usize) = if width <= 0.005 {
        (&GAUSS_LEGENDRE_2, 1)
    } else if width <= 0.2 {
        (&GAUSS_LEGENDRE_4, 1)
    } else {
        (&GAUSS_LEGENDRE_8, width.ceil() as usize)
    };
    let h = width / panels as f64;
    for p in 0..panels {
        let center = a + h * (p as f64 + 0.5);
        let half = 0.5 * h;
        for &(node, weight) in rule {
            f(center + half * node, half * weight);
        }
    }
}

/// 拡げる上限エネルギー（NJOY BROADR の既定値と同じ規則）。
///
/// 分離共鳴があればその上限、非分離共鳴だけならその下限とする（どちらも `MAX_LIMIT` まで）。
/// 共鳴パラメータがなければ `MAX_LIMIT` とし、1 eV より上から始まるしきい反応があれば、
/// そのしきい値の手前までにする。
pub fn broadening_limit<'a>(mf2: Option<&Mf2>, tables: impl Iterator<Item = &'a Tab1>) -> f64 {
    let ranges = mf2.map_or(&[][..], |m| &m.ranges);
    let resolved = ranges
        .iter()
        .filter(|r| r.kind.is_resolved())
        .map(|r| r.eh)
        .reduce(f64::max);
    let unresolved = ranges
        .iter()
        .filter(|r| matches!(r.kind, RangeKind::Unresolved(_)))
        .map(|r| r.el)
        .reduce(f64::min);
    if let Some(e) = resolved.or(unresolved) {
        return e.min(MAX_LIMIT);
    }
    tables
        .map(|t| t.x[0])
        .filter(|&e| e > LOW_THRESHOLD)
        .fold(MAX_LIMIT, |limit, e| limit.min(BELOW_THRESHOLD * e))
}

#[cfg(test)]
mod tests;
