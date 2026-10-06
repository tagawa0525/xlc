//! 非分離共鳴の幅のゆらぎ（平均 1 の χ² 分布）による期待値。
//!
//! 自由度 ν の χ² 分布の密度は k^k x^(k-1) e^(-kx) / Γ(k)（k = ν/2）。x = u² と置くと
//! 2 k^k u^(2k-1) e^(-k u²) / Γ(k) du になり、ν = 1 でも原点で有限になるので、u について適応積分する。

use crate::special::{erfcx, exp_e1};

/// u の積分の上限を k u² = TAIL とする。密度は e^(-40) ≈ 4e-18 倍になり、倍精度で無視できる
const TAIL: f64 = 40.0;

/// 適応積分で区間を分ける深さの上限（区間の幅は 2^-50 倍まで）
const MAX_DEPTH: u32 = 50;

/// これより大きい c では、閉じた式の代わりに 1/c の級数を使う。閉じた式は
/// E[x²/(x+c)] = 1 - c E[x/(x+c)] の引き算で相対誤差が eps c² 程度になり（c = 100 で 2e-12）、
/// 級数の項 m_n/c^n は c ≥ 100 なら 20 項ほどで倍精度の大きさを下回るため
const SERIES_ABOVE: f64 = 100.0;

/// 自由度 `dof` の χ² 分布での [E[x/(x + c)], E[x²/(x + c)]]（c > 0）。
///
/// x/(x+c) = 1 - c/(x+c) と x²/(x+c) = x - c x/(x+c) により E[1/(x+c)] に帰着させる。
/// E[1/(x+c)] = k^k c^(k-1) e^(kc) Γ(1-k, kc)（k = ν/2）は、ν = 1, 3 では erfc、ν = 2, 4 では E1 で書ける。
pub(super) fn width_ratios(dof: u32, c: f64) -> [f64; 2] {
    if dof == 0 {
        let r = 1.0 / (1.0 + c);
        return [r, r];
    }
    if c > SERIES_ABOVE {
        return width_ratio_series(dof, c);
    }
    let k = f64::from(dof) / 2.0;
    let sqrt_pi = std::f64::consts::PI.sqrt();
    let inverse = match dof {
        1 => (std::f64::consts::PI * k / c).sqrt() * erfcx((k * c).sqrt()),
        2 => exp_e1(c),
        3 => 2.0 * k - 2.0 * sqrt_pi * k.powf(1.5) * c.sqrt() * erfcx((k * c).sqrt()),
        4 => k - k * k * c * exp_e1(k * c),
        _ => panic!("no χ² density for {dof} degrees of freedom"),
    };
    let r = 1.0 - c * inverse;
    [r, 1.0 - c * r]
}

/// 大きな c での [E[x/(x+c)], E[x²/(x+c)]] の漸近級数。t_n = m_n/c^n（m_n = Π_{i<n} (1 + 2i/ν) は
/// 平均 1 の χ²_ν の積率）として、E[x/(x+c)] = Σ_{n≥1} (-1)^(n+1) t_n、E[x²/(x+c)] = Σ_{n≥1} (-1)^(n+1) c t_(n+1)。
fn width_ratio_series(dof: u32, c: f64) -> [f64; 2] {
    const MAX_TERMS: u32 = 60;
    let nu = f64::from(dof);
    let (mut first, mut second) = (0.0, 0.0);
    let mut t = 1.0 / c;
    let mut sign = 1.0;
    for n in 1..=MAX_TERMS {
        let next = t * (1.0 + 2.0 * f64::from(n) / nu) / c;
        first += sign * t;
        second += sign * c * next;
        if t <= f64::EPSILON * first.abs() && c * next <= f64::EPSILON * second.abs() {
            return [first, second];
        }
        t = next;
        sign = -sign;
    }
    panic!("width ratio series did not converge for c = {c}");
}

/// 自由度 `dof` の χ² 分布での f(x) の期待値（相対許容誤差 `tolerance`）。
/// 自由度 0 は幅がゆらがない（x = 1）ことを表す。
///
/// 期待値を入れ子にするときは、外側の許容誤差を内側より十分大きくする。内側の結果は、分割の仕方が
/// 切り替わるところで許容誤差の大きさだけ飛ぶので、外側が同じ精度を求めると分割が止まらない。
pub(super) fn expectation(dof: u32, tolerance: f64, f: impl Fn(f64) -> f64) -> f64 {
    let k = f64::from(dof) / 2.0;
    // Γ(k)（k = 1/2, 1, 3/2, 2）
    let gamma = match dof {
        0 => return f(1.0),
        1 => std::f64::consts::PI.sqrt(),
        2 | 4 => 1.0,
        3 => 0.5 * std::f64::consts::PI.sqrt(),
        _ => panic!("no χ² density for {dof} degrees of freedom"),
    };
    let scale = 2.0 * k.powf(k) / gamma;
    let integrand = |u: f64| {
        let x = u * u;
        f(x) * u.powf(2.0 * k - 1.0) * (-k * x).exp()
    };
    scale * adaptive_simpson(integrand, 0.0, (TAIL / k).sqrt(), tolerance)
}

/// [a, b] の適応 Simpson 積分。全体の積分の `tolerance` 倍まで、区間を分ける。
fn adaptive_simpson(f: impl Fn(f64) -> f64, a: f64, b: f64, tolerance: f64) -> f64 {
    // 細い構造を最初の標本で見逃さないよう、はじめに区間を等分しておく
    const START: usize = 16;
    let h = (b - a) / START as f64;
    let mut panels: Vec<Panel> = (0..START)
        .map(|i| Panel::new(&f, a + h * i as f64, a + h * (i + 1) as f64, 0))
        .collect();
    let rough: f64 = panels.iter().map(|p| p.whole).sum();
    let mut eps = tolerance * rough.abs() / START as f64;
    if eps == 0.0 {
        eps = f64::MIN_POSITIVE;
    }
    let mut sum = 0.0;
    while let Some(p) = panels.pop() {
        let (left, right) = p.split(&f);
        let delta = left.whole + right.whole - p.whole;
        // 誤差の目安は |delta|/15。許容誤差は区間の幅に比例して配る
        let allowed = eps * (p.b - p.a) / h;
        if delta.abs() <= 15.0 * allowed {
            sum += left.whole + right.whole + delta / 15.0;
        } else {
            assert!(
                p.depth < MAX_DEPTH,
                "adaptive integration did not converge on [{}, {}]",
                p.a,
                p.b
            );
            panels.push(left);
            panels.push(right);
        }
    }
    sum
}

/// Simpson 則の1区間（両端と中点の値を持つ）。
struct Panel {
    a: f64,
    b: f64,
    fa: f64,
    fm: f64,
    fb: f64,
    whole: f64,
    depth: u32,
}

impl Panel {
    fn new(f: &impl Fn(f64) -> f64, a: f64, b: f64, depth: u32) -> Self {
        Self::with_ends(f, a, b, f(a), f(b), depth)
    }

    fn with_ends(f: &impl Fn(f64) -> f64, a: f64, b: f64, fa: f64, fb: f64, depth: u32) -> Self {
        let fm = f(0.5 * (a + b));
        Self {
            a,
            b,
            fa,
            fm,
            fb,
            whole: (b - a) / 6.0 * (fa + 4.0 * fm + fb),
            depth,
        }
    }

    fn split(&self, f: &impl Fn(f64) -> f64) -> (Self, Self) {
        let m = 0.5 * (self.a + self.b);
        let left = Self::with_ends(f, self.a, m, self.fa, self.fm, self.depth + 1);
        let right = Self::with_ends(f, m, self.b, self.fm, self.fb, self.depth + 1);
        (left, right)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::special::erfc;

    const SQRT_PI: f64 = 1.772_453_850_905_516;
    const INNER: f64 = 1.0e-10;
    const OUTER: f64 = 1.0e-8;

    #[test]
    fn reproduces_moments_of_chi_square() {
        // 平均 1 の χ²_ν の k 次の積率は Π_{i<k} (1 + 2i/ν)
        for dof in 1..=4 {
            let nu = f64::from(dof);
            let mut expected = 1.0;
            for k in 0..=6 {
                let moment = expectation(dof, INNER, |x| x.powi(k));
                assert!(
                    (moment - expected).abs() <= 1e-9 * expected,
                    "ν={dof}, k={k}: {moment} vs {expected}"
                );
                expected *= 1.0 + 2.0 * f64::from(k) / nu;
            }
        }
    }

    #[test]
    fn matches_closed_form_near_the_pole() {
        // E[1/(x+c)] = k^k c^(k-1) e^(kc) Γ(1-k, kc)。
        // Γ(1/2, z) = √π erfc(√z)、Γ(-1/2, z) = 2 e^(-z)/√z - 2√π erfc(√z)
        for c in [1.0e-4_f64, 1.0e-2, 1.0, 10.0] {
            let z = 0.5 * c;
            let nu1 = 0.5_f64.sqrt() / c.sqrt() * z.exp() * SQRT_PI * erfc(z.sqrt());
            let z = 1.5 * c;
            let gamma = 2.0 * (-z).exp() / z.sqrt() - 2.0 * SQRT_PI * erfc(z.sqrt());
            let nu3 = 1.5_f64.powf(1.5) * c.sqrt() * z.exp() * gamma;
            for (dof, exact) in [(1, nu1), (3, nu3)] {
                let actual = expectation(dof, INNER, |x| 1.0 / (x + c));
                assert!(
                    (actual - exact).abs() <= 1e-8 * exact,
                    "ν={dof}, c={c}: {actual} vs {exact}"
                );
            }
        }
    }

    #[test]
    fn nested_expectations_converge_in_either_order() {
        // 内側の期待値は許容誤差の範囲でばらつくので、外側は許容誤差を緩めないと収束しない
        let f = |x: f64, y: f64| x * y / (x + 0.3 * y + 0.01);
        let x_inside = expectation(1, OUTER, |y| expectation(1, INNER, |x| f(x, y)));
        let y_inside = expectation(1, OUTER, |x| expectation(1, INNER, |y| f(x, y)));
        assert!(
            (x_inside - y_inside).abs() <= 1e-7 * x_inside,
            "{x_inside} vs {y_inside}"
        );
    }

    #[test]
    fn width_ratios_match_integration() {
        // c が大きいところ（Γn ≪ Γγ、たとえば低エネルギーの l=2）まで確かめる
        for dof in 0..=4 {
            for c in [
                1.0e-6, 1.0e-3, 0.1, 1.0, 10.0, 99.0, 101.0, 1.0e3, 1.0e6, 1.0e9,
            ] {
                let numeric = [
                    expectation(dof, 1.0e-12, |x| x / (x + c)),
                    expectation(dof, 1.0e-12, |x| x * x / (x + c)),
                ];
                let closed = width_ratios(dof, c);
                for i in 0..2 {
                    assert!(
                        (closed[i] - numeric[i]).abs() <= 1e-9 * numeric[i],
                        "ν={dof}, c={c}, [{i}]: {} vs {}",
                        closed[i],
                        numeric[i]
                    );
                }
            }
        }
    }

    #[test]
    fn zero_degrees_of_freedom_do_not_fluctuate() {
        assert_eq!(expectation(0, INNER, |x| 3.0 * x + 1.0), 4.0);
    }
}
