//! 非分離共鳴の幅のゆらぎ（平均 1 の χ² 分布）による期待値。
//!
//! 自由度 ν の χ² 分布の密度は k^k x^(k-1) e^(-kx) / Γ(k)（k = ν/2）。x = u² と置くと
//! 2 k^k u^(2k-1) e^(-k u²) / Γ(k) du になり、ν = 1 でも原点で有限になるので、u について適応積分する。

/// 期待値の相対許容誤差。群定数の許容誤差（1e-3）より十分小さくとる
const RELATIVE_TOLERANCE: f64 = 1.0e-10;

/// u の積分の上限を k u² = TAIL とする。密度は e^(-40) ≈ 4e-18 倍になり、倍精度で無視できる
const TAIL: f64 = 40.0;

/// 適応積分で区間を分ける深さの上限（区間の幅は 2^-50 倍まで）
const MAX_DEPTH: u32 = 50;

/// 自由度 `dof` の χ² 分布での f(x) の期待値。自由度 0 は幅がゆらがない（x = 1）ことを表す。
pub(super) fn expectation(dof: u32, f: impl Fn(f64) -> f64) -> f64 {
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
    scale * adaptive_simpson(integrand, 0.0, (TAIL / k).sqrt())
}

/// [a, b] の適応 Simpson 積分。全体の積分の `RELATIVE_TOLERANCE` 倍まで、区間を分ける。
fn adaptive_simpson(f: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    // 細い構造を最初の標本で見逃さないよう、はじめに区間を等分しておく
    const START: usize = 16;
    let h = (b - a) / START as f64;
    let mut panels: Vec<Panel> = (0..START)
        .map(|i| Panel::new(&f, a + h * i as f64, a + h * (i + 1) as f64, 0))
        .collect();
    let rough: f64 = panels.iter().map(|p| p.whole).sum();
    let mut eps = RELATIVE_TOLERANCE * rough.abs() / START as f64;
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

    #[test]
    fn reproduces_moments_of_chi_square() {
        // 平均 1 の χ²_ν の k 次の積率は Π_{i<k} (1 + 2i/ν)
        for dof in 1..=4 {
            let nu = f64::from(dof);
            let mut expected = 1.0;
            for k in 0..=6 {
                let moment = expectation(dof, |x| x.powi(k));
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
                let actual = expectation(dof, |x| 1.0 / (x + c));
                assert!(
                    (actual - exact).abs() <= 1e-8 * exact,
                    "ν={dof}, c={c}: {actual} vs {exact}"
                );
            }
        }
    }

    #[test]
    fn zero_degrees_of_freedom_do_not_fluctuate() {
        assert_eq!(expectation(0, |x| 3.0 * x + 1.0), 4.0);
    }
}
