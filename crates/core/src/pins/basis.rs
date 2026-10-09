//! ノードの中の ∇²ψ = μψ の解の 8 つの関数（ADR-0096）。
//!
//! x、y と2本の対角 u = (x + y)/√2、v = (x - y)/√2 の向きの平面波 C(ξ) = cosh(√μ ξ)、S(ξ) = sinh(√μ ξ)/√μ
//! （μ < 0 では cos と sin）の8つは、どれも ∇²ψ = μψ を満たす。μ → 0 ではどれも 1 か ξ に近づいて区別できなくなるので、
//! 同じ8つの関数が張る空間を、μ → 0 でも独立な組み合わせで表す:
//!
//! - 偶関数: C_x、(C_y - C_x)/μ、(C_u - C_v)/μ、(C_u + C_v - C_x - C_y)/μ²
//! - 奇関数: S_x、S_y、(S_u - (S_x + S_y)/√2)/μ、(S_v - (S_x - S_y)/√2)/μ
//!
//! μ → 0 の極限は 1、(y² - x²)/2、xy、4次式と x、y、3次式2つ。組み合わせの分子は μ のべきの低い項が打ち消し合うので、
//! |μ| (a² + b²) ≤ 1（a、b はノードの半幅）では、打ち消し合う項を除いた Taylor 級数で求める。それより大きい μ では
//! 閉じた式で求め、桁落ちは小さい。どちらも同じ関数なので、境で値は跳ばない。

use std::f64::consts::FRAC_1_SQRT_2;

/// 面の平均の Gauss-Legendre の点の数。集合体の幅で κ L ≈ 10 の指数関数を 1e-12 より細かく積分する
pub(crate) const FACE_POINTS: usize = 24;

/// C（偶）か S（奇）
#[derive(Debug, Clone, Copy, PartialEq)]
enum Kind {
    Even,
    Odd,
}

/// Σ coef F(cx x + cy y) / μ^power。F の Taylor 級数の μ^n の項は、n < power で打ち消し合う
#[derive(Debug, Clone)]
struct Combo {
    kind: Kind,
    /// (coef, cx, cy)
    terms: Vec<(f64, f64, f64)>,
    power: i32,
}

impl Combo {
    fn new(kind: Kind, terms: &[(f64, f64, f64)], power: i32) -> Self {
        Self {
            kind,
            terms: terms.to_vec(),
            power,
        }
    }

    /// x（axis 0）か y（axis 1）での微分。dC/dξ = μ S、dS/dξ = C
    fn derivative(&self, axis: usize, mu: f64) -> Self {
        let along = |&(c, cx, cy): &(f64, f64, f64)| (c * if axis == 0 { cx } else { cy }, cx, cy);
        let terms: Vec<_> = self.terms.iter().map(along).collect();
        match self.kind {
            Kind::Odd => Self::new(Kind::Even, &terms, self.power),
            // μ S / μ^p = S / μ^(p-1)。p = 0 なら μ を係数に掛ける
            Kind::Even if self.power > 0 => Self::new(Kind::Odd, &terms, self.power - 1),
            Kind::Even => {
                let scaled: Vec<_> = terms.iter().map(|&(c, cx, cy)| (mu * c, cx, cy)).collect();
                Self::new(Kind::Odd, &scaled, 0)
            }
        }
    }

    fn eval(&self, mu: f64, series: bool, x: f64, y: f64) -> f64 {
        if series {
            self.series(mu, x, y)
        } else {
            let root = mu.abs().sqrt();
            let f = |xi: f64| match (self.kind, mu > 0.0) {
                (Kind::Even, true) => (root * xi).cosh(),
                (Kind::Even, false) => (root * xi).cos(),
                (Kind::Odd, true) => (root * xi).sinh() / root,
                (Kind::Odd, false) => (root * xi).sin() / root,
            };
            let sum: f64 = self
                .terms
                .iter()
                .map(|&(c, cx, cy)| c * f(cx * x + cy * y))
                .sum();
            sum / mu.powi(self.power)
        }
    }

    /// Σ_{n ≥ power} μ^(n - power) Σ coef ξ^(2n + s) / (2n + s)!（s は偶で 0、奇で 1）
    fn series(&self, mu: f64, x: f64, y: f64) -> f64 {
        let s = match self.kind {
            Kind::Even => 0,
            Kind::Odd => 1,
        };
        let start = self.power.max(0);
        // 項は 4 つまで（`Basis::new`）。評価のたびに確保しないよう固定長の配列にする
        let count = self.terms.len().min(4);
        let mut xis = [(0.0, 0.0); 4];
        for (slot, &(c, cx, cy)) in xis.iter_mut().zip(&self.terms) {
            *slot = (c, cx * x + cy * y);
        }
        let xis = &xis[..count];
        // n = start の項の ξ^(2n+s)/(2n+s)!
        let first_order = 2 * start + s;
        let mut factorial = 1.0;
        for k in 1..=first_order {
            factorial *= f64::from(k);
        }
        let mut powers = [0.0; 4];
        for (p, &(_, xi)) in powers.iter_mut().zip(xis) {
            *p = xi.powi(first_order) / factorial;
        }
        let powers = &mut powers[..count];
        let mut sum = 0.0;
        let mut mu_power = 1.0;
        for n in start..start + 60 {
            let term: f64 = mu_power
                * xis
                    .iter()
                    .zip(powers.iter())
                    .map(|(&(c, _), p)| c * p)
                    .sum::<f64>();
            sum += term;
            if n > start && term.abs() <= 1e-17 * sum.abs() {
                break;
            }
            // 次の n: ξ² / ((2n+s+1)(2n+s+2)) を掛ける
            let order = f64::from(2 * n + s);
            for (p, &(_, xi)) in powers.iter_mut().zip(xis) {
                *p *= xi * xi / ((order + 1.0) * (order + 2.0));
            }
            mu_power *= mu;
        }
        sum
    }
}

/// ノード（半幅 a、b、中心が原点）の中の、固有値 μ の 8 つの関数と、一様な線源の特解。
#[derive(Debug, Clone)]
pub(crate) struct Basis {
    mu: f64,
    half: [f64; 2],
    series: bool,
    functions: Vec<Combo>,
    gradients: Vec<[Combo; 2]>,
    particular: Combo,
    particular_gradient: [Combo; 2],
}

impl Basis {
    pub(crate) fn new(mu: f64, a: f64, b: f64) -> Self {
        let r = FRAC_1_SQRT_2;
        let (x, y, u, v) = ((1.0, 0.0), (0.0, 1.0), (r, r), (r, -r));
        let t = |c: f64, (cx, cy): (f64, f64)| (c, cx, cy);
        let functions = vec![
            Combo::new(Kind::Even, &[t(1.0, x)], 0),
            Combo::new(Kind::Even, &[t(1.0, y), t(-1.0, x)], 1),
            Combo::new(Kind::Even, &[t(1.0, u), t(-1.0, v)], 1),
            Combo::new(
                Kind::Even,
                &[t(1.0, u), t(1.0, v), t(-1.0, x), t(-1.0, y)],
                2,
            ),
            Combo::new(Kind::Odd, &[t(1.0, x)], 0),
            Combo::new(Kind::Odd, &[t(1.0, y)], 0),
            Combo::new(Kind::Odd, &[t(1.0, u), t(-r, x), t(-r, y)], 1),
            Combo::new(Kind::Odd, &[t(1.0, v), t(-r, x), t(r, y)], 1),
        ];
        let gradients = functions
            .iter()
            .map(|f| [f.derivative(0, mu), f.derivative(1, mu)])
            .collect();
        // P = (C_x - 1)/μ は ∇²P = C_x = μP + 1
        let particular = Combo::new(Kind::Even, &[t(1.0, x), t(-1.0, (0.0, 0.0))], 1);
        let particular_gradient = [particular.derivative(0, mu), particular.derivative(1, mu)];
        Self {
            mu,
            half: [a, b],
            series: mu.abs() * (a * a + b * b) <= 1.0,
            functions,
            gradients,
            particular,
            particular_gradient,
        }
    }

    /// (x, y) での 8 つの関数の値
    pub(crate) fn values(&self, x: f64, y: f64) -> [f64; 8] {
        std::array::from_fn(|k| self.functions[k].eval(self.mu, self.series, x, y))
    }

    /// (x, y) での 8 つの関数の勾配 [∂/∂x, ∂/∂y]
    pub(crate) fn gradients(&self, x: f64, y: f64) -> [[f64; 2]; 8] {
        std::array::from_fn(|k| {
            let [gx, gy] = &self.gradients[k];
            [
                gx.eval(self.mu, self.series, x, y),
                gy.eval(self.mu, self.series, x, y),
            ]
        })
    }

    /// 面（-x、+x、-y、+y の順の番号）の上で平均した、外向きの法線微分。`rule` は [-1, 1] の Gauss-Legendre の点と重み
    pub(crate) fn face_derivatives(&self, face: usize, rule: &[(f64, f64)]) -> [f64; 8] {
        let (axis, sign) = (face / 2, if face.is_multiple_of(2) { -1.0 } else { 1.0 });
        let across = self.half[axis];
        let along = self.half[1 - axis];
        let mut out = [0.0; 8];
        for &(t, w) in rule {
            let (x, y) = if axis == 0 {
                (sign * across, t * along)
            } else {
                (t * along, sign * across)
            };
            for (o, g) in out.iter_mut().zip(self.gradients(x, y)) {
                // [-1, 1] の重みの和は 2 なので、平均は重み / 2
                *o += 0.5 * w * sign * g[axis];
            }
        }
        out
    }

    /// 長方形 [x0, x1] × [y0, y1] の上の 8 つの関数の平均（Gauss-Legendre の積）
    pub(crate) fn averages(
        &self,
        x0: f64,
        x1: f64,
        y0: f64,
        y1: f64,
        rule: &[(f64, f64)],
    ) -> [f64; 8] {
        let mut out = [0.0; 8];
        for (x, y, w) in points(x0, x1, y0, y1, rule) {
            for (o, v) in out.iter_mut().zip(self.values(x, y)) {
                *o += w * v;
            }
        }
        out
    }

    /// 特解 P の勾配 [∂/∂x, ∂/∂y]
    pub(crate) fn particular_gradient(&self, x: f64, y: f64) -> [f64; 2] {
        self.particular_gradient
            .each_ref()
            .map(|g| g.eval(self.mu, self.series, x, y))
    }

    /// 面（-x、+x、-y、+y）の上で平均した、特解 P の外向きの法線微分
    pub(crate) fn particular_face_derivative(&self, face: usize, rule: &[(f64, f64)]) -> f64 {
        let (axis, sign) = (face / 2, if face.is_multiple_of(2) { -1.0 } else { 1.0 });
        let (across, along) = (self.half[axis], self.half[1 - axis]);
        rule.iter()
            .map(|&(t, w)| {
                let (x, y) = if axis == 0 {
                    (sign * across, t * along)
                } else {
                    (t * along, sign * across)
                };
                0.5 * w * sign * self.particular_gradient(x, y)[axis]
            })
            .sum()
    }

    /// 長方形の上の特解 P の平均
    pub(crate) fn particular_average(
        &self,
        x0: f64,
        x1: f64,
        y0: f64,
        y1: f64,
        rule: &[(f64, f64)],
    ) -> f64 {
        points(x0, x1, y0, y1, rule)
            .map(|(x, y, w)| w * self.particular(x, y))
            .sum()
    }

    /// ∇²P - μP = 1 を満たす特解 P（ノードの中で一様な線源）
    pub(crate) fn particular(&self, x: f64, y: f64) -> f64 {
        self.particular.eval(self.mu, self.series, x, y)
    }
}

/// 長方形の上の Gauss-Legendre の積の点 (x, y) と、平均の重み（和が 1）
fn points(
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    rule: &[(f64, f64)],
) -> impl Iterator<Item = (f64, f64, f64)> + '_ {
    let (cx, hx, cy, hy) = (
        0.5 * (x0 + x1),
        0.5 * (x1 - x0),
        0.5 * (y0 + y1),
        0.5 * (y1 - y0),
    );
    rule.iter().flat_map(move |&(s, ws)| {
        rule.iter()
            .map(move |&(t, wt)| (cx + hx * s, cy + hy * t, 0.25 * ws * wt))
    })
}

#[cfg(test)]
mod tests;
