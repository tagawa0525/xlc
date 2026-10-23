//! 特殊関数（標準ライブラリにない誤差関数など）。

use num_complex::Complex64;

use crate::numeric::NumericError;

/// 誤差関数 erf(x)。
pub fn erf(x: f64) -> f64 {
    if x < 0.0 {
        return -erf(-x);
    }
    if x < SERIES_LIMIT {
        erf_series(x)
    } else {
        1.0 - erfc(x)
    }
}

/// これより小さい x では級数、大きい x では連分数を使う
const SERIES_LIMIT: f64 = 2.5;

/// erf(x) = 2/√π exp(-x²) Σ 2ⁿ x^(2n+1) / (1·3·…·(2n+1))。項がすべて正なので桁落ちしない。
fn erf_series(x: f64) -> f64 {
    let x2 = x * x;
    let mut term = x;
    let mut sum = x;
    let mut n = 0.0;
    while term > sum * f64::EPSILON {
        n += 1.0;
        term *= 2.0 * x2 / (2.0 * n + 1.0);
        sum += term;
    }
    2.0 / std::f64::consts::PI.sqrt() * (-x2).exp() * sum
}

/// 相補誤差関数 erfc(x) = 1 - erf(x)。
pub fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        2.0 - erfc(-x)
    } else if x < SERIES_LIMIT {
        1.0 - erf_series(x)
    } else {
        erfcx(x) * (-x * x).exp()
    }
}

/// スケールした相補誤差関数 erfcx(x) = exp(x²) erfc(x)（x ≥ 0）。
/// exp(x²) と erfc(x) を別に計算すると、大きな x で桁あふれ・桁落ちする。
pub fn erfcx(x: f64) -> f64 {
    if x < SERIES_LIMIT {
        return (x * x).exp() * (1.0 - erf_series(x));
    }
    // 連分数 erfcx(x) = 1/√π · 1/(x + (1/2)/(x + 1/(x + (3/2)/(x + …)))) を後ろから評価する。
    // x ≥ 2.5 では 40 項で倍精度の精度になる
    const TERMS: u32 = 40;
    let mut t = x;
    for k in (1..=TERMS).rev() {
        t = x + 0.5 * f64::from(k) / t;
    }
    1.0 / (std::f64::consts::PI.sqrt() * t)
}

/// スケールした指数積分 e^z E1(z)（z > 0）。E1(z) = ∫_z^∞ e^(-t)/t dt。
///
/// z ≤ 1 では級数 E1 = -γ - ln z - Σ (-z)^n/(n n!)、z > 1 では連分数
/// E1 = e^(-z) / (z + 1 - 1²/(z + 3 - 2²/(z + 5 - …)))（修正 Lentz 法）で計算する。
pub fn exp_e1(z: f64) -> Result<f64, NumericError> {
    if !(z > 0.0 && z.is_finite()) {
        return Err(NumericError::Domain(format!("E1({z})")));
    }
    const MAX_TERMS: u32 = 1000;
    if z <= 1.0 {
        let (mut sum, mut term) = (0.0, 1.0);
        for n in 1..=MAX_TERMS {
            let n = f64::from(n);
            term *= -z / n;
            sum += term / n;
            if (term / n).abs() <= f64::EPSILON * sum.abs() {
                return Ok((-std::f64::consts::EULER_GAMMA - z.ln() - sum) * z.exp());
            }
        }
    } else {
        let tiny = f64::MIN_POSITIVE;
        let mut b = z + 1.0;
        let (mut c, mut d) = (1.0 / tiny, 1.0 / b);
        let mut h = d;
        for i in 1..=MAX_TERMS {
            let a = -f64::from(i * i);
            b += 2.0;
            d = 1.0 / (a * d + b);
            c = b + a / c;
            let delta = c * d;
            h *= delta;
            if (delta - 1.0).abs() <= f64::EPSILON {
                return Ok(h);
            }
        }
    }
    Err(NumericError::NotConverged(format!("E1({z})")))
}

/// 正則化した下側不完全ガンマ関数 P(3/2, x) = γ(3/2, x)/Γ(3/2)（x ≥ 0）。
/// Maxwell 分布 √E e^(-E/θ) のうち、E/θ < x の割合。
pub fn gamma_p_3_2(x: f64) -> f64 {
    if x < GAMMA_SERIES_LIMIT {
        // Γ(5/2) = 3√π/4
        x * x.sqrt() * (-x).exp() * lower_gamma_series(1.5, x)
            / (0.75 * std::f64::consts::PI.sqrt())
    } else {
        erf(x.sqrt()) - 2.0 * (x / std::f64::consts::PI).sqrt() * (-x).exp()
    }
}

/// 正則化した下側不完全ガンマ関数 P(2, x) = 1 - e^(-x)(1 + x)（x ≥ 0）。
/// 蒸発スペクトル E e^(-E/θ) のうち、E/θ < x の割合。
pub fn gamma_p_2(x: f64) -> f64 {
    if x < GAMMA_SERIES_LIMIT {
        // Γ(3) = 2
        0.5 * x * x * (-x).exp() * lower_gamma_series(2.0, x)
    } else {
        1.0 - (-x).exp() * (1.0 + x)
    }
}

/// これより小さい x では級数、大きい x では閉じた式を使う。閉じた式は小さい x で 1 に近い
/// 2つの項の差になって桁落ちし、級数は大きい x で収束が遅くなる
const GAMMA_SERIES_LIMIT: f64 = 1.0;

/// P(a, x) = x^a e^(-x)/Γ(a+1) Σ xⁿ/((a+1)…(a+n)) の和の部分。項がすべて正なので桁落ちしない。
fn lower_gamma_series(a: f64, x: f64) -> f64 {
    let (mut term, mut sum, mut n) = (1.0, 1.0, 0.0);
    while term > sum * f64::EPSILON {
        n += 1.0;
        term *= x / (a + n);
        sum += term;
    }
    sum
}

/// Faddeeva 関数 w(z) = e^(-z²) erfc(-iz)（z = x + iy）。
///
/// y ≥ 0 では Weideman (1994, SIAM J. Numer. Anal. 31) の有理近似 w(z) ≈ 2 p(Z)/(L - iz)² + 1/(√π (L - iz))、
/// Z = (L + iz)/(L - iz) を使う。p は `FADDEEVA_TERMS` 次の多項式で、係数は `faddeeva_coefficients`。
/// y < 0 では w(z) = 2 e^(-z²) - w(-z) で上半面に移す。
pub fn faddeeva(z: Complex64) -> Complex64 {
    let (x, y) = (z.re, z.im);
    if y < 0.0 {
        let w = faddeeva(-z);
        // e^(-z²) = e^(y² - x²) (cos 2xy - i sin 2xy)
        let scale = 2.0 * (y * y - x * x).exp();
        let phase = 2.0 * x * y;
        return Complex64::new(scale * phase.cos() - w.re, -scale * phase.sin() - w.im);
    }
    let a = faddeeva_coefficients();
    let l = faddeeva_scale();
    // L - iz = (L + y) - ix、L + iz = (L - y) + ix
    let den = Complex64::new(l + y, -x);
    let ratio = Complex64::new(l - y, x) / den;
    // Horner 法で p(Z) = Σ a_m Z^(m-1)
    let mut p = Complex64::ZERO;
    for &coefficient in a.iter().rev() {
        p = p * ratio + coefficient;
    }
    let first = p * 2.0 / (den * den);
    let second = Complex64::new(1.0 / std::f64::consts::PI.sqrt(), 0.0) / den;
    first + second
}

/// `faddeeva` を、同じ y ≥ 0 と並んだ x について求め、`w` に書く。
///
/// 点ごとの演算は `faddeeva` と同じ式の同じ順序で、結果のビットも同じ。Horner 法の項ごとに前の項に依るので、点を内側の
/// ループにして、独立な点を並べて（AVX-512 などのベクトル命令で）計算できるようにする（ADR-0188）。
pub fn faddeeva_upper(y: f64, x: &[f64], w: &mut [Complex64]) {
    assert!(y >= 0.0, "faddeeva_upper needs y >= 0, got {y}");
    assert_eq!(x.len(), w.len());
    let a = faddeeva_coefficients();
    let l = faddeeva_scale();
    let mut den = [Complex64::ZERO; FADDEEVA_CHUNK];
    let mut ratio = [Complex64::ZERO; FADDEEVA_CHUNK];
    let mut p = [Complex64::ZERO; FADDEEVA_CHUNK];
    for (x, w) in x.chunks(FADDEEVA_CHUNK).zip(w.chunks_mut(FADDEEVA_CHUNK)) {
        let n = x.len();
        let (den, ratio, p) = (&mut den[..n], &mut ratio[..n], &mut p[..n]);
        for (((&x, den), ratio), p) in x
            .iter()
            .zip(den.iter_mut())
            .zip(ratio.iter_mut())
            .zip(p.iter_mut())
        {
            *den = Complex64::new(l + y, -x);
            *ratio = Complex64::new(l - y, x) / *den;
            *p = Complex64::ZERO;
        }
        for &coefficient in a.iter().rev() {
            for (p, &ratio) in p.iter_mut().zip(ratio.iter()) {
                *p = *p * ratio + coefficient;
            }
        }
        for ((w, &p), &den) in w.iter_mut().zip(p.iter()).zip(den.iter()) {
            let first = p * 2.0 / (den * den);
            let second = Complex64::new(1.0 / std::f64::consts::PI.sqrt(), 0.0) / den;
            *w = first + second;
        }
    }
}

/// `faddeeva_upper` が一度に並べる点の数（スタックの作業配列の大きさ）
const FADDEEVA_CHUNK: usize = 64;

/// Weideman の近似の項数 N。N = 32 で、定義の数値積分と相対 1e-10、w(iy) = erfcx(y) と 1e-12、
/// 実軸上の e^(-x²) と絶対 1e-13 で一致する（テストを参照）
const FADDEEVA_TERMS: usize = 32;

/// Weideman の近似の尺度 L = (N/√2)^(1/2)
fn faddeeva_scale() -> f64 {
    (FADDEEVA_TERMS as f64 / std::f64::consts::SQRT_2).sqrt()
}

/// p の係数 a_1..a_N。f(t) = e^(-t²)(L² + t²) を t_k = L tan(kπ/(2M))（M = 2N, |k| < M）で標本化し、
/// a_m = (1/2M) Σ_k f(t_k) cos(π m k / M)（Weideman の FFT を偶関数として直接計算したもの）。
fn faddeeva_coefficients() -> &'static [f64; FADDEEVA_TERMS] {
    static COEFFICIENTS: std::sync::OnceLock<[f64; FADDEEVA_TERMS]> = std::sync::OnceLock::new();
    COEFFICIENTS.get_or_init(|| {
        let m = 2 * FADDEEVA_TERMS;
        let l = faddeeva_scale();
        let samples: Vec<(f64, f64)> = (1 - m as i64..m as i64)
            .map(|k| {
                let theta = k as f64 * std::f64::consts::PI / m as f64;
                let t = l * (0.5 * theta).tan();
                (k as f64, (-t * t).exp() * (l * l + t * t))
            })
            .collect();
        std::array::from_fn(|i| {
            let order = (i + 1) as f64;
            let sum: f64 = samples
                .iter()
                .map(|&(k, f)| f * (std::f64::consts::PI * order * k / m as f64).cos())
                .sum();
            sum / (2 * m) as f64
        })
    })
}

#[cfg(test)]
mod tests;
