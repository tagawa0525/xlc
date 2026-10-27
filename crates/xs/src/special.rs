//! 特殊関数（標準ライブラリにない誤差関数など）。

use num_complex::Complex64;

use crate::numeric::NumericError;

/// 誤差関数 erf(x)。
pub fn erf(x: f64) -> f64 {
    if x < 0.0 {
        return -erf(-x);
    }
    if x < CODY_LOWER_LIMIT {
        erf_series(x)
    } else {
        // erfc の exp_square でなく exp(−x*x) で、lattice の `erf` と同じビットにする。x*x の丸めの効きは erfc(x) 倍に縮む
        1.0 - erfcx(x) * (-x * x).exp()
    }
}

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
    } else if x < CODY_LOWER_LIMIT {
        1.0 - erf_series(x)
    } else {
        erfcx(x) * exp_square(x, -1.0)
    }
}

/// スケールした相補誤差関数 erfcx(x) = exp(x²) erfc(x)。
/// exp(x²) と erfc(x) を別に計算すると、大きな x で桁あふれ・桁落ちする。
pub fn erfcx(x: f64) -> f64 {
    if x < 0.0 {
        // erfc(−x) = 2 − erfc(x)
        return 2.0 * exp_square(x, 1.0) - erfcx(-x);
    }
    if x < CODY_LOWER_LIMIT {
        // erf(x) < 0.5 なので 1 − erf は桁落ちしない
        return (x * x).exp() * (1.0 - erf_series(x));
    }
    if x <= CODY_MIDDLE_LIMIT {
        compensated_rational(&CODY_MIDDLE_NUMERATOR, &CODY_MIDDLE_DENOMINATOR, x)
    } else {
        // erfcx(x) = (1/√π − w R(w))/x、w = 1/x²。x² が桁あふれする x では w = 0 で 1/(√π x) になる
        let w = 1.0 / (x * x);
        let tail = w * horner(&CODY_TAIL_NUMERATOR, w) / horner(&CODY_TAIL_DENOMINATOR, w);
        (0.5 * std::f64::consts::FRAC_2_SQRT_PI - tail) / x
    }
}

/// exp(sign·x²) を、x² の丸めを指数に入れずに求める。x を 1/16 の倍数 s に切り捨てると s² は丸めなしで表せ（|x| < 2²² で。
/// それより大きい x の exp(±x²) は桁あふれか 0）、残りの x² − s² = (x − s)(x + s) は 1 より小さいので、その丸めは指数で増えない
/// （CALERF と同じ分け方）。exp(±x*x) と書くと、x*x の丸めが x² 倍で効く（x = 24.6 で 510 ulp。ADR-0211）
fn exp_square(x: f64, sign: f64) -> f64 {
    let s = (x * 16.0).trunc() / 16.0;
    let rest = (x - s) * (x + s);
    (sign * s * s).exp() * (sign * rest).exp()
}

// Cody の有理近似の係数（ADR-0211）。W. J. Cody, "Rational Chebyshev approximations for the error function",
// Math. Comp. 23 (1969) 631–637 の近似を、Netlib の specfun の CALERF（倍精度の DATA 文。1990-03-19 版）からそのままコピーし、
// 高次の係数から並べた。コピーした値は、多倍長（mpmath）で有理式を評価し、erfcx との相対差が 0.46875 ≤ x ≤ 4 で 1.2e−18、
// 4 ≤ x ≤ 4e4 で 5.9e−19 以下になることで確かめた。出典と照らせるように、f64 で表せる桁より長い値も丸めずに残す

/// これより小さい x では erf を級数で、erfcx をその erf から求め、これ以上の x では erfcx を Cody の有理近似で、erf を
/// 1 − erfcx(x) exp(−x²) で求める（CALERF の区間の境 THRESH）。級数の 1 − erf は erf が 1 に近いほど桁落ちし、級数の erf は
/// 掛ける exp(−x²) に x² の丸めが x² 倍で効く（ADR-0211）
const CODY_LOWER_LIMIT: f64 = 0.46875;

/// これ以下の x では erfcx を 8 次 / 8 次の有理式 P(x)/Q(x) で、より大きい x では 1/x² の 5 次 / 5 次の有理式で求める
/// （CALERF の区間の境）
const CODY_MIDDLE_LIMIT: f64 = 4.0;

/// P(x) の係数（CALERF の C(9), C(1), …, C(8)）
#[allow(clippy::excessive_precision)]
const CODY_MIDDLE_NUMERATOR: [f64; 9] = [
    2.15311535474403846e-8,
    5.64188496988670089e-1,
    8.88314979438837594e0,
    6.61191906371416295e1,
    2.98635138197400131e2,
    8.81952221241769090e2,
    1.71204761263407058e3,
    2.05107837782607147e3,
    1.23033935479799725e3,
];

/// Q(x) の係数（1, CALERF の D(1), …, D(8)）
#[allow(clippy::excessive_precision)]
const CODY_MIDDLE_DENOMINATOR: [f64; 9] = [
    1.0,
    1.57449261107098347e1,
    1.17693950891312499e2,
    5.37181101862009858e2,
    1.62138957456669019e3,
    3.29079923573345963e3,
    4.36261909014324716e3,
    3.43936767414372164e3,
    1.23033935480374942e3,
];

/// x > 4 の R(w) = P(w)/Q(w) の分子の係数（CALERF の P(6), P(1), …, P(5)）
#[allow(clippy::excessive_precision)]
const CODY_TAIL_NUMERATOR: [f64; 6] = [
    1.63153871373020978e-2,
    3.05326634961232344e-1,
    3.60344899949804439e-1,
    1.25781726111229246e-1,
    1.60837851487422766e-2,
    6.58749161529837803e-4,
];

/// R(w) の分母の係数（1, CALERF の Q(1), …, Q(5)）
#[allow(clippy::excessive_precision)]
const CODY_TAIL_DENOMINATOR: [f64; 6] = [
    1.0,
    2.56852019228982242e0,
    1.87295284992346047e0,
    5.27905102951428412e-1,
    6.05183413124413191e-2,
    2.33520497626869185e-3,
];

/// 高次から並べた係数の多項式の Horner 法。積和は融合積和（mul_add）で、FMA の無い CPU でも同じビットになる
fn horner(coefficients: &[f64], x: f64) -> f64 {
    coefficients[1..]
        .iter()
        .fold(coefficients[0], |sum, &c| sum.mul_add(x, c))
}

/// 高次から並べた係数の 2 つの多項式の比 P(x)/Q(x) を、補償 Horner 法（Graillat, Langlois, Louvet, "Compensated Horner
/// scheme", 2005）で求める。各段の積と和の丸めの誤差を誤差なしの変換（積は融合積和、和は Knuth の TwoSum）で取り出し、別の
/// Horner 法で足し戻す。係数と x が正の多項式では、P と Q はそれぞれ倍精度の丸め 1 回ほどの誤差になる。普通の Horner 法では
/// 丸めが積もって erfcx の誤差が 5.5 ulp に達する（ADR-0211）。P と Q を 1 つのループで進めるのは、独立な 2 つの依存の列を
/// 並べて待ちを隠すため（別々のループの 1.7 倍の速さ）
fn compensated_rational<const N: usize>(
    numerator: &[f64; N],
    denominator: &[f64; N],
    x: f64,
) -> f64 {
    let step = |sum: f64, error: f64, c: f64| {
        let product = sum * x;
        let product_error = sum.mul_add(x, -product);
        let next = product + c;
        let z = next - product;
        let sum_error = (product - (next - z)) + (c - z);
        (next, error.mul_add(x, product_error + sum_error))
    };
    let (mut p, mut p_error) = (numerator[0], 0.0);
    let (mut q, mut q_error) = (denominator[0], 0.0);
    for (&a, &b) in numerator[1..].iter().zip(&denominator[1..]) {
        (p, p_error) = step(p, p_error, a);
        (q, q_error) = step(q, q_error, b);
    }
    (p + p_error) / (q + q_error)
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

/// `faddeeva` を、同じ y と並んだ x について求め、`w` の後ろに足す。
///
/// 点ごとの演算は `faddeeva` と同じ式の同じ順序で、結果のビットも同じ。y ≥ 0 では、Horner 法の項ごとに前の項に依るので、点を
/// 内側のループにして、独立な点を並べて（AVX-512 などのベクトル命令で）計算できるようにする（ADR-0188）。y < 0 は点ごとに
/// `faddeeva` で求める。
pub fn faddeeva_along_x(y: f64, x: &[f64], w: &mut Vec<Complex64>) {
    if y < 0.0 {
        w.extend(x.iter().map(|&x| faddeeva(Complex64::new(x, y))));
        return;
    }
    w.reserve(x.len());
    let a = faddeeva_coefficients();
    let l = faddeeva_scale();
    let inverse_root_pi = 1.0 / std::f64::consts::PI.sqrt();
    // 複素数の実部と虚部を別の配列に置き（num_complex の積と商と同じ式を成分で書く）、点ごとに独立な Horner 法のチェーンを
    // FADDEEVA_CHUNK 本並べて、前の項を待つ間に別の点を進める
    for x in x.chunks(FADDEEVA_CHUNK) {
        let n = x.len();
        let mut xs = [0.0; FADDEEVA_CHUNK];
        xs[..n].copy_from_slice(x);
        // L − iz = (L + y) − ix、L + iz = (L − y) + ix
        let (den_re, den_im) = (l + y, xs.map(|x| -x));
        let mut ratio_re = [0.0; FADDEEVA_CHUNK];
        let mut ratio_im = [0.0; FADDEEVA_CHUNK];
        for i in 0..FADDEEVA_CHUNK {
            let (a_re, a_im) = (l - y, xs[i]);
            let norm = den_re * den_re + den_im[i] * den_im[i];
            ratio_re[i] = (a_re * den_re + a_im * den_im[i]) / norm;
            ratio_im[i] = (a_im * den_re - a_re * den_im[i]) / norm;
        }
        // Horner 法で p(Z) = Σ a_m Z^(m-1)
        let mut p_re = [0.0; FADDEEVA_CHUNK];
        let mut p_im = [0.0; FADDEEVA_CHUNK];
        for &coefficient in a.iter().rev() {
            for i in 0..FADDEEVA_CHUNK {
                let re = p_re[i] * ratio_re[i] - p_im[i] * ratio_im[i];
                let im = p_re[i] * ratio_im[i] + p_im[i] * ratio_re[i];
                p_re[i] = re + coefficient;
                p_im[i] = im;
            }
        }
        for i in 0..n {
            // first = p * 2 / (den den)
            let (twice_re, twice_im) = (p_re[i] * 2.0, p_im[i] * 2.0);
            let square_re = den_re * den_re - den_im[i] * den_im[i];
            let square_im = den_re * den_im[i] + den_im[i] * den_re;
            let norm = square_re * square_re + square_im * square_im;
            let first_re = (twice_re * square_re + twice_im * square_im) / norm;
            let first_im = (twice_im * square_re - twice_re * square_im) / norm;
            // second = (1/√π + 0i) / den
            let norm = den_re * den_re + den_im[i] * den_im[i];
            let second_re = (inverse_root_pi * den_re + 0.0 * den_im[i]) / norm;
            let second_im = (0.0 * den_re - inverse_root_pi * den_im[i]) / norm;
            w.push(Complex64::new(first_re + second_re, first_im + second_im));
        }
    }
}

/// `faddeeva_along_x` が一度に並べる点の数（AVX-512 の 8 つの f64 の 4 本。Horner 法の 1 項の遅延を隠せる数）
const FADDEEVA_CHUNK: usize = 32;

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
