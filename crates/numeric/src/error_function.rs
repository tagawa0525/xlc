//! 誤差関数 erf、相補誤差関数 erfc、スケールした相補誤差関数 erfcx（ADR-0211）。xs（共鳴の揺らぎ、自由ガスの核、核分裂スペクトル）と
//! lattice（UFG の自由ガスの核）が使う。

/// 誤差関数 erf(x)。
pub fn erf(x: f64) -> f64 {
    erf_given(x, || (-x * x).exp())
}

/// 誤差関数 erf(x) を、呼び出し側が求めた `gaussian` = `(-x * x).exp()` を使って求める（[`erf`] とビットで同じ値）。同じ
/// exp(−x²) をほかの項にも掛ける呼び出し側が、指数関数を 1 度で済ませるため。符号を反転した x の (−x)·x は符号だけが違う
/// 積なので、x の符号によらず同じビットになる
pub fn erf_with_gaussian(x: f64, gaussian: f64) -> f64 {
    erf_given(x, || gaussian)
}

/// erf(x)。exp(−x²) は要る区間でだけ `gaussian` で求める（x ≥ [`SATURATION`] では求めない）
fn erf_given(x: f64, gaussian: impl FnOnce() -> f64) -> f64 {
    if x < 0.0 {
        return -erf_given(-x, gaussian);
    }
    if x < CODY_LOWER_LIMIT {
        erf_series(x, gaussian())
    } else if x >= SATURATION {
        1.0
    } else {
        // exp(−x*x) の x*x の丸めの効きは erfc(x) 倍に縮むので、erfc の exp_square で分けない。
        // NaN もここで NaN になる
        1.0 - erfcx(x) * gaussian()
    }
}

/// これより大きい x の erf(x) は倍精度で 1（erfc(6) ≈ 2e−17 は 1 の丸めの半分より小さいので、1 − erfc の計算も 1 になる）。
/// erfcx と指数関数を省く
const SATURATION: f64 = 6.0;

/// erf(x) = 2/√π exp(-x²) Σ 2ⁿ x^(2n+1) / (1·3·…·(2n+1))。項がすべて正なので桁落ちしない。`gaussian` は exp(−x²)
fn erf_series(x: f64, gaussian: f64) -> f64 {
    let x2 = x * x;
    let mut term = x;
    let mut sum = x;
    let mut n = 0.0;
    while term > sum * f64::EPSILON {
        n += 1.0;
        term *= 2.0 * x2 / (2.0 * n + 1.0);
        sum += term;
    }
    2.0 / std::f64::consts::PI.sqrt() * gaussian * sum
}

/// 相補誤差関数 erfc(x) = 1 - erf(x)。
pub fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        2.0 - erfc(-x)
    } else if x < CODY_LOWER_LIMIT {
        1.0 - erf_series(x, (-x * x).exp())
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
        return (x * x).exp() * (1.0 - erf_series(x, (-x * x).exp()));
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

/// exp(sign·x²) を、x² の丸めを指数に入れずに求める。x を 2⁻²¹ の倍数 s に切り捨てると、|x| < 32 で s の仮数は 26 ビット以下で s² は
/// 丸めなしで表せ、残りの x² − s² = (x − s)(x + s) は 2⁻²¹ × 2|x| より小さいので、その丸めは指数で増えない。|x| ≥ 32 では
/// exp(±s²) が 0 か桁あふれになり、そのまま返す（残りを掛けると 0 × ∞ になりうる）。exp(±x*x) と書かないのは、x*x の丸めが
/// 指数で x² 倍になるため。切る細かさは ADR-0211
fn exp_square(x: f64, sign: f64) -> f64 {
    const SCALE: f64 = 2_097_152.0; // 2²¹
    let s = (x * SCALE).trunc() / SCALE;
    let head = (sign * s * s).exp();
    if head == 0.0 || head.is_infinite() {
        return head;
    }
    let rest = (x - s) * (x + s);
    head * (sign * rest).exp()
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
/// 丸めが積もる（ADR-0211）。P と Q を 1 つのループで進めるのは、独立な 2 つの依存の列を並べて待ちを隠すため
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

#[cfg(test)]
mod tests;
