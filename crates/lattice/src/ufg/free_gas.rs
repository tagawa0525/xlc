//! UFG の減速計算の、標的の熱運動を温度 T の自由ガスで扱う散乱核。分離共鳴を持たない核種（ADR-0203）と、選んだときの分離共鳴を
//! 持つ核種（ADR-0204）に使う。選んだときの S(α, β) の核種の、S(α, β) の上限より下の入射は、ライブラリの S(α, β) の移行で
//! 扱う（[`SabSource`]、ADR-0205）。
//!
//! 入射 E'（レサジー u'）、出射 E（u）のレサジーあたりの散乱核は R(u'→u) = E S(E'→E)/N(E')。S は断面積を一定とした
//! 自由ガスの核（[`kernel`]）、N は温度に拡げた断面積と自由原子の断面積の比（[`effective_ratio`]）で、ライブラリの
//! ポイントワイズの弾性散乱断面積（温度に拡げたもの）に R を掛けると、自由原子の断面積 × S になる。範囲の中で
//!
//!   R = R_static/N(E') + D、R_static = e^{−(u−u')}/(1 − α)（u − Δ ≤ u' ≤ u）
//!
//! と分けると、R_static/N は静止した標的の窓の積分（弾性散乱断面積を σs/N にする）で扱え、D は静止した核の両端の近くに
//! だけ残る: 上端 E' = E の近く（上方散乱の e^{−(E−E')/kT} の裾と、下方散乱の側の幅 4kT/A ほどの折れ）と、下端 E' = E/α の
//! 近く（erf の引数 θ√(E/kT) − ζ√(E'/kT) が 0 を横切るぼけ。E' の幅 2√(E' kT)/ζ）。D は、その近くの区間（[`Gas::cells`]）
//! ごとに 3 点の 2 次式で近似し、中性子束の累積の積率（[`Prefix`]）で積分する。区間の外の D は、上端の側で e^{−[`TAIL`]}、
//! 下端の側で erfc([`EDGE`]) 程度より小さい（`correction_vanishes_outside_its_cells`）。

use std::f64::consts::PI;
use std::time::{Duration, Instant};

use rayon::prelude::*;
use xlc_format::library::SabTransfer;

use super::UfgError;
use xlc_numeric::quadrature::{GAUSS_LEGENDRE_2, GAUSS_LEGENDRE_4};

/// ボルツマン定数 (eV/K)。CODATA 2018 の値（xs の `physics::BK` と同じ）
pub(super) const BOLTZMANN: f64 = 8.617_333_262e-5;

/// 上端の側の区間の幅（kT の倍数）。上方散乱の裾 e^{−(E−E')/kT} がここで e^{−28} ≈ 7e−13 に下がる
pub(super) const TAIL: f64 = 28.0;

/// D の散乱源を並列に積分する区間の塊の大きさ（1 つの区間は核の評価 2 回とリングごとの積分で 1 µs ほど。塊を小さくすると
/// 並列の手間が勝つ）
const CELLS_PER_TASK: usize = 32;

/// 下端の側の区間の半幅（ぼけの幅 2√(E' kT)/ζ の倍数）。erfc(6) ≈ 2e−17
pub(super) const EDGE: f64 = 6.0;

/// 移行を数えるとき（[`Gas::share_cells`]）の、[`Gas::shares`] の区間の細かさ。核を直接積分した値に 1e−10 で合う細かさ
/// （`correction_shares_split_the_kernel_at_an_outgoing_energy`。減速計算の速い方法の 2 では、H の 4 eV で 6e−10 違う）
pub(super) const SHARE_CELLS: usize = 8;

/// 上端の下方散乱の側（E' > E）の D の大きさ e^{−E/(A kT)} の指数の上限。e^{−37} ≈ 9e−17 より小さい D は倍精度で R に
/// 埋もれるので、E/(A kT) がこれを超える格子点ではその側に区間を置かない
const DOWNSIDE: f64 = 37.0;

/// 上端の側の区間の幅を、端から y（kT の倍数）離れたところで e^{y/GRADE} 倍に広げる。裾は e^{−y} で減るので、
/// 2 次式の補間の誤差（区間の幅の 4 乗 × e^{−y}）が端から離れても増えず、区間の数は y によらず 2 GRADE × 区間の細かさで済む
const GRADE: f64 = 4.0;

/// 誤差関数 erf(x)。xs の `special::erf` と同じ計算（ステージのクレートは互いに依存せず、共通のクレートに移すと xs と
/// lattice の実行ファイルが変わり、すべてのステージが古くなるので、ここに複製する。ADR-0203）。
pub(super) fn erf(x: f64) -> f64 {
    if x < 0.0 {
        return -erf(-x);
    }
    if x < SERIES_LIMIT {
        erf_series(x)
    } else if x < SATURATION {
        1.0 - erfcx(x) * (-x * x).exp()
    } else {
        1.0
    }
}

/// これより大きい x の erf(x) は倍精度で 1（erfc(6) ≈ 2e−17 は 1 の丸めの半分より小さいので、1 − erfc の計算も 1 になる）。
/// erfcx と指数関数を省く
const SATURATION: f64 = 6.0;

/// これより小さい x の erf は級数で、大きい x の erf は [`SATURATION`] まで 1 − erfcx(x) exp(−x²) で求める
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
    2.0 / PI.sqrt() * (-x2).exp() * sum
}

/// スケールした相補誤差関数 erfcx(x) = exp(x²) erfc(x)（x ≥ 0）。xs の `special::erfcx` と同じ計算。
pub(super) fn erfcx(x: f64) -> f64 {
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

// Cody の有理近似の係数（ADR-0211）。W. J. Cody, "Rational Chebyshev approximations for the error function",
// Math. Comp. 23 (1969) 631–637 の近似を、Netlib の specfun の CALERF（倍精度の DATA 文。1990-03-19 版）からそのままコピーし、
// 高次の係数から並べた。コピーした値は、多倍長（mpmath）で有理式を評価し、erfcx との相対差が 0.46875 ≤ x ≤ 4 で 1.2e−18、
// 4 ≤ x ≤ 4e4 で 5.9e−19 以下になることで確かめた。出典と照らせるように、f64 で表せる桁より長い値も丸めずに残す

/// これより小さい x では erfcx を級数の erf から、これ以上の x では Cody の有理近似で求める（CALERF の区間の境 THRESH）。
/// 級数の 1 − erf は x が大きいほど桁落ちし、x = 2.5 の直下で 1 万 ulp を超える（ADR-0211）
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

/// 自由原子の散乱断面積を 1 としたときの、断面積を一定とした自由ガスのエネルギー移行の確率密度 S(E→E') (1/eV)。`awr` は
/// 標的の質量比、`kt` は標的の温度の kT (eV)。重心系で等方な散乱の、自由ガスの散乱則
///   S(α, β) = e^{−(α+β)²/(4α)}/√(4πα)（α, β は kT で測る。ADR-0133 の式と同じ）
/// を散乱角について積分した S(E→E') = ((A + 1)/A)² A/(4E) ∫ e^{−β/2} S(α, β) dα（α は (√E ∓ √E')²/(A kT) の間、
/// β = (E' − E)/kT。xs の `sab` の P0 の核と同じ形）の閉じた式
///   S = θ²/(2E) [erf(θ√ε' − ζ√ε) ± erf(θ√ε' + ζ√ε) + e^{ε−ε'} (erf(θ√ε − ζ√ε') ∓ erf(θ√ε + ζ√ε'))]、
///   θ = (A + 1)/(2√A)、ζ = (A − 1)/(2√A)、ε = E/kT、上の符号は E' ≤ E
/// （散乱則の α の積分と一致することは `kernel_is_the_integral_of_the_free_gas_law` で確かめる）。xs の `free_gas_kernel` と
/// 同じ計算で、複製する理由は [`erf`] と同じ。
pub(super) fn kernel(e: f64, e_out: f64, awr: f64, kt: f64) -> f64 {
    if e <= 0.0 || e_out <= 0.0 {
        return 0.0;
    }
    let sqrt_a = awr.sqrt();
    let theta = (awr + 1.0) / (2.0 * sqrt_a);
    let zeta = (awr - 1.0) / (2.0 * sqrt_a);
    let (s, s_out) = ((e / kt).sqrt(), (e_out / kt).sqrt());
    // e^(-y²) erfcx(x) = e^(-y²+x²) erfc(x)。指数を先にまとめて桁あふれを避ける。e^(-y²) が 0 に落ちるなら erfcx を
    // 求めない（erfcx(x) ≤ 1 なので積は 0）。`pair` は同じ指数因子の 2 つの erfcx の差
    let pair = |y: f64, a: f64, b: f64| {
        let factor = (-y * y).exp();
        if factor == 0.0 {
            0.0
        } else {
            factor * (erfcx(a) - erfcx(b))
        }
    };
    let scaled = |y: f64, x: f64| {
        let factor = (-y * y).exp();
        if factor == 0.0 {
            0.0
        } else {
            factor * erfcx(x)
        }
    };
    let bracket = if e_out <= e {
        // erf(θ√ε' - ζ√ε) + erf(θ√ε' + ζ√ε) + e^(ε-ε') [erfc(θ√ε + ζ√ε') - erfc(θ√ε - ζ√ε')]
        let (x1, x2) = (theta * s + zeta * s_out, theta * s - zeta * s_out);
        let (u, v) = (zeta * s - theta * s_out, zeta * s + theta * s_out);
        if u > 0.0 {
            // 遠い下方散乱: O(1) の項が打ち消し合うので、同じ指数因子をもつ項どうしをまとめる
            pair(u, u, x2) + pair(v, x1, v)
        } else {
            erf(-u) + erf(v) + scaled(v, x1) - scaled(u, x2)
        }
    } else {
        // erf(θ√ε' - ζ√ε) - erf(θ√ε' + ζ√ε) + e^(ε-ε') [erf(θ√ε - ζ√ε') + erf(θ√ε + ζ√ε')]
        let (p, q) = (theta * s_out - zeta * s, theta * s_out + zeta * s);
        let (r, t) = (theta * s - zeta * s_out, theta * s + zeta * s_out);
        if r < 0.0 {
            // 遠い上方散乱: 同様に、同じ指数因子をもつ項どうしをまとめる
            pair(p, -r, p) + pair(q, q, t)
        } else {
            scaled(q, q) - scaled(p, p) + ((e - e_out) / kt).exp() * (erf(r) + erf(t))
        }
    };
    // 打ち消し合いの丸めで負になった裾は 0 にする
    (theta * theta / (2.0 * e) * bracket).max(0.0)
}

/// ∫ S(E→E') dE' = N(E) = (1 + 1/(2a²)) erf(a) + e^{−a²}/(√π a)、a² = A E/kT。断面積を一定とした自由原子の断面積を
/// 温度に拡げた比（ドップラー拡がりの式で、断面積が一定のときの値）。
pub(super) fn effective_ratio(e: f64, awr: f64, kt: f64) -> f64 {
    let a2 = awr * e / kt;
    let a = a2.sqrt();
    (1.0 + 0.5 / a2) * erf(a) + (-a2).exp() / (PI.sqrt() * a)
}

/// 自由ガスの核で散乱する 1 つの核種（名前と温度の組）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Gas {
    awr: f64,
    kt: f64,
    /// 1/(1 − α)
    gain: f64,
    /// Δ = ln(1/α)
    width: f64,
    /// |ζ| = |A − 1|/(2√A)
    zeta: f64,
    /// この入射エネルギー (eV) 未満で自由ガスの核を使う
    cap: f64,
    /// この入射エネルギー (eV) 未満では、自由ガスの核の代わりに S(α, β) の移行を使う（[`Gas::below_sab`]。0 は使わない）
    floor: f64,
}

impl Gas {
    /// 質量比 `awr`、温度 `temperature` (K) の核種。`upper` は自由ガスの核を使う入射エネルギーの上限（kT の倍数）で、
    /// `None` は上限なし。
    pub(super) fn new(awr: f64, temperature: f64, upper: Option<f64>) -> Self {
        let kt = BOLTZMANN * temperature;
        let alpha = ((awr - 1.0) / (awr + 1.0)).powi(2);
        Self {
            awr,
            kt,
            gain: 1.0 / (1.0 - alpha),
            width: -alpha.ln(),
            zeta: (awr - 1.0).abs() / (2.0 * awr.sqrt()),
            cap: upper.map_or(f64::INFINITY, |n| n * kt),
            floor: 0.0,
        }
    }

    /// S(α, β) の移行の表 `table` の上限（入射の最後の点）より下の入射で、自由ガスの核の代わりにその表を使う核種（ADR-0205）。
    /// 静止した標的の窓は上限より下で 0 にし（[`Gas::window_scale`]）、自由ガスの核の差 D の区間は上限より上の入射に限る
    pub(super) fn below_sab(self, table: &SabTransfer) -> Self {
        Self {
            floor: table.incident[table.incident.len() - 1],
            ..self
        }
    }

    /// S(α, β) の移行を使う入射エネルギーの上限 (eV)。使わなければ `None`
    pub(super) fn floor(&self) -> Option<f64> {
        (self.floor > 0.0).then_some(self.floor)
    }

    /// 自由ガスの核か S(α, β) の移行を使う入射エネルギーの上限 (eV)
    pub(super) fn top(&self) -> f64 {
        self.cap.max(self.floor)
    }

    /// 自由ガスの核を使う入射エネルギーの上限 (eV)
    pub(super) fn cap(&self) -> f64 {
        self.cap
    }

    /// kT (eV)
    pub(super) fn kt(&self) -> f64 {
        self.kt
    }

    /// 入射エネルギー `e_in` の静止した標的の窓に掛ける 1/N(E')（自由ガスの核を使わない入射エネルギーでは 1、S(α, β) の移行を使う
    /// 入射エネルギーでは 0。`below` は不連続点の低い側の極限で、S(α, β) の上限ちょうどの低い側も 0）
    pub(super) fn window_scale(&self, e_in: f64, below: bool) -> f64 {
        if e_in < self.floor || (below && e_in == self.floor) {
            0.0
        } else if e_in < self.cap {
            1.0 / effective_ratio(e_in, self.awr, self.kt)
        } else {
            1.0
        }
    }

    /// D(u'→u) = R(u'→u) − R_static(u'→u)/N(E')（入射 `x` = u'、出射 `u` とそのエネルギー `e_out`、`high` は範囲の上端）。
    /// `inside` は u' を静止した核の窓（u − Δ < u' < u）の中として扱う（窓の端では、区間の側の極限を選ぶ）。
    pub(super) fn correction(&self, high: f64, x: f64, u: f64, e_out: f64, inside: bool) -> f64 {
        let e_in = high * (-x).exp();
        let ratio = effective_ratio(e_in, self.awr, self.kt);
        let thermal = e_out * kernel(e_in, e_out, self.awr, self.kt) / ratio;
        if inside {
            thermal - self.gain * (-(u - x)).exp() / ratio
        } else {
            thermal
        }
    }

    /// 出射のレサジー `u`（エネルギー `e`）で D が 0 でない区間を、昇順に `out` に書く（`bounds` は作業用）。上端の側は
    /// u の両側に、E' = E ± y kT（y は `grades` の段、[`grades`]）の点を置き、下端の側は u − Δ を中心に、ぼけの幅の
    /// 1/`cells` の等間隔で ±[`EDGE`] 幅を覆う。両方を合わせ、u、u − Δ、自由ガスの核の上限、S(α, β) の上限、範囲の上端
    /// （u' = 0）、`bottom`（格子の最後の点）で区切る（D と中性子束の扱いが、区間の中で変わらないように）。範囲の上（u' < 0）と、
    /// 自由ガスの核を使わない入射エネルギーの区間（上限より上は静止した標的、S(α, β) の上限より下は S(α, β) の移行で扱う）は落とす。
    #[allow(clippy::too_many_arguments)]
    pub(super) fn cells(
        &self,
        high: f64,
        u: f64,
        e: f64,
        grades: &[f64],
        cells: usize,
        bottom: f64,
        bounds: &mut Vec<f64>,
        out: &mut Vec<(f64, f64)>,
    ) {
        bounds.clear();
        out.clear();
        let scale = self.kt / e;
        // E' = E ± y kT のレサジー（上方散乱の裾は E − E' で減るので、エネルギーで測る）。E' が 0 に届く段は打ち切る
        let below = |y: f64| u - (-y * scale).ln_1p();
        let above = |y: f64| u - (y * scale).ln_1p();
        // E' = E − y kT が 0 に近づく段（y kT が E の 0.999 倍を超える）は置かない（ln(1 − y kT/E) が発散する）。範囲の下限 4 eV
        // では、TAIL kT が 4 eV を超える約 1660 K より高い温度で、この打ち切りが裾を落とす
        let reach = TAIL.min(0.999 / scale);
        // 下方散乱の側の折れの大きさは e^{−E/(A kT)} ほどで、DOWNSIDE より小さければ区間を置かない
        let downside = if e < DOWNSIDE * self.awr * self.kt {
            TAIL
        } else {
            0.0
        };
        let mut supports = [(above(downside), below(reach)), (f64::NAN, f64::NAN)];
        for &y in grades {
            if y <= downside {
                bounds.push(above(y));
            }
            if y <= reach {
                bounds.push(below(y));
            }
        }
        bounds.push(below(reach));
        let center = u - self.width;
        let e_edge = high * (-center).exp();
        let blur = 2.0 * (self.kt / e_edge).sqrt() / self.zeta;
        if blur.is_finite() && center + EDGE * blur > 0.0 {
            let n = (2.0 * EDGE * cells as f64).ceil() as usize;
            let step = 2.0 * EDGE * blur / n as f64;
            for j in 0..=n {
                bounds.push(center - EDGE * blur + j as f64 * step);
            }
            supports[1] = (center - EDGE * blur, center + EDGE * blur);
        }
        let covered = |x: f64| supports.iter().any(|&(lo, hi)| lo <= x && x <= hi);
        let floor = if self.floor > 0.0 {
            (high / self.floor).ln()
        } else {
            f64::INFINITY
        };
        let breaks = [center, (high / self.cap).ln(), floor, 0.0, bottom];
        bounds.extend(breaks.into_iter().filter(|&x| covered(x)));
        bounds.sort_by(f64::total_cmp);
        // 丸めで重なった境界（下端の側の等間隔の点と u − Δ など）は 1 つにまとめる。幅がほぼ 0 の区間の 2 次式は、3 点の
        // 値の丸めの差を幅の 2 乗で割って発散する。まとめると、区切りが区間の幅の 1e−9 倍（u の丸めの 1e3 倍より狭い区間も）
        // だけずれる
        let tiny = (1e-9 * scale.min(if blur.is_finite() { blur } else { scale }) / cells as f64)
            .max(1e3 * f64::EPSILON * u.abs());
        bounds.dedup_by(|b, a| *b - *a <= tiny);
        for w in bounds.windows(2) {
            let (a, b) = (w[0], w[1]);
            let mid = 0.5 * (a + b);
            let e_mid = high * (-mid).exp();
            if b > a && mid >= 0.0 && e_mid < self.cap && e_mid >= self.floor && covered(mid) {
                out.push((a, b));
            }
        }
    }

    /// 出射のレサジー `u`（エネルギー `e`）で、[`Gas::cells`] の区間がどれも上限より上の入射にあって落ちるか。上端の側の
    /// 区間の最も低い入射は E − TAIL kT、下端の側は u − Δ + EDGE × ぼけの幅（[`Gas::cells`] と同じ式。O-16 のように下端の側の
    /// ぼけが窓の幅を超える核種では、上方散乱の裾より下まで届く）
    pub(super) fn idle(&self, high: f64, u: f64, e: f64) -> bool {
        if e - TAIL * self.kt < self.cap {
            return false;
        }
        let center = u - self.width;
        let blur = 2.0 * (self.kt / (high * (-center).exp())).sqrt() / self.zeta;
        !blur.is_finite() || high * (-(center + EDGE * blur)).exp() >= self.cap
    }

    /// 入射 `e_in` の D を、出射のエネルギー `e_x` の下と上に分けて積分した [∫_{E<e_x} D_E dE, ∫_{E>e_x} D_E dE]。D_E は出射の
    /// エネルギーあたりの D（D_E(E) = (S(E'→E) − S_static(E))/N(E')、S_static = 1/((1 − α)E') は αE' ≤ E ≤ E' の静止した核）で、
    /// 2 つの和は 1 − 1/N（∫ S dE = N、∫ S_static dE = 1）。D_E が 0 でないのは、上端の側 E' ± TAIL kT と、下端の側の核の erf の
    /// 引数 θ√(E/kT) − |ζ|√(E'/kT) が ±EDGE の間（A ≈ 1 の H では E = 0 から 36 kT ほど。分母が θ なので、[`Gas::cells`] の入射の
    /// 側のぼけの幅と違う）だけ。上端の側は [`grades`] の段（`grades` は `cells` で作ったもの）、下端の側は √E の等間隔の
    /// 2 EDGE `cells` 個の区間で区切り、αE'、E'、`e_x` でも区切る（S_static の飛びと核の折れ）。区間ごとに E = s² kT で s について
    /// 4 点の Gauss-Legendre 則で積分する（H の核は E → 0 で √E のように立ち上がり、s については滑らか）
    pub(super) fn shares(&self, e_in: f64, e_x: f64, grades: &[f64], cells: usize) -> [f64; 2] {
        let kt = self.kt;
        let ratio = effective_ratio(e_in, self.awr, kt);
        let alpha = 1.0 - 1.0 / self.gain;
        let edge = alpha * e_in;
        let level = self.gain / e_in;
        let density = |e: f64| {
            let s_static = if edge <= e && e <= e_in { level } else { 0.0 };
            (kernel(e_in, e, self.awr, kt) - s_static) / ratio
        };
        let reach = TAIL.min(e_in / kt);
        let mut bounds: Vec<f64> = Vec::with_capacity(2 * grades.len() + 8);
        for &y in grades {
            bounds.push(e_in + y * kt);
            if y <= reach {
                bounds.push(e_in - y * kt);
            }
        }
        bounds.push(e_in - reach * kt);
        let upper = (e_in - reach * kt, e_in + TAIL * kt);
        let theta = (self.awr + 1.0) / (2.0 * self.awr.sqrt());
        let center = self.zeta * (e_in / kt).sqrt();
        let (s0, s1) = (((center - EDGE) / theta).max(0.0), (center + EDGE) / theta);
        let pieces = (2.0 * EDGE * cells as f64).ceil() as usize;
        for j in 0..=pieces {
            let s = s0 + (s1 - s0) * j as f64 / pieces as f64;
            bounds.push(s * s * kt);
        }
        let lower = (s0 * s0 * kt, s1 * s1 * kt);
        bounds.extend([edge, e_in, e_x]);
        bounds.retain(|b| b.is_finite() && *b >= 0.0);
        bounds.sort_by(f64::total_cmp);
        bounds.dedup();
        let covered = |x: f64| (upper.0 <= x && x <= upper.1) || (lower.0 <= x && x <= lower.1);
        let mut shares = [0.0; 2];
        for w in bounds.windows(2) {
            let (a, b) = (w[0], w[1]);
            let mid = 0.5 * (a + b);
            if b <= a || !covered(mid) {
                continue;
            }
            let (sa, sb) = ((a / kt).sqrt(), (b / kt).sqrt());
            let (center, half) = (0.5 * (sa + sb), 0.5 * (sb - sa));
            let value: f64 = GAUSS_LEGENDRE_4
                .iter()
                .map(|&(x, w)| {
                    let s = center + half * x;
                    w * half * density(s * s * kt) * 2.0 * s * kt
                })
                .sum();
            shares[usize::from(mid > e_x)] += value;
        }
        shares
    }

    /// 出射の境界（エネルギー `e_b`、レサジー `x`。`high` は範囲の上端）について、[`Gas::shares`] が 0 でない入射のレサジーを、
    /// 自由ガスの核を使う入射のレサジー `incident`（[始め, 終わり]）の中で区間に分け、区間ごとに割合を 2 次式で近似したもの
    /// [境界より上の入射の、境界より下への割合, 境界より下の入射の、境界より上への割合]。D の山が境界にかかる入射は、上端の側が
    /// E' = `e_b` ± y kT（y は `grades` の段）、下端の側が erf の引数が境界で ±EDGE になる √E' の間（2 EDGE `cells` 個の等間隔）で、
    /// そこで区切る。それ以外の入射の割合は、下端の側の山の全体の積分か 0 で、入射とともに滑らかに変わるので、幅 `width` 以下の
    /// レサジーの等間隔で区切る
    #[allow(clippy::too_many_arguments)]
    pub(super) fn share_cells(
        &self,
        high: f64,
        x: f64,
        e_b: f64,
        incident: [f64; 2],
        grades: &[f64],
        cells: usize,
        width: f64,
    ) -> [Vec<ShareCell>; 2] {
        let kt = self.kt;
        let theta = (self.awr + 1.0) / (2.0 * self.awr.sqrt());
        let s_b = (e_b / kt).sqrt();
        // 下端の側の山 [((|ζ|√ε' − EDGE)⁺/θ)², ((|ζ|√ε' + EDGE)/θ)²] kT が境界の下にかかる入射は √ε' < (θ s_b + EDGE)/|ζ|、上に
        // かかる入射は √ε' > (θ s_b − EDGE)/|ζ|
        let below_reach = if self.zeta > 0.0 {
            ((theta * s_b + EDGE) / self.zeta).powi(2) * kt
        } else {
            f64::INFINITY
        };
        let above_reach = if theta * s_b <= EDGE {
            0.0
        } else if self.zeta > 0.0 {
            ((theta * s_b - EDGE) / self.zeta).powi(2) * kt
        } else {
            f64::INFINITY
        };
        let lethargy = |e: f64| {
            if e > 0.0 {
                (high / e).ln()
            } else {
                f64::INFINITY
            }
        };
        let spans = [
            [
                incident[0].max(lethargy((e_b + TAIL * kt).max(below_reach))),
                incident[1].min(x),
            ],
            [
                incident[0].max(x),
                incident[1].min(lethargy((e_b - TAIL * kt).min(above_reach))),
            ],
        ];
        let mut breaks: Vec<f64> = Vec::with_capacity(2 * grades.len() + 4 * EDGE as usize * cells);
        for &y in grades {
            breaks.push(lethargy(e_b + y * kt));
            breaks.push(lethargy(e_b - y * kt));
        }
        if self.zeta > 0.0 {
            let (lo, hi) = (
                ((theta * s_b - EDGE) / self.zeta).max(0.0),
                (theta * s_b + EDGE) / self.zeta,
            );
            let pieces = (2.0 * EDGE * cells as f64).ceil() as usize;
            for j in 0..=pieces {
                let s = lo + (hi - lo) * j as f64 / pieces as f64;
                breaks.push(lethargy(s * s * kt));
            }
        }
        let shares = |v: f64, side: usize| self.shares(high * (-v).exp(), e_b, grades, cells)[side];
        [0, 1].map(|side| {
            let [a, b] = spans[side];
            let mut out = Vec::new();
            if a.is_nan() || b.is_nan() || a >= b {
                return out;
            }
            let mut nodes: Vec<f64> = breaks
                .iter()
                .copied()
                .filter(|&v| a < v && v < b)
                .chain([a, b])
                .collect();
            nodes.sort_by(f64::total_cmp);
            nodes.dedup();
            let mut start = (a, shares(a, side));
            for w in nodes.windows(2) {
                let (lo, hi) = (w[0], w[1]);
                let n = ((hi - lo) / width).ceil().max(1.0) as usize;
                for j in 1..=n {
                    let end = if j == n {
                        hi
                    } else {
                        lo + (hi - lo) * j as f64 / n as f64
                    };
                    let mid = 0.5 * (start.0 + end);
                    let finish = (end, shares(end, side));
                    out.push(ShareCell {
                        a: start.0,
                        b: end,
                        quadratic: Quadratic::through(
                            start.0,
                            end,
                            start.1,
                            shares(mid, side),
                            finish.1,
                        ),
                    });
                    start = finish;
                }
            }
            out
        })
    }

    /// 区間 [a, b] を静止した核の窓の中として扱うか
    pub(super) fn inside(&self, u: f64, a: f64, b: f64) -> bool {
        let mid = 0.5 * (a + b);
        u - self.width < mid && mid < u
    }
}

/// 上端の側の区間の段: 0 から幅 1/`cells`（kT の倍数）で始め、端から y 離れたところで e^{y/GRADE} 倍に広げて [`TAIL`] まで。
pub(super) fn grades(cells: usize) -> Vec<f64> {
    let h = 1.0 / cells as f64;
    let mut y = 0.0;
    let mut out = vec![0.0];
    while y < TAIL {
        y = (y + h * (y / GRADE).exp()).min(TAIL);
        out.push(y);
    }
    out
}

/// 倍精度の 2 つの和で表した数（累積の積率の桁落ちを防ぐ）。
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(super) struct Dd {
    hi: f64,
    lo: f64,
}

fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}

fn fast_two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    (s, b - (s - a))
}

impl Dd {
    fn add(self, x: f64) -> Self {
        let (s, e) = two_sum(self.hi, x);
        let (hi, lo) = fast_two_sum(s, e + self.lo);
        Self { hi, lo }
    }

    fn add_dd(self, other: Self) -> Self {
        let (s, e) = two_sum(self.hi, other.hi);
        let (hi, lo) = fast_two_sum(s, e + self.lo + other.lo);
        Self { hi, lo }
    }

    fn sub(self, other: Self) -> Self {
        self.add_dd(Self {
            hi: -other.hi,
            lo: -other.lo,
        })
    }

    fn mul(self, x: f64) -> Self {
        let p = self.hi * x;
        let e = self.hi.mul_add(x, -p);
        let (hi, lo) = fast_two_sum(p, e + self.lo * x);
        Self { hi, lo }
    }

    fn value(self) -> f64 {
        self.hi + self.lo
    }
}

/// [u_a, u_b] で線形な F（両端の値 fa, fb）の ∫ F s^j ds（j = 0, 1, 2）。2 点の Gauss-Legendre 則は 3 次まで厳密
fn segment_moments(ua: f64, ub: f64, fa: f64, fb: f64) -> [f64; 3] {
    let (mid, half) = (0.5 * (ua + ub), 0.5 * (ub - ua));
    let mut m = [0.0; 3];
    if half == 0.0 {
        return m;
    }
    for &(x, w) in &GAUSS_LEGENDRE_2 {
        let s = mid + half * x;
        let f = fa + (fb - fa) * 0.5 * (x + 1.0);
        let v = w * half * f;
        m[0] += v;
        m[1] += v * s;
        m[2] += v * s * s;
    }
    m
}

/// 1 つのリングの中性子束の累積の積率 P_j(u_q) = ∫_{u_start}^{u_q} F s^j ds（j = 0, 1, 2。F は格子点の間で u について線形）。
#[derive(Debug, Clone, Default, PartialEq)]
pub(super) struct Prefix {
    start: usize,
    values: Vec<[Dd; 3]>,
}

impl Prefix {
    pub(super) fn new(start: usize) -> Self {
        Self {
            start,
            values: Vec::new(),
        }
    }

    /// 格子点 m（= 次の点）の累積の積率を足す（`flux` は m までの値を持つ）
    pub(super) fn push(&mut self, u: &[f64], flux: &[f64], m: usize) {
        if m == self.start {
            self.values.clear();
            self.values.push([Dd::default(); 3]);
            return;
        }
        let last = self.values[m - 1 - self.start];
        let s = segment_moments(u[m - 1], u[m], flux[m - 1], flux[m]);
        self.values
            .push([last[0].add(s[0]), last[1].add(s[1]), last[2].add(s[2])]);
    }

    /// ∫_{x0}^{x1} F (s − c)^j ds（j = 0, 1, 2）。`q0`, `q1` は u_q ≤ x の格子点。端の区間の部分は c のまわりで直接
    /// 積分し、間の格子の区間だけを累積の積率の差にする（幅の小さい部分を、大きな累積の差で求めると桁落ちする）
    #[allow(clippy::too_many_arguments)]
    fn moments(
        &self,
        u: &[f64],
        flux: &[f64],
        q0: usize,
        x0: f64,
        q1: usize,
        x1: f64,
        c: f64,
    ) -> [f64; 3] {
        if q0 >= q1 {
            return local_moments(u, flux, q0, x0, x1, c);
        }
        let head = local_moments(u, flux, q0, x0, u[q0 + 1], c);
        let tail = if x1 > u[q1] {
            local_moments(u, flux, q1, u[q1], x1, c)
        } else {
            [0.0; 3]
        };
        let (p0, p1) = (
            self.values[q0 + 1 - self.start],
            self.values[q1 - self.start],
        );
        let d = [p1[0].sub(p0[0]), p1[1].sub(p0[1]), p1[2].sub(p0[2])];
        let m0 = d[0];
        let m1 = d[1].sub(d[0].mul(c));
        let m2 = d[2].sub(d[1].mul(2.0 * c)).add_dd(d[0].mul(c).mul(c));
        [
            head[0] + m0.value() + tail[0],
            head[1] + m1.value() + tail[1],
            head[2] + m2.value() + tail[2],
        ]
    }
}

/// 格子の区間 [u_q, u_{q+1}] の中の [x0, x1] での ∫ F (s − c)^j ds（j = 0, 1, 2。F は区間の両端の値の線形補間）。2 点の
/// Gauss-Legendre 則は 3 次まで厳密
fn local_moments(u: &[f64], flux: &[f64], q: usize, x0: f64, x1: f64, c: f64) -> [f64; 3] {
    let mut m = [0.0; 3];
    let width = u[q + 1] - u[q];
    if x1 <= x0 || width <= 0.0 {
        return m;
    }
    let (mid, half) = (0.5 * (x0 + x1), 0.5 * (x1 - x0));
    for &(x, w) in &GAUSS_LEGENDRE_2 {
        let s = mid + half * x;
        let f = flux[q] + (flux[q + 1] - flux[q]) * (s - u[q]) / width;
        let d = s - c;
        let v = w * half * f;
        m[0] += v;
        m[1] += v * d;
        m[2] += v * d * d;
    }
    m
}

/// 区間 [a, b] の 2 次式 f(s) = a0 + a1 (s − c) + a2 (s − c)²（両端と中点の値から）。
#[derive(Debug, Clone, Copy)]
pub(super) struct Quadratic {
    c: f64,
    coefficients: [f64; 3],
}

impl Quadratic {
    pub(super) fn through(a: f64, b: f64, fa: f64, fm: f64, fb: f64) -> Self {
        let c = 0.5 * (a + b);
        let h = 0.5 * (b - a);
        Self {
            c,
            coefficients: [
                fm,
                (fb - fa) / (2.0 * h),
                (fa + fb - 2.0 * fm) / (2.0 * h * h),
            ],
        }
    }

    pub(super) fn at(&self, s: f64) -> f64 {
        let d = s - self.c;
        let [a0, a1, a2] = self.coefficients;
        a0 + d * (a1 + d * a2)
    }

    /// ∫_{x0}^{x1} f ds
    fn integral(&self, x0: f64, x1: f64) -> f64 {
        let [a0, a1, a2] = self.coefficients;
        let (d0, d1) = (x0 - self.c, x1 - self.c);
        a0 * (x1 - x0) + a1 * 0.5 * (d1 * d1 - d0 * d0) + a2 * (d1 * d1 * d1 - d0 * d0 * d0) / 3.0
    }

    /// ∫ f F ds（積率 `m` は (s − c)^j の積分）
    fn against(&self, m: [f64; 3]) -> f64 {
        let [a0, a1, a2] = self.coefficients;
        a0 * m[0] + a1 * m[1] + a2 * m[2]
    }

    /// [x0, x1] で線形な F（両端の重み）について、∫ f F ds の F(x0) と F(x1) の係数
    fn linear(&self, x0: f64, x1: f64) -> (f64, f64) {
        let (mid, half) = (0.5 * (x0 + x1), 0.5 * (x1 - x0));
        let (mut wa, mut wb) = (0.0, 0.0);
        if half <= 0.0 {
            return (0.0, 0.0);
        }
        for &(x, w) in &GAUSS_LEGENDRE_2 {
            let f = w * half * self.at(mid + half * x);
            let t = 0.5 * (x + 1.0);
            wa += f * (1.0 - t);
            wb += f * t;
        }
        (wa, wb)
    }
}

/// 入射のレサジーの区間 [a, b] で、D を出射の境界の下（または上）に積分した割合を近似する 2 次式（[`Gas::share_cells`]）。
#[derive(Debug, Clone, Copy)]
pub(super) struct ShareCell {
    pub(super) a: f64,
    pub(super) b: f64,
    pub(super) quadratic: Quadratic,
}

/// 1 つのリングの、格子点 m での D の散乱源の計算に使う中性子束。
pub(super) struct History<'a> {
    /// 今のスイープの中性子束（m より前の格子点の値）と累積の積率
    pub(super) flux: &'a [f64],
    pub(super) prefix: &'a Prefix,
    /// 前のスイープの中性子束と累積の積率（全体）。最初のスイープは `None` で、m より後の中性子束を F_m とする
    pub(super) previous: Option<(&'a [f64], &'a Prefix)>,
}

/// 格子点 m の中性子束を F_m として、区間 [a, b] の ∫ f F ds を（既知の部分, F_m の係数）に分ける。`u` は格子の全体の
/// レサジー（昇順）、`qa`, `qb` は u_q ≤ a、u_q ≤ b の格子点（格子の最後の点より後は最後の点）。格子の最後の点より後の
/// 中性子束は、最後の点の値とする。
#[allow(clippy::too_many_arguments)]
pub(super) fn integrate(
    f: &Quadratic,
    u: &[f64],
    m: usize,
    a: f64,
    b: f64,
    qa: usize,
    qb: usize,
    history: &History,
) -> (f64, f64) {
    let last = u.len() - 1;
    let (mut known, mut own) = (0.0, 0.0);
    // 今のスイープの、u_{m−1} までの部分
    if m >= 1 {
        let x1 = b.min(u[m - 1]);
        if x1 > a {
            let m1 = history
                .prefix
                .moments(u, history.flux, qa, a, qb.min(m - 1), x1, f.c);
            known += f.against(m1);
        }
        // [u_{m−1}, u_m]: F_{m−1} は既知、F_m は未知
        let (x0, x1) = (a.max(u[m - 1]), b.min(u[m]));
        if x1 > x0 {
            let (wa, wb) = f.linear(x0, x1);
            // 区間の端の F の重みを、格子点 m−1 と m の重みに直す
            let (ta, tb) = (
                (x0 - u[m - 1]) / (u[m] - u[m - 1]),
                (x1 - u[m - 1]) / (u[m] - u[m - 1]),
            );
            let fm1 = history.flux[m - 1];
            known += fm1 * (wa * (1.0 - ta) + wb * (1.0 - tb));
            own += wa * ta + wb * tb;
        }
    }
    if m < last {
        // [u_m, u_{m+1}]: F_m は未知、F_{m+1} は前のスイープ（最初のスイープでは F_m）
        let (x0, x1) = (a.max(u[m]), b.min(u[m + 1]));
        if x1 > x0 {
            let (wa, wb) = f.linear(x0, x1);
            let (ta, tb) = (
                (x0 - u[m]) / (u[m + 1] - u[m]),
                (x1 - u[m]) / (u[m + 1] - u[m]),
            );
            match history.previous {
                Some((flux, _)) => {
                    own += wa * (1.0 - ta) + wb * (1.0 - tb);
                    known += flux[m + 1] * (wa * ta + wb * tb);
                }
                None => own += wa + wb,
            }
        }
        // [u_{m+1}, u_last]: 前のスイープ
        let (x0, x1) = (a.max(u[m + 1]), b.min(u[last]));
        if x1 > x0 {
            match history.previous {
                Some((flux, prefix)) => {
                    let q0 = qa.max(m + 1);
                    let q1 = qb.min(last);
                    known += f.against(prefix.moments(u, flux, q0, x0, q1, x1, f.c));
                }
                None => own += f.integral(x0, x1),
            }
        }
    }
    // 格子の最後の点より後: 最後の点の中性子束
    let x0 = a.max(u[last]);
    if b > x0 {
        let integral = f.integral(x0, b);
        match history.previous {
            Some((flux, _)) if m < last => known += flux[last] * integral,
            _ => own += integral,
        }
    }
    (known, own)
}

/// 格子の昇順のレサジー `u` で、u_q ≤ x となる最後の格子点（x が最初の点より前なら 0）
fn locate(u: &[f64], x: f64) -> usize {
    u.partition_point(|&v| v <= x).saturating_sub(1)
}

/// [`locate`] と同じ格子点を、u_hint ≤ x なら `hint` から間隔を倍にしながら先へ探す。区間の端は昇順に並び、前の端の格子点の
/// 近くにあるので、格子の全体（数十万点）を 2 分探索するより、キャッシュに載らない点を読まずに済む
fn locate_from(u: &[f64], x: f64, hint: usize) -> usize {
    if !(hint < u.len() && u[hint] <= x) {
        return locate(u, x);
    }
    // u[low] ≤ x を保ち、u[low + step] > x か格子の終わりになるまで進める
    let (mut low, mut step) = (hint, 1);
    loop {
        let high = (low + step).min(u.len());
        if high == u.len() || u[high] > x {
            return low + u[low..high].partition_point(|&v| v <= x) - 1;
        }
        low = high;
        step *= 2;
    }
}

/// 格子点の値 `values` を、レサジー x で u について線形に補間する（最後の点より後は最後の点の値）
fn interpolate(u: &[f64], values: &[f64], q: usize, x: f64) -> f64 {
    let last = u.len() - 1;
    if q >= last || x <= u[q] {
        return values[q.min(last)];
    }
    let width = u[q + 1] - u[q];
    if width <= 0.0 {
        return values[q + 1];
    }
    values[q] + (values[q + 1] - values[q]) * (x - u[q]) / width
}

/// 上限 `cap` (eV) より下の入射を補間する最初の格子点（`energies` はエネルギーの降順）。区間は上限で区切られ、上限より下の入射だけが
/// 残る（[`Gas::cells`]）ので、補間する格子の区間は上限を含む区間から先にある。上限が範囲の下限より下なら最後の格子点（範囲の
/// 下限より下の入射は、最後の格子点の値を使う）
pub(super) fn first_below(energies: &[f64], cap: f64) -> usize {
    let last = energies.len().saturating_sub(1);
    energies
        .iter()
        .position(|&e| e < cap)
        .map_or(last, |k| k.saturating_sub(1))
}

/// 減速計算のスイープの間の、自由ガスの核の D の散乱源の状態。
#[derive(Debug, Clone)]
pub(super) struct Thermal<'a> {
    /// [自由ガスの核種] の核
    gases: Vec<Gas>,
    /// [自由ガスの核種] の（最初の格子点、そこからの [格子点] のポイントワイズの弾性散乱断面積）（D に掛ける。最初の格子点は
    /// [`first_below`]）
    sigma: Vec<(usize, Vec<f64>)>,
    /// [自由ガスの核種] の（リング, 個数密度）
    rings: Vec<Vec<(usize, f64)>>,
    /// [自由ガスの核種] の、S(α, β) の上限より下の入射で使う S(α, β) の移行の表（[`Gas::below_sab`] の核種だけ）
    tables: Vec<Option<&'a SabTransfer>>,
    /// 上端の側の区間の段（[`grades`]）と、下端の側の区間の細かさ
    grades: Vec<f64>,
    cells: usize,
    /// D が 0 でない格子点のうち最初のものの 1 つ前（ここから掃き直す。これより前の格子点の D は 0）
    start: usize,
    /// [リング] の今のスイープの累積の積率
    prefix: Vec<Prefix>,
    /// 前のスイープの [リング] の中性子束と累積の積率
    previous: Option<(Vec<Vec<f64>>, Vec<Prefix>)>,
    /// 作業用の区間の境界と区間
    bounds: Vec<f64>,
    pieces: Vec<(f64, f64)>,
    /// 作業用の、すべての自由ガスの核種の（核種, 区間の両端）
    all: Vec<(usize, f64, f64)>,
    /// 作業用の、S(α, β) の移行の小区間と、その重み
    sab: SabWork,
    /// 散乱源の経過時間（ログの内訳）: [自由ガスの核の区間, S(α, β) の移行]
    elapsed: [Duration; 2],
}

impl<'a> Thermal<'a> {
    /// 自由ガスの核種の核 `gases`、その [格子点] の弾性散乱断面積 `sigma`、（リング, 個数密度）`rings`、S(α, β) の移行の表
    /// `tables`、区間の細かさ `cells` の状態。`energies` は格子点のエネルギー（降順）、`n` はリングの数。自由ガスの核種がなければ
    /// `None`。
    pub(super) fn new(
        gases: Vec<Gas>,
        sigma: Vec<(usize, Vec<f64>)>,
        rings: Vec<Vec<(usize, f64)>>,
        tables: Vec<Option<&'a SabTransfer>>,
        cells: usize,
        energies: &[f64],
        n: usize,
    ) -> Option<Self> {
        if gases.is_empty() {
            return None;
        }
        // 上限より上の入射エネルギーは静止した標的なので、上方散乱の裾が上限の下から届く格子点より前は D が 0。S(α, β) の移行は、
        // 表の上限から表の移行の最大まで上に届く
        let reach = gases
            .iter()
            .zip(&tables)
            .map(|(l, table)| {
                let sab = table
                    .as_ref()
                    .map_or(0.0, |t| l.top() + t.transfer[t.transfer.len() - 1]);
                (l.cap() + TAIL * l.kt()).max(sab)
            })
            .fold(0.0, f64::max);
        let first = energies.iter().position(|&e| e <= reach)?;
        let start = first.saturating_sub(1);
        Some(Self {
            gases,
            sigma,
            rings,
            tables,
            grades: grades(cells),
            cells,
            start,
            prefix: vec![Prefix::new(start); n],
            previous: None,
            bounds: Vec::new(),
            pieces: Vec::new(),
            all: Vec::new(),
            sab: SabWork::default(),
            elapsed: [Duration::ZERO; 2],
        })
    }

    pub(super) fn start(&self) -> usize {
        self.start
    }

    /// 自由ガスの核種 `l` の核と、[`Thermal::new`] に渡したその弾性散乱断面積の列（最初の格子点と、そこからの [格子点] の値）、
    /// S(α, β) の移行の表
    pub(super) fn gas(&self, l: usize) -> (&Gas, usize, &[f64], Option<&SabTransfer>) {
        let (offset, sigma) = &self.sigma[l];
        (&self.gases[l], *offset, sigma, self.tables[l])
    }

    /// 格子点 m（エネルギー `e`）の D の散乱源を、リングごとに既知の部分を `known` に、F_m の係数を `own` に足す。`flux` は
    /// [リング] の今のスイープの中性子束（m より前の格子点の値）、`u` は格子の全体のレサジー、`high` は範囲の上端。
    #[allow(clippy::too_many_arguments)]
    pub(super) fn source(
        &mut self,
        m: usize,
        u: &[f64],
        high: f64,
        e: f64,
        flux: &[Vec<f64>],
        known: &mut [f64],
        own: &mut [f64],
    ) -> Result<(), UfgError> {
        if m < self.start {
            return Ok(());
        }
        let clock = Instant::now();
        let last = u.len() - 1;
        // すべての自由ガスの核種の区間を 1 つの並びにする（並列の領域を格子点ごとに 1 つにし、核種ごとに分けた小さな仕事を起こす
        // 手間を避ける）
        self.all.clear();
        for (l, gas) in self.gases.iter().enumerate() {
            // 分離共鳴を持つ核種の 400 kT のように上限が範囲の上端より低い核種で、区間を作って捨てる手間を省く
            if gas.idle(high, u[m], e) {
                continue;
            }
            gas.cells(
                high,
                u[m],
                e,
                &self.grades,
                self.cells,
                u[last],
                &mut self.bounds,
                &mut self.pieces,
            );
            self.all.extend(self.pieces.iter().map(|&(a, b)| (l, a, b)));
        }
        let n = known.len();
        let (gases, sigmas, rings) = (&self.gases, &self.sigma, &self.rings);
        let prefix = &self.prefix;
        let previous = self.previous.as_ref();
        // 区間の塊ごとに並列に積分し、塊の順に足す（並列の分け方によらず同じ和にする）
        let parts: Vec<Vec<(f64, f64)>> = self
            .all
            .par_chunks(CELLS_PER_TASK)
            .map(|chunk| {
                let mut sums = vec![(0.0, 0.0); n];
                // 隣の区間と共有する端の値（核種、端の位置、窓の中として扱ったか、格子点、値）
                let mut shared: Option<(usize, f64, bool, usize, f64)> = None;
                // 前の区間の上端の格子点（[`locate_from`] の起点。塊の最初の区間は格子の全体を探す）
                let mut hint = usize::MAX;
                for &(l, a, b) in chunk {
                    let (gas, (offset, sigma)) = (&gases[l], &sigmas[l]);
                    let inside = gas.inside(u[m], a, b);
                    let f = |x: f64, q: usize| {
                        interpolate(&u[*offset..], sigma, q - offset, x)
                            * gas.correction(high, x, u[m], e, inside)
                    };
                    let (qa, fa) = match shared {
                        Some((k, x, side, q, value)) if k == l && x == a && side == inside => {
                            (q, value)
                        }
                        _ => {
                            let q = locate_from(u, a, hint);
                            (q, f(a, q))
                        }
                    };
                    let mid = 0.5 * (a + b);
                    let q_mid = locate_from(u, mid, qa);
                    let qb = locate_from(u, b, q_mid);
                    hint = qb;
                    let fb = f(b, qb);
                    shared = Some((l, b, inside, qb, fb));
                    let quadratic = Quadratic::through(a, b, fa, f(mid, q_mid), fb);
                    for &(i, density) in &rings[l] {
                        let history = History {
                            flux: &flux[i],
                            prefix: &prefix[i],
                            previous: previous.map(|(f, p)| (f[i].as_slice(), &p[i])),
                        };
                        let (k, o) = integrate(&quadratic, u, m, a, b, qa, qb, &history);
                        sums[i].0 += density * k;
                        sums[i].1 += density * o;
                    }
                }
                sums
            })
            .collect();
        for part in parts {
            for ((k, o), (known, own)) in part.into_iter().zip(known.iter_mut().zip(own.iter_mut()))
            {
                *known += k;
                *own += o;
            }
        }
        let gas_done = Instant::now();
        // S(α, β) の上限より下の入射からの、S(α, β) の移行による散乱源
        for ((table, (offset, sigma)), rings) in self.tables.iter().zip(sigmas).zip(rings) {
            if let Some(table) = *table {
                let source = SabSource {
                    table,
                    u,
                    high,
                    offset: *offset,
                    sigma,
                };
                let previous = previous.map(|(f, _)| f.as_slice());
                source.add(m, e, rings, flux, previous, known, own, &mut self.sab)?;
            }
        }
        self.elapsed[0] += gas_done - clock;
        self.elapsed[1] += gas_done.elapsed();
        Ok(())
    }

    /// 散乱源の経過時間の和 (s): [自由ガスの核の区間, S(α, β) の移行]
    pub(super) fn elapsed(&self) -> [f64; 2] {
        self.elapsed.map(|t| t.as_secs_f64())
    }

    /// 格子点 m の中性子束が決まったあと、自由ガスの核種を含むリングの累積の積率を進める
    pub(super) fn record(&mut self, m: usize, u: &[f64], flux: &[Vec<f64>]) {
        if m < self.start {
            return;
        }
        for rings in &self.rings {
            for &(i, _) in rings {
                // 同じリングを複数の核種が持つときは 1 度だけ進める
                if self.prefix[i].values.len() == m - self.start {
                    self.prefix[i].push(u, &flux[i], m);
                }
            }
        }
    }

    /// スイープを終え、前のスイープとの中性子束の差（start 以降の格子点の相対差の最大。最初のスイープは `None`）を返し、今のスイープを
    /// 次のスイープの「前のスイープ」にする。
    pub(super) fn finish(&mut self, flux: &[Vec<f64>]) -> Option<f64> {
        let change = self.previous.as_ref().map(|(old, _)| {
            flux.iter()
                .zip(old)
                .flat_map(|(f, o)| {
                    f[self.start..]
                        .iter()
                        .zip(&o[self.start..])
                        .map(|(a, b)| ((a - b) / a).abs())
                })
                .fold(0.0, f64::max)
        });
        // 今のスイープの累積の積率（`record` が自由ガスの核種を含むリングで進めたもの）を、次のスイープの「前のスイープ」にする
        let fresh = vec![Prefix::new(self.start); self.prefix.len()];
        let prefixes = std::mem::replace(&mut self.prefix, fresh);
        self.previous = Some((flux.to_vec(), prefixes));
        change
    }
}

/// S(α, β) の移行の表による散乱源（ADR-0205）。表は入射 E' と移行 t = E − E' の累積の割合 C(E', t) を持ち、表の区間の中の核を
/// 一定とする。出射 E（レサジー u）の散乱源 ∫ σs(u') F(u') R(u'→u) du'（R = E p(E'→E)、p は出射のエネルギーあたりの割合）を、
/// 入射を表の移行の格子点と減速計算の格子点で区切った小区間の和にする。小区間では、核の割合 C(E'_c, t_b) − C(E'_c, t_a)
/// （E'_c は小区間の中点の入射）に E/E'_c を掛けたものが ∫ R du' で、σs と F は中点の値（F は格子点の値の線形補間）とする。
/// 核は E' = E に対数の山（準弾性散乱）を持つので、点の値で求積せず、表の累積の差で割合を求める。
pub(super) struct SabSource<'a> {
    pub(super) table: &'a SabTransfer,
    /// 格子の全体のレサジー（昇順）と範囲の上端 (eV)
    pub(super) u: &'a [f64],
    pub(super) high: f64,
    /// 弾性散乱断面積の列の最初の格子点と、そこからの [格子点] の値
    pub(super) offset: usize,
    pub(super) sigma: &'a [f64],
}

/// [`SabSource::add`] の作業用の配列（格子点ごとに再利用する）
#[derive(Debug, Clone, Default)]
pub(super) struct SabWork {
    /// 小区間の移行の両端
    cells: Vec<(f64, f64)>,
    /// 小区間ごとの重み
    terms: Vec<SabTerm>,
}

/// S(α, β) の移行の小区間の（中性子束を補間する格子の区間の始まり、∫ σs R du'、区間の中の位置）。割合が 0 以下の小区間は `None`、
/// 表の外は `Err`（移行の下端）
type SabTerm = Result<Option<(usize, f64, f64)>, f64>;

/// 並列の 1 つの仕事にする小区間の数の下限（格子点ごとの小区間は、表の移行の格子点の数 2389 と同じ桁）
const SAB_CELLS_PER_TASK: usize = 256;

impl SabSource<'_> {
    /// 格子点 m（出射 `e`）の散乱源を、[`Thermal::source`] と同じく、リング（`rings` は（リング, 個数密度））ごとに既知の部分を
    /// `known` に、F_m の係数を `own` に足す。`flux` は今のスイープの [リング] の中性子束（m より前）、`previous` は前のスイープの
    /// 中性子束（最初のスイープは `None` で、m より後の中性子束を F_m とする）。格子の最後の点より後（範囲の下限より下の入射）の
    /// 中性子束は最後の点の値。`work` は作業用。
    ///
    /// 小区間の区切りを順に決め、小区間ごとの重みを並列に求めてから、小区間の順に足す（足す順が直列と同じなので、和は並列の分け方に
    /// よらない）
    #[allow(clippy::too_many_arguments)]
    pub(super) fn add(
        &self,
        m: usize,
        e: f64,
        rings: &[(usize, f64)],
        flux: &[Vec<f64>],
        previous: Option<&[Vec<f64>]>,
        known: &mut [f64],
        own: &mut [f64],
        work: &mut SabWork,
    ) -> Result<(), UfgError> {
        let (table, u, high) = (self.table, self.u, self.high);
        let last = u.len() - 1;
        let floor = table.incident[table.incident.len() - 1];
        // 入射は表の上限から（t = e − 上限）、表の入射の最初の点か移行の最大まで
        let t_low = e - floor;
        let t_high = table.transfer[table.transfer.len() - 1].min(e - table.incident[0]);
        if t_low >= t_high {
            return Ok(());
        }
        // 格子点 q の移行 t_q = e − E'_q（q とともに増える）。表の上限以下の入射の最初の格子点から
        let transfer_at = |q: usize| e - high * (-u[q]).exp();
        let mut q = u.partition_point(|&x| high * (-x).exp() > floor);
        let mut j = table.transfer.partition_point(|&t| t <= t_low);
        let cells = &mut work.cells;
        cells.clear();
        let mut a = t_low;
        while a < t_high {
            // 次の区切り: 表の移行の格子点か、減速計算の格子点の移行の小さい方
            while q <= last && transfer_at(q) <= a {
                q += 1;
            }
            let next_grid = if q <= last {
                transfer_at(q)
            } else {
                f64::INFINITY
            };
            let next_table = table.transfer.get(j).copied().unwrap_or(f64::INFINITY);
            let b = next_grid.min(next_table).min(t_high);
            if next_table <= b {
                j += 1;
            }
            cells.push((a, b));
            // 区切りは a より先にあり（格子点の q と表の j は、a 以下の点を飛ばしてある）、小区間は空にならない。空になれば、その
            // 小区間を誤りとして残し（前の小区間の誤りが先に報告されるように）、ここで止める
            if b <= a {
                break;
            }
            a = b;
        }
        // 小区間の核の割合、σs、中性子束の補間の位置は、小区間の両端だけで決まる
        let term = |hint: &mut usize, &(a, b): &(f64, f64)| {
            if b <= a {
                return Err(a);
            }
            let center = 0.5 * (a + b);
            let incident = e - center;
            let x = (high / incident).ln();
            // 入射と移行は表の中にある（入射は表の最初の点から上限まで、移行は上限から範囲の下限への下方散乱から表の最大まで。
            // 表がそれを覆うことは、核種のデータを引くときに確かめてある）
            let row = table.row(incident).ok_or(a)?;
            let (Some(cb), Some(ca)) = (row.cumulative(b), row.cumulative(a)) else {
                return Err(a);
            };
            let share = cb - ca;
            if share <= 0.0 {
                return Ok(None);
            }
            // 中点を含む格子の区間の始まり（u_cell ≤ x の最後の格子点）。x は小区間とともに増える
            let cell = locate_from(u, x, *hint);
            *hint = cell;
            let sigma = interpolate(
                &u[self.offset..],
                self.sigma,
                cell.max(self.offset) - self.offset,
                x,
            );
            let weight = sigma * share * e / incident;
            // F は格子の区間 [u_cell, u_cell+1] の線形補間（最後の点より後は最後の点の値）
            Ok(Some(if cell >= last || x >= u[last] {
                (last, weight, 0.0)
            } else {
                (cell, weight, (x - u[cell]) / (u[cell + 1] - u[cell]))
            }))
        };
        cells
            .par_iter()
            .with_min_len(SAB_CELLS_PER_TASK)
            .map_init(|| usize::MAX, term)
            .collect_into_vec(&mut work.terms);
        let assign = |idx: usize, c: f64, i: usize, known: &mut [f64], own: &mut [f64]| {
            if idx < m {
                known[i] += c * flux[i][idx];
            } else if idx == m {
                own[i] += c;
            } else {
                match previous {
                    Some(p) => known[i] += c * p[i][idx],
                    None => own[i] += c,
                }
            }
        };
        for term in &work.terms {
            let Some((lower, weight, tau)) = term.map_err(|transfer| UfgError::SabTable {
                energy: e,
                transfer,
            })?
            else {
                continue;
            };
            for &(i, density) in rings {
                let c = density * weight;
                assign(lower, c * (1.0 - tau), i, known, own);
                if tau > 0.0 {
                    assign(lower + 1, c * tau, i, known, own);
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
