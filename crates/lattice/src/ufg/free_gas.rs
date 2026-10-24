//! UFG の減速計算の、分離共鳴を持たない核種（軽い核種）の標的の熱運動を、温度 T の自由ガスで扱う散乱核（ADR-0203）。
//!
//! 入射 E'（レサジー u'）、出射 E（u）のレサジーあたりの散乱核は R(u'→u) = E S(E'→E)/N(E')。S は断面積を一定とした
//! 自由ガスの核（[`kernel`]）、N は温度に拡げた断面積と自由原子の断面積の比（[`effective_ratio`]）で、ライブラリの
//! ポイントワイズの弾性散乱断面積（温度に拡げたもの）に R を掛けると、自由原子の断面積 × S になる。範囲の中で
//!
//!   R = R_static/N(E') + D、R_static = e^{−(u−u')}/(1 − α)（u − Δ ≤ u' ≤ u）
//!
//! と分けると、R_static/N は静止した標的の窓の積分（弾性散乱断面積を σs/N にする）で扱え、D は静止した核の両端の近くに
//! だけ残る: 上端 E' = E の近く（上方散乱の e^{−(E−E')/kT} の裾と、下方散乱の側の幅 4kT/A ほどの折れ）と、下端 E' = E/α の
//! 近く（erf の引数 θ√(E/kT) − ζ√(E'/kT) が 0 を横切るぼけ。E' の幅 2√(E' kT)/ζ）。D は、その近くの区間（[`Light::cells`]）
//! ごとに 3 点の 2 次式で近似し、中性子束の累積の積率（[`Prefix`]）で積分する。区間の外の D は、上端の側で e^{−[`TAIL`]}、
//! 下端の側で erfc([`EDGE`]) 程度より小さい（`correction_vanishes_outside_its_cells`）。

use std::f64::consts::PI;

use xlc_numeric::quadrature::GAUSS_LEGENDRE_2;

/// ボルツマン定数 (eV/K)。CODATA 2018 の値（xs の `physics::BK` と同じ）
pub(super) const BOLTZMANN: f64 = 8.617_333_262e-5;

/// 上端の側の区間の幅（kT の倍数）。上方散乱の裾 e^{−(E−E')/kT} がここで e^{−28} ≈ 7e−13 に下がる
pub(super) const TAIL: f64 = 28.0;

/// 下端の側の区間の半幅（ぼけの幅 2√(E' kT)/ζ の倍数）。erfc(6) ≈ 2e−17
pub(super) const EDGE: f64 = 6.0;

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
/// 連分数の 40 項を省く
const SATURATION: f64 = 6.0;

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
    2.0 / PI.sqrt() * (-x2).exp() * sum
}

/// スケールした相補誤差関数 erfcx(x) = exp(x²) erfc(x)（x ≥ 0）。xs の `special::erfcx` と同じ計算。
pub(super) fn erfcx(x: f64) -> f64 {
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
    1.0 / (PI.sqrt() * t)
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
    // 求めない（erfcx(x) ≤ 1 なので積は 0）
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
            let pair = |y: f64, a: f64, b: f64| {
                let factor = (-y * y).exp();
                if factor == 0.0 {
                    0.0
                } else {
                    factor * (erfcx(a) - erfcx(b))
                }
            };
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
            let pair = |y: f64, a: f64, b: f64| {
                let factor = (-y * y).exp();
                if factor == 0.0 {
                    0.0
                } else {
                    factor * (erfcx(a) - erfcx(b))
                }
            };
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

/// 自由ガスの核で散乱する 1 つの軽い核種（名前と温度の組）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Light {
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
}

impl Light {
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
        }
    }

    /// 自由ガスの核を使う入射エネルギーの上限 (eV)
    pub(super) fn cap(&self) -> f64 {
        self.cap
    }

    /// kT (eV)
    pub(super) fn kt(&self) -> f64 {
        self.kt
    }

    /// 入射エネルギー `e_in` の静止した標的の窓に掛ける 1/N(E')（自由ガスの核を使わない入射エネルギーでは 1）
    pub(super) fn window_scale(&self, e_in: f64) -> f64 {
        if e_in < self.cap {
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
    /// 1/`cells` の等間隔で ±[`EDGE`] 幅を覆う。両方を合わせ、u、u − Δ、自由ガスの核の上限、範囲の上端（u' = 0）、`bottom`
    /// （格子の最後の点）で区切る（D と中性子束の扱いが、区間の中で変わらないように）。範囲の上（u' < 0）と、自由ガスの
    /// 核を使わない入射エネルギーの区間は落とす（静止した標的で扱う）。
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
        let reach = TAIL.min(0.5 / scale);
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
        let breaks = [center, (high / self.cap).ln(), 0.0, bottom];
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
            if b > a && mid >= 0.0 && high * (-mid).exp() < self.cap && covered(mid) {
                out.push((a, b));
            }
        }
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

    /// 格子点 `start`, `start` + 1, … の中性子束 `flux`（全体の添字）の累積の積率
    pub(super) fn of(u: &[f64], flux: &[f64], start: usize) -> Self {
        let mut prefix = Self::new(start);
        for m in start..flux.len() {
            prefix.push(u, flux, m);
        }
        prefix
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

    /// 戻す: 格子点 `len` 以降を捨てる
    pub(super) fn truncate(&mut self, len: usize) {
        self.values.truncate(len.saturating_sub(self.start));
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

    fn at(&self, s: f64) -> f64 {
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

/// 減速計算のスイープの間の、自由ガスの核の D の散乱源の状態。
#[derive(Debug, Clone)]
pub(super) struct Thermal {
    /// [軽い核種] の核
    lights: Vec<Light>,
    /// [軽い核種][格子点] のポイントワイズの弾性散乱断面積（D に掛ける）
    sigma: Vec<Vec<f64>>,
    /// [軽い核種] の（リング, 個数密度）
    rings: Vec<Vec<(usize, f64)>>,
    /// 上端の側の区間の段（[`grades`]）と、下端の側の区間の細かさ
    grades: Vec<f64>,
    cells: usize,
    /// D が 0 でありうる最初の格子点の 1 つ前（ここから掃き直す）
    start: usize,
    /// [リング] の今のスイープの累積の積率
    prefix: Vec<Prefix>,
    /// 前のスイープの [リング] の中性子束と累積の積率
    previous: Option<(Vec<Vec<f64>>, Vec<Prefix>)>,
    /// 作業用の区間の境界と区間
    bounds: Vec<f64>,
    pieces: Vec<(f64, f64)>,
}

impl Thermal {
    /// 軽い核種の核 `lights`、その [格子点] の弾性散乱断面積 `sigma`、（リング, 個数密度）`rings`、区間の細かさ `cells`
    /// の状態。`energies` は格子点のエネルギー（降順）、`n` はリングの数。軽い核種がなければ `None`。
    pub(super) fn new(
        lights: Vec<Light>,
        sigma: Vec<Vec<f64>>,
        rings: Vec<Vec<(usize, f64)>>,
        cells: usize,
        energies: &[f64],
        n: usize,
    ) -> Option<Self> {
        if lights.is_empty() {
            return None;
        }
        // 上限より上の入射エネルギーは静止した標的なので、上方散乱の裾が上限の下から届く格子点より前は D が 0
        let reach = lights
            .iter()
            .map(|l| l.cap() + TAIL * l.kt())
            .fold(0.0, f64::max);
        let first = energies.iter().position(|&e| e <= reach)?;
        let start = first.saturating_sub(1);
        Some(Self {
            lights,
            sigma,
            rings,
            grades: grades(cells),
            cells,
            start,
            prefix: vec![Prefix::new(start); n],
            previous: None,
            bounds: Vec::new(),
            pieces: Vec::new(),
        })
    }

    pub(super) fn start(&self) -> usize {
        self.start
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
    ) {
        if m < self.start {
            return;
        }
        let last = u.len() - 1;
        for ((light, sigma), rings) in self.lights.iter().zip(&self.sigma).zip(&self.rings) {
            light.cells(
                high,
                u[m],
                e,
                &self.grades,
                self.cells,
                u[last],
                &mut self.bounds,
                &mut self.pieces,
            );
            // 隣の区間と共有する端の値（端の位置、窓の中として扱ったか、格子点、値）
            let mut shared: Option<(f64, bool, usize, f64)> = None;
            for &(a, b) in &self.pieces {
                let inside = light.inside(u[m], a, b);
                let f = |x: f64, q: usize| {
                    interpolate(u, sigma, q, x) * light.correction(high, x, u[m], e, inside)
                };
                let (qa, fa) = match shared {
                    Some((x, side, q, value)) if x == a && side == inside => (q, value),
                    _ => {
                        let q = locate(u, a);
                        (q, f(a, q))
                    }
                };
                let qb = locate(u, b);
                let fb = f(b, qb);
                shared = Some((b, inside, qb, fb));
                let mid = 0.5 * (a + b);
                let quadratic = Quadratic::through(a, b, fa, f(mid, locate(u, mid)), fb);
                for &(i, density) in rings {
                    let history = History {
                        flux: &flux[i],
                        prefix: &self.prefix[i],
                        previous: self
                            .previous
                            .as_ref()
                            .map(|(f, p)| (f[i].as_slice(), &p[i])),
                    };
                    let (k, o) = integrate(&quadratic, u, m, a, b, qa, qb, &history);

                    known[i] += density * k;
                    own[i] += density * o;
                }
            }
        }
    }

    /// 格子点 m の中性子束が決まったあと、軽い核種を含むリングの累積の積率を進める
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
    pub(super) fn finish(&mut self, u: &[f64], flux: &[Vec<f64>]) -> Option<f64> {
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
        // 累積の積率は、軽い核種を含むリングだけ
        let prefixes = flux
            .iter()
            .enumerate()
            .map(|(i, f)| {
                if self.rings.iter().flatten().any(|&(r, _)| r == i) {
                    Prefix::of(u, f, self.start)
                } else {
                    Prefix::new(self.start)
                }
            })
            .collect();
        self.previous = Some((flux.to_vec(), prefixes));
        for p in &mut self.prefix {
            p.truncate(self.start);
        }
        change
    }
}

#[cfg(test)]
mod tests;
