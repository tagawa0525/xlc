//! 均質な無限媒質の B1 の基本モード（臨界スペクトル）と、そこから求める拡散係数（ADR-0104）。
//!
//! 中性子束を平面波 φ(r) = φ e^{iB·r} とした多群の B1 の方程式（実数の形）:
//!
//! - Σt,g φ_g − Σ_g' Σs0(g'→g) φ_g' + B² d_g = χ_g Σ_g' νΣf,g' φ_g' / k
//! - 3 γ(B²/Σt,g²) Σt,g d_g − 3 Σ_g' Σs1(g'→g) d_g' = φ_g
//!
//! d_g = J_g / B は中性子流を B で割ったもので、拡散係数は D_g = J_g / (B φ_g) = d_g / φ_g になる。
//! γ(x²) = x² arctan(x) / (3 (x − arctan x))（x = B/Σt）は、P1 を超える高次の角度の効果を入れる
//! 関数である（B → 0 で 1 になり、P1 の方程式に戻る）。B² < 0（x = iy）では解析接続
//! γ = y² artanh(y) / (3 (artanh(y) − y)) を使い、|B| < Σt,g のすべての群で定まる。
//! 出典: R. J. J. Stamm'ler, M. J. Abbate, "Methods of Steady-State Reactor Physics in Nuclear
//! Design" (1983) の B1 の方程式、A. Hébert, "Applied Reactor Physics" (2009) の均質な B1 の漏れの
//! モデル（同じ形）。Σs1 は散乱角の余弦の Legendre の1次のモーメント（(2l + 1) を含まない）で、
//! ライブラリの散乱行列（ADR-0076）と同じ定義。

use xlc_numeric::lu::Lu;

/// γ を級数で計算する |B²/Σt²| の上限（これより小さいと閉じた式は桁落ちする）
const SERIES_LIMIT: f64 = 1.0e-3;
/// 臨界のバックリングの探索の、k − 1 の許容値
const K_TOLERANCE: f64 = 1.0e-12;
/// 探索が区間の幅で止まったときに、臨界とみなす k − 1 の上限
const CRITICAL_TOLERANCE: f64 = 1.0e-9;
/// 臨界のバックリングの探索の反復の上限
const MAX_ITERATIONS: usize = 200;
/// 正のバックリングの上限を探すときに倍にする回数の上限
const MAX_DOUBLINGS: usize = 200;
/// 負のバックリングの下限の、−min(Σt)² に対する割合（−min(Σt)² では γ の解析接続が端に達する）
const NEGATIVE_LIMIT: f64 = 1.0 - 1.0e-12;

/// 細かい群のまま均質化した無限媒質の断面積 (1/cm)。
#[derive(Debug, Clone, PartialEq)]
pub struct Medium {
    pub total: Vec<f64>,
    /// P0 の散乱行列 [入射の群][出射の群]
    pub scatter: Vec<Vec<f64>>,
    /// P1 の散乱行列 [入射の群][出射の群]
    pub scatter1: Vec<Vec<f64>>,
    pub nu_fission: Vec<f64>,
    /// 核分裂スペクトル（和が 1。核分裂のない媒質では 0）
    pub chi: Vec<f64>,
}

/// あるバックリングでの B1 の基本モード。
#[derive(Debug, Clone, PartialEq)]
pub struct Mode {
    /// B² (1/cm²)。臨界のモードでは、k∞ > 1 なら正、k∞ < 1 なら負（B = iκ で、中性子束は空間で
    /// 指数関数的に変わり、漏れは負になる）
    pub buckling: f64,
    /// 漏れを入れた増倍率
    pub k: f64,
    /// [群] の中性子束。単位の χ の源（Σ χ_g = 1）で解いたもので、νΣf·φ = k
    pub flux: Vec<f64>,
    /// [群] の中性子流を B で割った d_g = J_g / B（中性子束と同じ規格化）
    pub current: Vec<f64>,
}

impl Mode {
    /// [群] の拡散係数 D_g = J_g / (B φ_g)。
    pub fn diffusion(&self) -> Vec<f64> {
        self.current
            .iter()
            .zip(&self.flux)
            .map(|(d, f)| d / f)
            .collect()
    }

    /// 臨界の探索で、物理的な基本モードとして使えるか（k と中性子束と中性子流が有限、k > 0、中性子束が負でない）。
    /// B² < 0 で極を越えると、中性子束の符号が変わるか k が負になる
    fn physical(&self) -> bool {
        let finite = |v: &[f64]| v.iter().all(|x| x.is_finite());
        self.k.is_finite()
            && self.k > 0.0
            && finite(&self.flux)
            && finite(&self.current)
            && self.flux.iter().all(|&f| f >= 0.0)
    }
}

/// B1 の計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum B1Error {
    #[error("B1 cross section arrays must all have {0} groups")]
    Groups(usize),
    #[error("group {group} has total cross section {total}")]
    Total { group: usize, total: f64 },
    #[error("B1 cross sections must be finite")]
    NotFinite,
    #[error("the medium has no fission source")]
    NoFission,
    #[error("buckling {buckling} is below −Σt² = {limit} of some group")]
    Imaginary { buckling: f64, limit: f64 },
    #[error("the B1 equations are singular at buckling {0}")]
    Singular(f64),
    #[error("the B1 mode at buckling {0} is not physical")]
    NotPhysical(f64),
    #[error("no critical buckling: k is {k} at buckling {buckling}")]
    NoCriticalBuckling { buckling: f64, k: f64 },
    #[error(
        "k∞ < 1 but the mode at buckling {buckling} near −min(Σt)² is still physical with k = {k} < 1"
    )]
    NoNegativeBracket { buckling: f64, k: f64 },
}

/// γ(s)（s = B²/Σt²）。s ≤ −1 では定まらない（呼ぶ側で除く）。
pub(crate) fn gamma(s: f64) -> f64 {
    if s.abs() < SERIES_LIMIT {
        // x² arctan(x) / (3 (x − arctan x)) の s = x² の級数（s < 0 でも同じ式）
        1.0 + s * (4.0 / 15.0 + s * (-12.0 / 175.0 + s * 92.0 / 2625.0))
    } else if s > 0.0 {
        let x = s.sqrt();
        let a = x.atan();
        s * a / (3.0 * (x - a))
    } else {
        let y = (-s).sqrt();
        let a = y.atanh();
        -s * a / (3.0 * (a - y))
    }
}

impl Medium {
    fn groups(&self) -> usize {
        self.total.len()
    }

    fn validate(&self) -> Result<(), B1Error> {
        let ng = self.groups();
        let square = |m: &[Vec<f64>]| m.len() == ng && m.iter().all(|row| row.len() == ng);
        if ng == 0
            || !square(&self.scatter)
            || !square(&self.scatter1)
            || self.nu_fission.len() != ng
            || self.chi.len() != ng
        {
            return Err(B1Error::Groups(ng));
        }
        let finite = |v: &[f64]| v.iter().all(|x| x.is_finite());
        if !finite(&self.total)
            || !finite(&self.nu_fission)
            || !finite(&self.chi)
            || !self.scatter.iter().chain(&self.scatter1).all(|r| finite(r))
        {
            return Err(B1Error::NotFinite);
        }
        if let Some((group, &total)) = self.total.iter().enumerate().find(|(_, t)| **t <= 0.0) {
            return Err(B1Error::Total { group, total });
        }
        Ok(())
    }

    fn multiplying(&self) -> bool {
        self.nu_fission.iter().any(|&v| v > 0.0) && self.chi.iter().any(|&c| c > 0.0)
    }

    /// B² の下限 −min(Σt)²（これ以下では γ が定まらない）
    fn negative_limit(&self) -> f64 {
        let min = self.total.iter().copied().fold(f64::INFINITY, f64::min);
        -min * min
    }

    /// バックリング B² での基本モード。核分裂源は χ の形に決まっているので、
    /// [φ; d] を χ の源で1回解き、k = νΣf·φ とする。
    pub fn mode(&self, buckling: f64) -> Result<Mode, B1Error> {
        self.validate()?;
        if !self.multiplying() {
            return Err(B1Error::NoFission);
        }
        let limit = self.negative_limit();
        if buckling.is_nan() || buckling <= limit || buckling == f64::INFINITY {
            return Err(B1Error::Imaginary { buckling, limit });
        }
        self.solve(buckling).ok_or(B1Error::Singular(buckling))
    }

    /// 検査済みの媒質の、B² でのモード。行列が特異なら `None`。
    fn solve(&self, buckling: f64) -> Option<Mode> {
        let ng = self.groups();
        let n = 2 * ng;
        // 行: 0..ng が中性子収支、ng..2ng が中性子流の式。列: 0..ng が φ、ng..2ng が d
        let mut a = vec![0.0; n * n];
        for g in 0..ng {
            let total = self.total[g];
            a[g * n + g] += total;
            a[g * n + ng + g] = buckling;
            a[(ng + g) * n + g] = -1.0;
            a[(ng + g) * n + ng + g] += 3.0 * gamma(buckling / (total * total)) * total;
            for from in 0..ng {
                a[g * n + from] -= self.scatter[from][g];
                a[(ng + g) * n + ng + from] -= 3.0 * self.scatter1[from][g];
            }
        }
        let lu = Lu::new(a, n)?;
        let mut x = vec![0.0; n];
        x[..ng].copy_from_slice(&self.chi);
        lu.solve(&mut x);
        let current = x.split_off(ng);
        let k = self.nu_fission.iter().zip(&x).map(|(v, f)| v * f).sum();
        Some(Mode {
            buckling,
            k,
            flux: x,
            current,
        })
    }

    /// B² = 0 の中性子流の式（P1 の式、γ = 1）3 Σt,g d_g − 3 Σ_g' Σs1(g'→g) d_g' = φ_g を、与えた
    /// 中性子束 φ で解いた d（= D φ）。核分裂のない媒質では臨界のバックリングが定まらないので、これを使う。
    pub fn zero_buckling_current(&self, flux: &[f64]) -> Result<Vec<f64>, B1Error> {
        self.validate()?;
        let ng = self.groups();
        if flux.len() != ng {
            return Err(B1Error::Groups(ng));
        }
        let mut a = vec![0.0; ng * ng];
        for g in 0..ng {
            a[g * ng + g] += 3.0 * self.total[g];
            for from in 0..ng {
                a[g * ng + from] -= 3.0 * self.scatter1[from][g];
            }
        }
        let lu = Lu::new(a, ng).ok_or(B1Error::Singular(0.0))?;
        let mut d = flux.to_vec();
        lu.solve(&mut d);
        Ok(d)
    }

    /// k = 1 になる臨界のバックリングを探す。物理的なモードの k は B² について単調に減るので、k(0) > 1 なら
    /// B² > 0 にある。k(0) < 1 なら B² < 0 にあり、B² を下げると漏れ（負）が除去を打ち消して k は極に
    /// 向かって増え、極を越えると物理的でなくなる（`Mode::physical`）。区間 (−min(Σt)², 0) の下端が
    /// 物理的で k < 1 なら、この前提が成り立たないのでエラーにする。
    ///
    /// 区間を保つ Illinois 法（はさみうち法の変形）で k − 1 = 0 を解く。物理的でない点（特異な行列、
    /// 有限でない k を含む）は「負に過ぎる」として下端を動かし、そのときは二分する。
    pub fn critical(&self) -> Result<Mode, B1Error> {
        self.validate()?;
        if !self.multiplying() {
            return Err(B1Error::NoFission);
        }
        let mut evaluations = 0;
        let mut trial = |buckling: f64| {
            evaluations += 1;
            self.solve(buckling).filter(Mode::physical)
        };
        let zero = trial(0.0).ok_or(B1Error::NotPhysical(0.0))?;
        if (zero.k - 1.0).abs() <= K_TOLERANCE {
            return Ok(zero);
        }
        // 下端（k > 1 か物理的でない）と上端（k < 1）。f = k − 1 は物理的な点だけで持つ
        let (mut low, mut f_low, mut high, mut f_high) = if zero.k > 1.0 {
            let mut high = self
                .total
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min)
                .powi(2);
            let mut doublings = 0;
            let f_high = loop {
                let mode = trial(high).ok_or(B1Error::NotPhysical(high))?;
                if mode.k < 1.0 {
                    break mode.k - 1.0;
                }
                doublings += 1;
                if doublings > MAX_DOUBLINGS {
                    return Err(B1Error::NoCriticalBuckling {
                        buckling: high,
                        k: mode.k,
                    });
                }
                high *= 2.0;
            };
            (0.0, Some(zero.k - 1.0), high, f_high)
        } else {
            let low = NEGATIVE_LIMIT * self.negative_limit();
            let f_low = match trial(low) {
                Some(mode) if mode.k < 1.0 => {
                    return Err(B1Error::NoNegativeBracket {
                        buckling: low,
                        k: mode.k,
                    });
                }
                Some(mode) => Some(mode.k - 1.0),
                None => None,
            };
            (low, f_low, 0.0, zero.k - 1.0)
        };
        let mut best = zero;
        // Illinois 法: 同じ側の端が2回続けて残ったら、その端の f を半分にする
        let mut kept: Option<bool> = None;
        for _ in 0..MAX_ITERATIONS {
            let middle = 0.5 * (low + high);
            let x = match f_low {
                Some(fl) => {
                    let x = high - f_high * (high - low) / (f_high - fl);
                    if x > low && x < high { x } else { middle }
                }
                None => middle,
            };
            let mode = trial(x);
            let too_negative = mode.as_ref().is_none_or(|m| m.k > 1.0);
            if too_negative {
                low = x;
                f_low = mode.as_ref().map(|m| m.k - 1.0);
                if kept == Some(true) {
                    f_high *= 0.5;
                }
                kept = Some(true);
            } else if let Some(m) = &mode {
                high = x;
                f_high = m.k - 1.0;
                if kept == Some(false)
                    && let Some(fl) = f_low.as_mut()
                {
                    *fl *= 0.5;
                }
                kept = Some(false);
            }
            if let Some(m) = mode
                && (m.k - 1.0).abs() < (best.k - 1.0).abs()
            {
                best = m;
            }
            if (best.k - 1.0).abs() <= K_TOLERANCE
                || high - low <= f64::EPSILON * low.abs().max(high.abs())
            {
                break;
            }
        }
        if (best.k - 1.0).abs() > CRITICAL_TOLERANCE {
            return Err(B1Error::NoCriticalBuckling {
                buckling: best.buckling,
                k: best.k,
            });
        }
        tracing::info!(
            buckling = best.buckling,
            evaluations,
            "B1 critical spectrum"
        );
        Ok(best)
    }
}

#[cfg(test)]
mod tests;
