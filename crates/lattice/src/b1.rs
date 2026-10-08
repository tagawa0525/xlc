//! 均質な無限媒質の B1 の基本モード（臨界スペクトル）と、そこから求める拡散係数（ADR-0098）。
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

use crate::lu::Lu;

/// γ を級数で計算する |B²/Σt²| の上限（これより小さいと閉じた式は桁落ちする）
const SERIES_LIMIT: f64 = 1.0e-3;
/// 臨界のバックリングの探索の、k − 1 の許容値
const K_TOLERANCE: f64 = 1.0e-12;
/// 二分法が区間の幅で止まったときに、臨界とみなす k − 1 の上限
const CRITICAL_TOLERANCE: f64 = 1.0e-9;
/// 二分法の反復の上限
const MAX_BISECTIONS: usize = 200;
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
    /// 核分裂スペクトル（和が 1）
    pub chi: Vec<f64>,
}

/// あるバックリングでの B1 の基本モード。
#[derive(Debug, Clone, PartialEq)]
pub struct Mode {
    /// B² (1/cm²)。負なら体系は未臨界で、中性子流は中性子束の勾配を上る
    pub buckling: f64,
    /// 漏れを入れた増倍率
    pub k: f64,
    /// [群] の中性子束（νΣf·φ = k に規格化）
    pub flux: Vec<f64>,
    /// [群] の中性子流を B で割った d_g = J_g / B
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
}

/// B1 の計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum B1Error {
    #[error("B1 cross section arrays must all have {0} groups")]
    Groups(usize),
    #[error("group {group} has total cross section {total}")]
    Total { group: usize, total: f64 },
    #[error("the medium has no fission source")]
    NoFission,
    #[error("buckling {buckling} is below −Σt² = {limit} of some group")]
    Imaginary { buckling: f64, limit: f64 },
    #[error("the B1 equations are singular at buckling {0}")]
    Singular(f64),
    #[error("no critical buckling: k is {k} at buckling {buckling}")]
    NoCriticalBuckling { buckling: f64, k: f64 },
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
        if let Some((group, &total)) = self
            .total
            .iter()
            .enumerate()
            .find(|(_, t)| !(**t > 0.0 && t.is_finite()))
        {
            return Err(B1Error::Total { group, total });
        }
        if self.nu_fission.iter().all(|&v| v <= 0.0) || self.chi.iter().all(|&c| c <= 0.0) {
            return Err(B1Error::NoFission);
        }
        Ok(())
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
        let limit = self.negative_limit();
        if buckling <= limit || !buckling.is_finite() {
            return Err(B1Error::Imaginary { buckling, limit });
        }
        let ng = self.groups();
        let n = 2 * ng;
        // 行: 0..ng が中性子の釣り合い、ng..2ng が中性子流の式。列: 0..ng が φ、ng..2ng が d
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
        let lu = Lu::new(a, n).ok_or(B1Error::Singular(buckling))?;
        let mut x = vec![0.0; n];
        x[..ng].copy_from_slice(&self.chi);
        lu.solve(&mut x);
        let current = x.split_off(ng);
        let k = self.nu_fission.iter().zip(&x).map(|(v, f)| v * f).sum();
        Ok(Mode {
            buckling,
            k,
            flux: x,
            current,
        })
    }

    /// k = 1 になる臨界のバックリングを二分法で探す。物理的な（中性子束がすべて正の）モードの k は B² に
    /// ついて単調に減るので、k(0) > 1 なら B² > 0 にある。k(0) < 1 なら B² < 0 にあり、B² を下げると
    /// 漏れ（負）が除去を打ち消して k は極に向かって増え、極を越えると中性子束が正でなくなる。
    /// 極は −min(Σt)² より上にあるので、(−min(Σt)², 0) を、k > 1 か物理的でないモードを「負に過ぎる」
    /// として二分する。
    pub fn critical(&self) -> Result<Mode, B1Error> {
        let zero = self.mode(0.0)?;
        if (zero.k - 1.0).abs() <= K_TOLERANCE {
            return Ok(zero);
        }
        let (mut low, mut high) = if zero.k > 1.0 {
            let mut high = self
                .total
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min)
                .powi(2);
            let mut doublings = 0;
            loop {
                let mode = self.mode(high)?;
                if mode.k < 1.0 {
                    break;
                }
                doublings += 1;
                if doublings > MAX_DOUBLINGS {
                    return Err(B1Error::NoCriticalBuckling {
                        buckling: high,
                        k: mode.k,
                    });
                }
                high *= 2.0;
            }
            (0.0, high)
        } else {
            (NEGATIVE_LIMIT * self.negative_limit(), 0.0)
        };
        let mut best = zero;
        for _ in 0..MAX_BISECTIONS {
            let middle = 0.5 * (low + high);
            let mode = self.mode(middle)?;
            let physical = mode.flux.iter().all(|&f| f > 0.0);
            if mode.k > 1.0 || !physical {
                low = middle;
            } else {
                high = middle;
            }
            if physical && (mode.k - 1.0).abs() < (best.k - 1.0).abs() {
                best = mode;
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
        Ok(best)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tolerance: f64) -> bool {
        (a - b).abs() <= tolerance * b.abs()
    }

    fn two_group(scatter1: [[f64; 2]; 2], nu_fission: [f64; 2]) -> Medium {
        Medium {
            total: vec![0.55, 1.4],
            scatter: vec![vec![0.48, 0.03], vec![0.001, 1.25]],
            scatter1: scatter1.iter().map(|r| r.to_vec()).collect(),
            nu_fission: nu_fission.to_vec(),
            chi: vec![1.0, 0.0],
        }
    }

    #[test]
    fn gamma_series_joins_the_closed_form() {
        for s in [SERIES_LIMIT, -SERIES_LIMIT] {
            let series = gamma(s * (1.0 - 1e-9));
            let x = s.abs().sqrt();
            let closed = if s > 0.0 {
                s * x.atan() / (3.0 * (x - x.atan()))
            } else {
                -s * x.atanh() / (3.0 * (x.atanh() - x))
            };
            assert!(close(series, closed, 1e-10), "{s}: {series} vs {closed}");
        }
        assert_eq!(gamma(0.0), 1.0);
    }

    #[test]
    fn isotropic_scattering_gives_one_third_total_at_zero_buckling() {
        let medium = two_group([[0.0; 2]; 2], [0.01, 0.3]);
        let mode = medium.mode(0.0).unwrap();
        for (d, t) in mode.diffusion().iter().zip(&medium.total) {
            assert!(close(*d, 1.0 / (3.0 * t), 1e-13), "{d} vs {t}");
        }
    }

    #[test]
    fn isotropic_scattering_follows_gamma_at_finite_buckling() {
        // Σs1 = 0 なら中性子流の式は群ごとに閉じ、D_g = 1/(3 γ(B²/Σt²) Σt)
        let medium = two_group([[0.0; 2]; 2], [0.01, 0.3]);
        for buckling in [0.05, -0.2] {
            let mode = medium.mode(buckling).unwrap();
            for (d, t) in mode.diffusion().iter().zip(&medium.total) {
                let expected = 1.0 / (3.0 * gamma(buckling / (t * t)) * t);
                assert!(close(*d, expected, 1e-13), "{buckling}: {d} vs {expected}");
            }
        }
    }

    #[test]
    fn diagonal_p1_scattering_gives_the_transport_cross_section() {
        let medium = two_group([[0.2, 0.0], [0.0, 0.6]], [0.01, 0.3]);
        let mode = medium.mode(0.0).unwrap();
        let expected = [1.0 / (3.0 * (0.55 - 0.2)), 1.0 / (3.0 * (1.4 - 0.6))];
        for (d, e) in mode.diffusion().iter().zip(expected) {
            assert!(close(*d, e, 1e-13), "{d} vs {e}");
        }
    }

    /// 1群の等方散乱の平面波の輸送方程式の分散関係から、臨界の x = B/Σt に対する c = (Σs + νΣf)/Σt を作る。
    /// 正の B²: c = x / arctan(x)。負の B²（x = iy）: c = y / artanh(y)。釣り合いから D = Σt (c − 1) / B²
    fn one_group(buckling_over_total2: f64) -> (Medium, f64) {
        let total = 0.8;
        let scatter = 0.5;
        let s: f64 = buckling_over_total2;
        let x = s.abs().sqrt();
        let c = if s > 0.0 { x / x.atan() } else { x / x.atanh() };
        let medium = Medium {
            total: vec![total],
            scatter: vec![vec![scatter]],
            scatter1: vec![vec![0.0]],
            nu_fission: vec![c * total - scatter],
            chi: vec![1.0],
        };
        let buckling = s * total * total;
        (medium, total * (c - 1.0) / buckling)
    }

    #[test]
    fn one_group_critical_mode_matches_the_transport_dispersion_relation() {
        // 1群の等方散乱では B1 は厳密で、臨界のバックリングと D が輸送方程式の解析解に一致する
        for s in [0.25, -0.25, 2.0e-4, -2.0e-4] {
            let (medium, expected) = one_group(s);
            let mode = medium.critical().unwrap();
            let buckling = s * 0.8 * 0.8;
            assert!(
                close(mode.buckling, buckling, 1e-8),
                "{} vs {buckling}",
                mode.buckling
            );
            assert!(close(mode.diffusion()[0], expected, 1e-8));
            let gamma_form = 1.0 / (3.0 * gamma(s) * 0.8);
            assert!(close(mode.diffusion()[0], gamma_form, 1e-8));
        }
    }

    #[test]
    fn critical_buckling_sign_follows_k_inf() {
        let p1 = [[0.15, 0.01], [0.0, 0.5]];
        for (nu_fission, positive) in [([0.01, 0.4], true), ([0.005, 0.2], false)] {
            let medium = two_group(p1, nu_fission);
            let k_inf = medium.mode(0.0).unwrap().k;
            assert_eq!(k_inf > 1.0, positive, "k∞ {k_inf}");
            let mode = medium.critical().unwrap();
            assert!((mode.k - 1.0).abs() < 1e-10, "k {}", mode.k);
            assert_eq!(mode.buckling > 0.0, positive, "B² {}", mode.buckling);
        }
    }

    #[test]
    fn no_fission_is_an_error() {
        let medium = two_group([[0.0; 2]; 2], [0.0, 0.0]);
        assert!(matches!(medium.critical(), Err(B1Error::NoFission)));
    }
}
