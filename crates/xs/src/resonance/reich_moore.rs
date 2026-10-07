//! Reich-Moore 形式の断面積（0 K）。
//!
//! NJOY2016 reconr の csrmat と同じ手順で計算する。捕獲は消去されたチャネルとして扱い、
//! 中性子と2つの核分裂チャネルの 3×3 の R 行列を解く。レベルシフトは使わない。

use std::f64::consts::PI;

use num_complex::Complex64;

use super::Sigma;
use super::factors::{channel_radius, penetrability, phase_shift};
use crate::endf::{OrbitalL, ReichMoore};
use crate::physics::wave_number_factor;

/// 断面積の計算に使う形に前処理した Reich-Moore パラメータ。
#[derive(Debug, Clone)]
pub struct ReichMooreXs {
    spi: f64,
    l_values: Vec<LValue>,
}

#[derive(Debug, Clone)]
struct LValue {
    l: OrbitalL,
    /// AWRI / (AWRI + 1)
    arat: f64,
    /// 透過因子に使うチャネル半径 (10^-12 cm)
    channel_radius: f64,
    /// 位相シフトに使う散乱半径 (10^-12 cm)
    scattering_radius: f64,
    resonances: Vec<Resonance>,
}

#[derive(Debug, Clone, Copy)]
struct Resonance {
    er: f64,
    aj: f64,
    /// Γn / P_l(ρ(|Er|))。E での中性子幅はこれに P_l(ρ(E)) を掛けて求める
    gn_per_penetrability: f64,
    gg: f64,
    gfa: f64,
    gfb: f64,
}

impl ReichMooreXs {
    /// 波数とチャネル半径は、l ごとの AWRI から求める。
    /// （NJOY は断面積の計算に最初の l の AWRI を使うが、各 l の AWRI が同じなら結果は同じ）
    pub fn new(data: &ReichMoore) -> Self {
        let cwaven = wave_number_factor();
        let l_values = data
            .l_values
            .iter()
            .map(|lv| {
                let scattering_radius = if lv.apl != 0.0 { lv.apl } else { data.ap };
                let channel_radius = if data.naps == 1 {
                    scattering_radius
                } else {
                    channel_radius(lv.awri)
                };
                let arat = lv.awri / (lv.awri + 1.0);
                let resonances = lv
                    .resonances
                    .iter()
                    .map(|r| {
                        let rho = cwaven * arat * r.er.abs().sqrt() * channel_radius;
                        Resonance {
                            er: r.er,
                            aj: r.aj,
                            gn_per_penetrability: r.gn / penetrability(lv.l, rho),
                            gg: r.gg,
                            gfa: r.gfa,
                            gfb: r.gfb,
                        }
                    })
                    .collect();
                LValue {
                    l: lv.l,
                    arat,
                    channel_radius,
                    scattering_radius,
                    resonances,
                }
            })
            .collect();
        Self {
            spi: data.spi,
            l_values,
        }
    }

    /// エネルギー e (eV) での断面積。
    pub fn cross_sections(&self, e: f64) -> Sigma {
        let gjd = 2.0 * (2.0 * self.spi + 1.0);
        let mut sum = Sigma::default();
        for lv in self.l_values.iter().filter(|lv| !lv.resonances.is_empty()) {
            let k = wave_number_factor() * lv.arat * e.abs().sqrt();
            let pifac = PI / (k * k);
            let pe = penetrability(lv.l, k * lv.channel_radius);
            let phi = phase_shift(lv.l, k * lv.scattering_radius);

            // l と標的スピンから作れる J を、チャネルスピン s = I ± 1/2 の両方について調べる
            let fl = f64::from(lv.l.number());
            let ajmin = ((self.spi - fl).abs() - 0.5).abs();
            let ajmax = self.spi + fl + 0.5;
            let numj = (ajmax - ajmin + 1.0).round() as usize;
            let jjl = if lv.l != OrbitalL::S && fl > self.spi - 0.5 && fl <= self.spi {
                0
            } else {
                1
            };
            for jj in 1..=numj {
                let ajc = ajmin + (jj - 1) as f64;
                let gj = (2.0 * ajc + 1.0) / gjd;
                let interior = jj > jjl && jj < numj;
                for channel in [ChannelSpin::Positive, ChannelSpin::Negative] {
                    let group = self.channel_terms(lv, e, pe, ajc, channel);
                    let hard_sphere = match channel.contribution(group.positive, group.negative) {
                        Contribution::None => continue,
                        Contribution::Resonances => false,
                        Contribution::ResonancesAndHardSphere => interior,
                    };
                    let mut terms = group.terms(gj, phi);
                    if hard_sphere {
                        // 硬球散乱 2 gJ (1 - cos 2φ)
                        let potential = 4.0 * gj * phi.sin().powi(2);
                        terms.elastic += potential;
                        terms.total += potential;
                    }
                    sum.total += pifac * terms.total;
                    sum.elastic += pifac * terms.elastic;
                    sum.fission += pifac * terms.fission;
                    sum.capture += pifac * (terms.total - terms.fission - terms.elastic);
                }
            }
        }
        sum
    }

    /// 1つの (J, チャネルスピン) に属する共鳴から R 行列を組み立てる。
    fn channel_terms(
        &self,
        lv: &LValue,
        e: f64,
        pe: f64,
        ajc: f64,
        channel: ChannelSpin,
    ) -> ChannelGroup {
        let mut group = ChannelGroup::default();
        for res in &lv.resonances {
            if (res.aj.abs() - ajc).abs() > 0.25 {
                continue;
            }
            if res.aj < 0.0 {
                group.negative += 1;
            }
            if res.aj > 0.0 {
                group.positive += 1;
            }
            if channel.excludes(res.aj) {
                continue;
            }
            // 幅の平方根に符号を付けたものが換算幅振幅
            let amplitude = |width: f64| width.abs().sqrt().copysign(width);
            let a = [
                (res.gn_per_penetrability * pe).sqrt(),
                amplitude(res.gfa),
                amplitude(res.gfb),
            ];
            let diff = res.er - e;
            let den = diff * diff + 0.25 * res.gg * res.gg;
            let de2 = 0.5 * diff / den;
            let gg4 = 0.25 * res.gg / den;
            let has_fission = res.gfa != 0.0 || res.gfb != 0.0;
            group.fission |= has_fission;
            let n = if has_fission { 3 } else { 1 };
            for i in 0..n {
                for j in i..n {
                    group.r[i][j] += gg4 * a[i] * a[j];
                    group.s[i][j] -= de2 * a[i] * a[j];
                }
            }
        }
        group
    }
}

/// NJOY が負の J で区別するチャネルスピン。
#[derive(Debug, Clone, Copy, PartialEq)]
enum ChannelSpin {
    /// J が正（または 0）の共鳴
    Positive,
    /// J が負（または 0）の共鳴
    Negative,
}

/// そのチャネルスピンの組で、断面積に足すもの。
enum Contribution {
    None,
    Resonances,
    /// 共鳴に加えて、共鳴のない側のチャネルスピンの硬球散乱も足す（J が両方のスピンで作れる場合）
    ResonancesAndHardSphere,
}

impl ChannelSpin {
    fn excludes(self, aj: f64) -> bool {
        match self {
            Self::Positive => aj < 0.0,
            Self::Negative => aj > 0.0,
        }
    }

    /// NJOY csrmat の kkkkkk の判定。
    fn contribution(self, positive: usize, negative: usize) -> Contribution {
        use Contribution::*;
        match (self, positive > 0, negative > 0) {
            (Self::Positive, _, false) => ResonancesAndHardSphere,
            (Self::Positive, true, true) => Resonances,
            (Self::Positive, false, true) => None,
            (Self::Negative, _, false) => None,
            (Self::Negative, true, true) => Resonances,
            (Self::Negative, false, true) => ResonancesAndHardSphere,
        }
    }
}

#[derive(Debug, Default)]
struct ChannelGroup {
    /// R 行列の実部（上三角）。単位行列は `terms` で足す
    r: [[f64; 3]; 3],
    /// R 行列の虚部（上三角）
    s: [[f64; 3]; 3],
    fission: bool,
    positive: usize,
    negative: usize,
}

/// gJ を掛けた断面積の項（π/k² を掛ける前）。
struct Terms {
    total: f64,
    elastic: f64,
    fission: f64,
}

impl ChannelGroup {
    /// 中性子チャネルの衝突行列 U_nn = e^{-2iφ} (2 (M^{-1})_nn - 1) から断面積の項を求める。
    /// M = I + K（K は共鳴の寄与）。
    ///
    /// 1 - Re U は低エネルギーで O(ρ²) と小さく、そのまま計算すると桁落ちする。
    /// I - M^{-1} = M^{-1} K と 1 - cos 2φ = 2 sin²φ を使って、差を取らずに計算する
    /// （NJOY は同じ目的で近似展開を使っている）。
    fn terms(&self, gj: f64, phi: f64) -> Terms {
        let n = if self.fission { 3 } else { 1 };
        // 上三角に積んだ K を対称に展開する。核分裂のないときは中性子チャネルだけを使う
        let k: [[Complex64; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let (a, b) = (i.min(j), i.max(j));
                Complex64::new(self.r[a][b], self.s[a][b])
            })
        });
        let m = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let identity = if i == j {
                    Complex64::ONE
                } else {
                    Complex64::ZERO
                };
                if i < n && j < n {
                    identity + k[i][j]
                } else {
                    identity
                }
            })
        });
        let inv = invert3(m);
        // x = (I - M^{-1})_nn = (M^{-1} K)_nn
        let x = (0..n).fold(Complex64::ZERO, |acc, j| acc + inv[0][j] * k[j][0]);
        let (p1, p2) = ((2.0 * phi).cos(), (2.0 * phi).sin());
        let one_minus_u = 2.0 * phi.sin().powi(2) + 2.0 * p1 * x.re + 2.0 * p2 * x.im;
        let u_imag = p2 * (2.0 * x.re - 1.0) - 2.0 * p1 * x.im;
        Terms {
            total: 2.0 * gj * one_minus_u,
            elastic: gj * (one_minus_u * one_minus_u + u_imag * u_imag),
            fission: 4.0 * gj * (inv[0][1].norm_sqr() + inv[0][2].norm_sqr()),
        }
    }
}

/// 複素 3×3 行列の逆行列（部分ピボットつき Gauss-Jordan 法）。
/// R 行列は単位行列に正定値の実部を足したものなので、正則である。
fn invert3(mut a: [[Complex64; 3]; 3]) -> [[Complex64; 3]; 3] {
    let mut inv = [[Complex64::ZERO; 3]; 3];
    for (i, row) in inv.iter_mut().enumerate() {
        row[i] = Complex64::ONE;
    }
    for col in 0..3 {
        // max_by と同じく、同じ大きさなら後の行を選ぶ
        let pivot = (col + 1..3).fold(col, |best, i| {
            if a[i][col]
                .norm_sqr()
                .total_cmp(&a[best][col].norm_sqr())
                .is_ge()
            {
                i
            } else {
                best
            }
        });
        a.swap(col, pivot);
        inv.swap(col, pivot);
        let p = a[col][col].inv();
        for j in 0..3 {
            a[col][j] *= p;
            inv[col][j] *= p;
        }
        for row in 0..3 {
            if row != col {
                let f = a[row][col];
                for j in 0..3 {
                    a[row][j] -= f * a[col][j];
                    inv[row][j] -= f * inv[col][j];
                }
            }
        }
    }
    inv
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::endf::{RmLValue, RmResonance};
    use crate::physics::wave_number_factor;
    use std::f64::consts::PI;

    pub(crate) const AWRI: f64 = 236.006;
    pub(crate) const AP: f64 = 0.948;

    pub(crate) fn single_level(gn: f64, gg: f64, gf: f64) -> ReichMoore {
        ReichMoore {
            spi: 0.0,
            ap: AP,
            naps: 1,
            l_values: vec![RmLValue {
                awri: AWRI,
                apl: 0.0,
                l: OrbitalL::S,
                resonances: vec![RmResonance {
                    er: 6.674,
                    aj: 0.5,
                    gn,
                    gg,
                    gfa: gf,
                    gfb: 0.0,
                }],
            }],
        }
    }

    /// スピン 0 の標的、l=0 の単一準位での厳密解（Breit-Wigner 型）。
    pub(crate) fn single_level_exact(e: f64, gn: f64, gg: f64, gf: f64) -> Sigma {
        let er: f64 = 6.674;
        let k = wave_number_factor() * AWRI / (AWRI + 1.0) * e.sqrt();
        let kr = wave_number_factor() * AWRI / (AWRI + 1.0) * er.sqrt();
        let pifac = PI / (k * k);
        let phi = k * AP;
        // l=0 の透過因子は ρ なので、中性子幅は ρ に比例する
        let gn_e = gn * (k * AP) / (kr * AP);
        let width = gn_e + gg + gf;
        let d = er - e;
        // U = e^{-2iφ} (1 + iΓn / (d - iΓ/2))
        let denom = d * d + width * width / 4.0;
        let (qr, qi) = (-gn_e * width / 2.0 / denom, gn_e * d / denom);
        let (wr, wi) = (1.0 + qr, qi);
        let (c, s) = ((2.0 * phi).cos(), -(2.0 * phi).sin());
        let ui = c * wi + s * wr;
        // 1 - Re U を差を取らずに計算する（低エネルギーでの桁落ちを避ける）
        let one_minus_ur =
            2.0 * phi.sin().powi(2) - (2.0 * phi).cos() * qr - (2.0 * phi).sin() * qi;
        let total = 2.0 * pifac * one_minus_ur;
        let elastic = pifac * (one_minus_ur.powi(2) + ui * ui);
        let fission = pifac * gn_e * gf / denom;
        Sigma {
            total,
            elastic,
            fission,
            capture: total - elastic - fission,
        }
    }

    fn assert_close(actual: f64, expected: f64, what: &str) {
        assert!(
            (actual - expected).abs() <= 1e-9 * expected.abs().max(1e-6),
            "{what}: actual {actual}, expected {expected}"
        );
    }

    fn assert_sigma(actual: Sigma, expected: Sigma, e: f64) {
        assert_close(actual.total, expected.total, &format!("total at {e}"));
        assert_close(actual.elastic, expected.elastic, &format!("elastic at {e}"));
        assert_close(actual.fission, expected.fission, &format!("fission at {e}"));
        assert_close(actual.capture, expected.capture, &format!("capture at {e}"));
    }

    #[test]
    fn singular_level_matrix_is_an_error() {
        // M = I + K が特異になる K。物理的な共鳴パラメータ（K の実部が半正定値）では起きない
        let mut group = ChannelGroup::default();
        group.r[0][0] = -1.0;
        assert!(group.terms(0.5, 0.1).is_err());
    }

    #[test]
    fn single_level_without_fission_matches_breit_wigner() {
        let xs = ReichMooreXs::new(&single_level(1.5e-3, 2.3e-2, 0.0));
        for e in [1.0e-3, 1.0, 6.0, 6.674, 6.69, 7.0, 100.0, 1.0e4] {
            assert_sigma(
                xs.cross_sections(e),
                single_level_exact(e, 1.5e-3, 2.3e-2, 0.0),
                e,
            );
        }
    }

    #[test]
    fn single_level_with_fission_matches_breit_wigner() {
        let xs = ReichMooreXs::new(&single_level(1.5e-3, 2.3e-2, 1.0e-1));
        for e in [1.0e-3, 1.0, 6.674, 6.8, 100.0] {
            assert_sigma(
                xs.cross_sections(e),
                single_level_exact(e, 1.5e-3, 2.3e-2, 1.0e-1),
                e,
            );
        }
    }

    #[test]
    fn negative_fission_amplitude_gives_same_single_level_result() {
        // 単一準位では核分裂幅の符号（振幅の符号）は結果に影響しない
        let xs = ReichMooreXs::new(&single_level(1.5e-3, 2.3e-2, -1.0e-1));
        let e = 6.7;
        assert_sigma(
            xs.cross_sections(e),
            single_level_exact(e, 1.5e-3, 2.3e-2, 1.0e-1),
            e,
        );
    }

    #[test]
    fn awri_of_another_l_does_not_change_result() {
        // NAPS=0（チャネル半径を AWRI から求める）で、共鳴のない l の AWRI が違っても、
        // 共鳴のある l の断面積は変わらない
        let mut alone = single_level(1.5e-3, 2.3e-2, 1.0e-2);
        alone.naps = 0;
        let mut with_other = alone.clone();
        with_other.l_values.insert(
            0,
            RmLValue {
                awri: 1.0,
                apl: 0.0,
                l: OrbitalL::P,
                resonances: vec![],
            },
        );
        for e in [1.0, 6.674, 100.0] {
            assert_eq!(
                ReichMooreXs::new(&with_other).cross_sections(e),
                ReichMooreXs::new(&alone).cross_sections(e),
                "at {e}"
            );
        }
    }

    #[test]
    fn p_wave_potential_scattering_sums_over_all_j() {
        // スピン 0 の標的の l=1 では J=1/2, 3/2 があり、gJ の和は 2l+1 = 3。
        // 共鳴から遠いエネルギーでは、どちらの J も硬球散乱 (4π/k²) gJ sin²φ1 だけになる
        let data = ReichMoore {
            spi: 0.0,
            ap: AP,
            naps: 1,
            l_values: vec![RmLValue {
                awri: AWRI,
                apl: 0.0,
                l: OrbitalL::P,
                resonances: vec![RmResonance {
                    er: 1.0e4,
                    aj: 1.5,
                    gn: 1.0e-9,
                    gg: 1.0e-9,
                    gfa: 0.0,
                    gfb: 0.0,
                }],
            }],
        };
        let e: f64 = 100.0;
        let k = wave_number_factor() * AWRI / (AWRI + 1.0) * e.sqrt();
        let rho = k * AP;
        let phi = rho - rho.atan();
        let expected = 4.0 * PI / (k * k) * 3.0 * phi.sin().powi(2);
        let sigma = ReichMooreXs::new(&data).cross_sections(e);
        assert!((sigma.elastic - expected).abs() <= 1e-6 * expected);
        assert!((sigma.total - expected).abs() <= 1e-6 * expected);
    }
}
