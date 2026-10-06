//! 非分離共鳴の無限希釈の平均断面積（NJOY RECONR csunr2 と同じ単準位 Breit-Wigner の式）。
//!
//! 幅のゆらぎの期待値は、NJOY の求積表ではなく適応積分で求める（`fluctuation`）。

use std::f64::consts::PI;

use super::Sigma;
use super::factors::{channel_radius, penetrability, phase_shift};
use super::fluctuation::{expectation, expectation_vec, width_ratios};
use super::line_shape::jk_interpolated;
use crate::endf::{Interpolation, Unresolved, UrrJValue, UrrPoint, WidthDof};
use crate::physics::{BK, wave_number_factor};

/// 非分離共鳴のパラメータから、無限希釈の平均断面積を計算する。
pub struct UnresolvedXs<'a> {
    urr: &'a Unresolved,
}

/// 非分離共鳴の計算のエラー。
#[derive(Debug, thiserror::Error)]
pub enum UrrError {
    #[error("NAPS={0} is not allowed in the unresolved range without NRO=1")]
    Naps(i64),
    #[error("self-shielding needs a positive temperature, got {0}")]
    Temperature(f64),
    #[error(
        "background {background} b at {energy} eV is not positive \
         (potential and interference terms {terms})"
    )]
    Background {
        background: f64,
        energy: f64,
        terms: f64,
    },
    #[error("β={beta:e}, θ={theta:e} is outside the J, K table")]
    JkTable { beta: f64, theta: f64 },
    #[error("unresolved range has no (l, J) sequences")]
    NoSequences,
    #[error("unresolved sequences use different interpolation laws ({first:?} and {other:?})")]
    MixedLaws {
        first: Interpolation,
        other: Interpolation,
    },
}

impl<'a> UnresolvedXs<'a> {
    /// チャネル半径の選び方（NAPS）を確かめる。NAPS=2 は散乱半径の表（NRO=1）があるときだけ使える。
    pub fn new(urr: &'a Unresolved) -> Result<Self, UrrError> {
        match (urr.naps, &urr.ap_energy) {
            (0 | 1, _) | (2, Some(_)) => Ok(Self { urr }),
            (naps, _) => Err(UrrError::Naps(naps)),
        }
    }

    /// エネルギー e (eV) での平均断面積。パラメータは lin-lin で補間する（NJOY と同じ）。
    ///
    /// 波数は l ごとの AWRI から求める（NJOY は最初の l の AWRI を使う。全 l で同じなら一致する）。
    pub fn cross_sections(&self, e: f64) -> Sigma {
        let urr = self.urr;
        let ay = self.scattering_radius(e);
        let mut sigma = Sigma::default();
        for lv in &urr.l_values {
            let aa = self.channel_radius(lv.awri, ay);
            let k = wave_number_factor() * lv.awri / (lv.awri + 1.0) * e.sqrt();
            let (rho, rhoc) = (k * aa, k * ay);
            let sin2 = phase_shift(lv.l, rhoc).sin().powi(2);
            // 1/k² に掛ける係数
            let pifac = PI / (k * k);
            if !lv.j_values.is_empty() {
                sigma.elastic += 4.0 * pifac * f64::from(2 * lv.l.number() + 1) * sin2;
            }
            for jv in &lv.j_values {
                let p = parameters(&jv.points, e);
                let gj = (2.0 * jv.aj + 1.0) / (4.0 * urr.spi + 2.0);
                // 平均中性子幅（換算幅 × 自由度 × 透過因子/ρ × √E）
                let gn =
                    p.gno * f64::from(jv.amun.number()) * penetrability(lv.l, rho) / rho * e.sqrt();
                let scale = 2.0 * PI * pifac * gj / p.d;
                let [f1, f2, f3] = fluctuation_integrals(gn, p.gf, p.gg, p.gx, jv);
                sigma.elastic += scale * (gn * gn * f1 - 2.0 * gn * sin2);
                sigma.capture += scale * gn * p.gg * f2;
                sigma.fission += scale * gn * p.gf * f3;
            }
        }
        sigma.total = sigma.elastic + sigma.fission + sigma.capture;
        sigma
    }
}

/// 非分離共鳴の自己遮蔽した断面積 (barn)。`current_total` は流れで重みづけした全断面積
/// （NJOY の MT261 に相当し、P1 の全断面積に使う）。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ShieldedSigma {
    pub total: f64,
    pub elastic: f64,
    pub fission: f64,
    pub capture: f64,
    pub current_total: f64,
}

impl UnresolvedXs<'_> {
    /// エネルギー e (eV)、温度 `temperature` (K) での、背景断面積 `sigma0` ごとの実効断面積
    /// （NJOY UNRESR unresl と同じ ETOX 法）。`background` は MF3 の背景断面積
    /// （LSSF=1 では MF3 が平均断面積そのものなので 0 を渡す）。
    pub fn shielded(
        &self,
        e: f64,
        temperature: f64,
        sigma0: &[f64],
        background: Sigma,
    ) -> Result<Vec<ShieldedSigma>, UrrError> {
        if temperature.is_nan() || temperature <= 0.0 {
            return Err(UrrError::Temperature(temperature));
        }
        let urr = self.urr;
        let ay = self.scattering_radius(e);
        // ポテンシャル散乱と干渉の補正、系列ごとの平均の幅
        let (mut spot, mut sint) = (0.0, 0.0);
        let mut sequences = Vec::new();
        for lv in &urr.l_values {
            let aa = self.channel_radius(lv.awri, ay);
            let k = wave_number_factor() * lv.awri / (lv.awri + 1.0) * e.sqrt();
            let (rho, rhoc) = (k * aa, k * ay);
            let sin2 = phase_shift(lv.l, rhoc).sin().powi(2);
            let ab = 4.0 * PI / (k * k);
            if !lv.j_values.is_empty() {
                spot += ab * f64::from(2 * lv.l.number() + 1) * sin2;
            }
            // ドップラー幅 Δ = 2 √(E kT / A)
            let delta = 2.0 * (e * BK * temperature / lv.awri).sqrt();
            for jv in &lv.j_values {
                let p = parameters(&jv.points, e);
                let gj = (2.0 * jv.aj + 1.0) / (4.0 * urr.spi + 2.0);
                let gn =
                    p.gno * f64::from(jv.amun.number()) * penetrability(lv.l, rho) / rho * e.sqrt();
                sint -= ab * PI * gj * gn * sin2 / p.d;
                sequences.push(Sequence {
                    gj,
                    ab,
                    gn,
                    gf: p.gf,
                    gg: p.gg,
                    gx: p.gx,
                    d: p.d,
                    delta,
                    jv,
                });
            }
        }
        let sigbt = background.total + spot + sint;
        let sigm: Vec<f64> = sigma0.iter().map(|s0| sigbt + s0).collect();
        if let Some(bad) = sigm.iter().find(|s| s.is_nan() || **s <= 0.0) {
            return Err(UrrError::Background {
                background: *bad,
                energy: e,
                terms: spot + sint,
            });
        }
        let integrals = sequences
            .iter()
            .map(|q| q.integrals(&sigm))
            .collect::<Result<Vec<_>, _>>()?;

        // 系列どうしの重なりを補正して足し合わせる（NJOY unresl と同じ式）
        Ok((0..sigma0.len())
            .map(|s| {
                let t = |q: usize, kx: usize| integrals[q][PARTS * s + kx];
                let tj = |q: usize| (0..4).map(|kx| t(q, kx)).sum::<f64>();
                let all_j: f64 = (0..integrals.len()).map(tj).sum();
                let all_k: f64 = (0..integrals.len()).map(|q| t(q, 4)).sum();
                let (mut yy, mut yj, mut yk) = ([0.0; 3], 0.0, 0.0);
                for q in 0..integrals.len() {
                    let others_j = all_j - tj(q);
                    let others_k = all_k - t(q, 4);
                    for (kx, y) in yy.iter_mut().enumerate() {
                        *y += t(q, kx) * (1.0 - others_j);
                    }
                    let own = tj(q) * (1.0 - others_j);
                    yj += own;
                    yk += (t(q, 4) - tj(q)) * (1.0 - others_k) + own;
                }
                let [fission, capture, elastic] = yy.map(|y| sigm[s] * y / (1.0 - yj));
                ShieldedSigma {
                    total: fission + capture + elastic + sigbt,
                    elastic: elastic + background.elastic + spot + sint,
                    fission: fission + background.fission,
                    capture: capture + background.capture,
                    current_total: sigbt * (1.0 - yj) / (1.0 - yk)
                        + sigma0[s] * (yk - yj) / (1.0 - yk),
                }
            })
            .collect())
    }

    /// 位相シフトに使う散乱半径（NRO=1 なら表から補間する）。
    fn scattering_radius(&self, e: f64) -> f64 {
        match &self.urr.ap_energy {
            Some(table) => table
                .eval(e)
                .expect("radius table covers the unresolved range"),
            None => self.urr.ap,
        }
    }

    /// 透過因子に使うチャネル半径（NAPS で選ぶ）。
    fn channel_radius(&self, awri: f64, ay: f64) -> f64 {
        match self.urr.naps {
            0 => channel_radius(awri),
            1 => ay,
            // NAPS=2（散乱半径の表があることは `new` で確かめてある）
            _ => self.urr.ap,
        }
    }
}

/// 自己遮蔽の積分で、σ0 ごとに求める量の数: J Γf, J Γγ, J Γn, J Γx, K Γ
const PARTS: usize = 5;

/// 幅のゆらぎの期待値の相対許容誤差（中性子幅、核分裂幅、競合幅の順）。外側は内側の 100 倍にする
/// （`expectation` を参照）。J, K の表の補間は連続なので、最も内側を緩める必要はない
const SHIELDING_TOLERANCE: [f64; 3] = [1.0e-10, 1.0e-8, 1.0e-6];

/// 1つの (l, J) 系列の、エネルギー e での平均の幅など。
struct Sequence<'a> {
    gj: f64,
    /// 4π/k²
    ab: f64,
    gn: f64,
    gf: f64,
    gg: f64,
    gx: f64,
    d: f64,
    /// ドップラー幅
    delta: f64,
    jv: &'a UrrJValue,
}

impl Sequence<'_> {
    /// σ0 ごとの [E[J Γf], E[J Γγ], E[J Γn], E[J Γx], E[K Γ]] / D。`sigm` は σ0 を含む背景断面積。
    /// J, K は表から補間する（`jk_interpolated`）。表の範囲を外れたら、その点を記録してエラーにする。
    fn integrals(&self, sigm: &[f64]) -> Result<Vec<f64>, UrrError> {
        let n = PARTS * sigm.len();
        // 平均断面積（`fluctuation_integrals`）と同じく、Γn または Γγ が 0 以下なら 0 とする
        if self.gn <= 0.0 || self.gg <= 0.0 {
            return Ok(vec![0.0; n]);
        }
        let dof_f = if self.gf > 0.0 {
            self.jv.amuf
        } else {
            WidthDof::Fixed
        };
        let dof_x = if self.gx > 0.0 {
            self.jv.amux
        } else {
            WidthDof::Fixed
        };
        let failure = std::cell::RefCell::new(None);
        let node = |x: f64, y: f64, z: f64, out: &mut [f64]| {
            let (gn, gf, gx) = (self.gn * x, self.gf * y, self.gx * z);
            let gamma = gn + gf + self.gg + gx;
            // 共鳴のピークの断面積と、θ = Γ/Δ
            let peak = self.ab * self.gj * gn / gamma;
            let theta = gamma / self.delta;
            for (i, s) in sigm.iter().enumerate() {
                let (j, k) = match jk_interpolated(s / peak, theta) {
                    Ok(jk) => jk,
                    Err(message) => {
                        failure.borrow_mut().get_or_insert(message);
                        (0.0, 0.0)
                    }
                };
                out[PARTS * i..PARTS * (i + 1)].copy_from_slice(&[
                    j * gf,
                    j * self.gg,
                    j * gn,
                    j * gx,
                    k * gamma,
                ]);
            }
        };
        let [tol_n, tol_f, tol_x] = SHIELDING_TOLERANCE;
        let values = expectation_vec(dof_x, tol_x, n, |z, out| {
            let inner = expectation_vec(dof_f, tol_f, n, |y, out| {
                let inner = expectation_vec(self.jv.amun, tol_n, n, |x, out| node(x, y, z, out));
                out.copy_from_slice(&inner);
            });
            out.copy_from_slice(&inner);
        });
        match failure.into_inner() {
            Some(message) => Err(message),
            None => Ok(values.into_iter().map(|v| v / self.d).collect()),
        }
    }
}

/// エネルギー e でのパラメータ（lin-lin 補間）。表の外では端の区間を延長する（NJOY と同じ）。
/// 同じエネルギーが続く段差の上では、右側の値を使う。
fn parameters(points: &[UrrPoint], e: f64) -> UrrPoint {
    if points.len() == 1 {
        return points[0];
    }
    let k = points
        .partition_point(|p| p.e <= e)
        .clamp(1, points.len() - 1);
    let (a, b) = (points[k - 1], points[k]);
    if b.e == a.e {
        return b;
    }
    let t = (e - a.e) / (b.e - a.e);
    let lerp = |x: f64, y: f64| x + t * (y - x);
    UrrPoint {
        e,
        d: lerp(a.d, b.d),
        gx: lerp(a.gx, b.gx),
        gno: lerp(a.gno, b.gno),
        gg: lerp(a.gg, b.gg),
        gf: lerp(a.gf, b.gf),
    }
}

/// 幅のゆらぎの積分 [E[x²/D], E[x/D], E[xy/D]]、D = Γn x + Γf y + Γγ + Γx z。
/// x, y, z は自由度 AMUN, AMUF, AMUX の χ² 分布に従う（幅が 0 ならその変数はゆらがない）。
/// 中性子幅 x の期待値は閉じた式（`width_ratios`）で、y と z は適応積分で求める。
///
/// NJOY と同じく、Γn または Γγ が 0 以下なら 0 とする（共鳴がない、または分母が 0 になる）。
fn fluctuation_integrals(gn: f64, gf: f64, gg: f64, gx: f64, jv: &UrrJValue) -> [f64; 3] {
    if gn <= 0.0 || gg <= 0.0 {
        return [0.0; 3];
    }
    let dof_f = if gf > 0.0 { jv.amuf } else { WidthDof::Fixed };
    let dof_x = if gx > 0.0 { jv.amux } else { WidthDof::Fixed };
    let ratios = |y: f64, z: f64| width_ratios(jv.amun, (gf * y + gg + gx * z) / gn);
    let outer = |g: &dyn Fn(f64, f64) -> f64| {
        expectation(dof_x, TOLERANCE[1], |z| {
            expectation(dof_f, TOLERANCE[0], |y| g(y, z))
        }) / gn
    };
    [
        outer(&|y, z| ratios(y, z)[1]),
        outer(&|y, z| ratios(y, z)[0]),
        if gf > 0.0 {
            outer(&|y, z| y * ratios(y, z)[0])
        } else {
            0.0
        },
    ]
}

/// 幅のゆらぎの期待値の相対許容誤差（内側の核分裂幅、外側の競合幅の順）。
/// 外側は内側の 100 倍にする（`expectation` を参照）。群定数の許容誤差 1e-3 より十分小さい
const TOLERANCE: [f64; 2] = [1.0e-10, 1.0e-8];

/// 平均断面積を計算する格子 [el, eh]。評価者のエネルギー（全系列の和）に、隣の点との比が
/// `WIDE` を超える区間では NJOY と同じ補助格子 `FILL` の点を加える。
/// NJOY は最初の系列のエネルギーだけを使うが、系列ごとに違う場合もパラメータの折れ点を落とさないよう和をとる。
pub fn unresolved_grid(urr: &Unresolved, el: f64, eh: f64) -> Vec<f64> {
    let mut energies: Vec<f64> = urr
        .l_values
        .iter()
        .flat_map(|lv| &lv.j_values)
        .flat_map(|jv| jv.points.iter().map(|p| p.e))
        .filter(|&e| e > el && e < eh)
        .chain([el, eh])
        .collect();
    energies.sort_by(f64::total_cmp);
    energies.dedup();
    let mut grid = vec![energies[0]];
    for w in energies.windows(2) {
        let (lo, hi) = (w[0], w[1]);
        if hi > WIDE * lo {
            grid.extend(
                FILL.iter()
                    .filter(|&&g| g > lo * (1.0 + FILL_MARGIN) && g < hi),
            );
        }
        grid.push(hi);
    }
    grid
}

/// 格子点の間の断面積の補間則。系列ごとに違うと、どれを使うかの根拠がないのでエラーにする。
pub fn interpolation(urr: &Unresolved) -> Result<Interpolation, UrrError> {
    let mut laws = urr
        .l_values
        .iter()
        .flat_map(|lv| &lv.j_values)
        .map(|jv| jv.int);
    let first = laws.next().ok_or(UrrError::NoSequences)?;
    match laws.find(|&law| law != first) {
        Some(other) => Err(UrrError::MixedLaws { first, other }),
        None => Ok(first),
    }
}

/// 隣り合う格子点のエネルギー比がこれを超えると、補助格子の点を加える（NJOY と同じ）
const WIDE: f64 = 1.26;

/// 補助格子の点は、区間の下端からこの割合以上離れたものだけを使う（NJOY と同じ）
const FILL_MARGIN: f64 = 1.0e-3;

/// NJOY RECONR の補助格子 (eV)
const FILL: [f64; 78] = [
    1.0e1, 1.25e1, 1.5e1, 1.7e1, 2.0e1, 2.5e1, 3.0e1, 3.5e1, 4.0e1, 5.0e1, 6.0e1, 7.2e1, 8.5e1,
    1.0e2, 1.25e2, 1.5e2, 1.7e2, 2.0e2, 2.5e2, 3.0e2, 3.5e2, 4.0e2, 5.0e2, 6.0e2, 7.2e2, 8.5e2,
    1.0e3, 1.25e3, 1.5e3, 1.7e3, 2.0e3, 2.5e3, 3.0e3, 3.5e3, 4.0e3, 5.0e3, 6.0e3, 7.2e3, 8.5e3,
    1.0e4, 1.25e4, 1.5e4, 1.7e4, 2.0e4, 2.5e4, 3.0e4, 3.5e4, 4.0e4, 5.0e4, 6.0e4, 7.2e4, 8.5e4,
    1.0e5, 1.25e5, 1.5e5, 1.7e5, 2.0e5, 2.5e5, 3.0e5, 3.5e5, 4.0e5, 5.0e5, 6.0e5, 7.2e5, 8.5e5,
    1.0e6, 1.25e6, 1.5e6, 1.7e6, 2.0e6, 2.5e6, 3.0e6, 3.5e6, 4.0e6, 5.0e6, 6.0e6, 7.2e6, 8.5e6,
];

#[cfg(test)]
mod tests {
    use super::*;
    use crate::endf::{OrbitalL, UrrLValue};
    use crate::special::erfc;

    const AWRI: f64 = 236.0;
    const AP: f64 = 0.948;

    fn point(e: f64, gno: f64) -> UrrPoint {
        UrrPoint {
            e,
            d: 20.0,
            gx: 0.0,
            gno,
            gg: 0.023,
            gf: 0.0,
        }
    }

    fn s_wave(points: Vec<UrrPoint>, int: Interpolation) -> Unresolved {
        Unresolved {
            lssf: 0,
            spi: 0.0,
            ap: AP,
            ap_energy: None,
            naps: 1,
            l_values: vec![UrrLValue {
                awri: AWRI,
                l: OrbitalL::S,
                j_values: vec![UrrJValue {
                    aj: 0.5,
                    int,
                    amux: WidthDof::Fixed,
                    amun: WidthDof::One,
                    amuf: WidthDof::Fixed,
                    points,
                }],
            }],
        }
    }

    /// 平均 1 の χ²_1 での E[1/(x+c)] = c^(-1/2) e^(c/2) √(π/2) erfc(√(c/2))。
    fn inverse_moment(c: f64) -> f64 {
        (PI / (2.0 * c)).sqrt() * (0.5 * c).exp() * erfc((0.5 * c).sqrt())
    }

    #[test]
    fn s_wave_without_fission_matches_closed_form() {
        // GNO は 5e3 eV で 1e-3、2e4 eV で 4e-3 なので、1e4 eV では lin-lin で 2e-3
        let urr = s_wave(
            vec![point(5.0e3, 1.0e-3), point(2.0e4, 4.0e-3)],
            Interpolation::LinLin,
        );
        let e = 1.0e4;
        let sigma = UnresolvedXs::new(&urr).unwrap().cross_sections(e);

        let k = wave_number_factor() * AWRI / (AWRI + 1.0) * e.sqrt();
        let (d, gg) = (20.0, 0.023);
        let gn = 2.0e-3 * e.sqrt();
        // 中性子幅 gn x、分母 gn x + gg の期待値（c = gg/gn）
        let c = gg / gn;
        let e_inv = inverse_moment(c);
        let e_x = (1.0 - c * e_inv) / gn;
        let e_x2 = (1.0 - c + c * c * e_inv) / gn;
        let scale = 2.0 * PI * PI / (k * k) / d;
        let phase = (k * AP).sin().powi(2);
        let capture = scale * gn * gg * e_x;
        let elastic =
            scale * gn * gn * e_x2 - scale * 2.0 * gn * phase + 4.0 * PI * phase / (k * k);
        for (actual, expected) in [(sigma.capture, capture), (sigma.elastic, elastic)] {
            assert!(
                (actual - expected).abs() <= 1e-8 * expected,
                "{actual} vs {expected}"
            );
        }
        assert_eq!(sigma.fission, 0.0);
        assert!((sigma.total - (sigma.elastic + sigma.capture)).abs() <= 1e-12 * sigma.total);
    }

    #[test]
    fn grid_fills_wide_steps_like_njoy() {
        let urr = s_wave(
            vec![
                point(1.5e3, 1.0e-3),
                point(3.0e3, 1.0e-3),
                point(1.0e4, 1.0e-3),
                point(1.0e5, 1.0e-3),
            ],
            Interpolation::LinLin,
        );
        let expected = vec![
            1.5e3, 1.7e3, 2.0e3, 2.5e3, 3.0e3, 3.5e3, 4.0e3, 5.0e3, 6.0e3, 7.2e3, 8.5e3, 1.0e4,
            1.25e4, 1.5e4, 1.7e4, 2.0e4, 2.5e4, 3.0e4, 3.5e4, 4.0e4, 5.0e4, 6.0e4, 7.2e4, 8.5e4,
            1.0e5,
        ];
        assert_eq!(unresolved_grid(&urr, 1.5e3, 1.0e5), expected);
    }

    #[test]
    fn fluctuating_fission_width_converges() {
        // U-234（ENDF/B-VIII.0）の l=1, J=1/2 の系列。中性子幅と核分裂幅の入れ子の積分で、
        // 内側の積分のばらつきのために外側が収束しなかった
        let point = |e, gno, gf| UrrPoint {
            e,
            d: 10.6,
            gx: 0.0,
            gno,
            gg: 0.02,
            gf,
        };
        let urr = Unresolved {
            lssf: 0,
            spi: 0.0,
            ap: 0.893,
            ap_energy: None,
            naps: 0,
            l_values: vec![UrrLValue {
                awri: 232.029,
                l: OrbitalL::P,
                j_values: vec![UrrJValue {
                    aj: 0.5,
                    int: Interpolation::LinLin,
                    amux: WidthDof::Fixed,
                    amun: WidthDof::One,
                    amuf: WidthDof::One,
                    points: vec![
                        point(4.0e4, 1.5872e-3, 6.7737e-4),
                        point(6.0e4, 1.4566e-3, 6.0163e-4),
                    ],
                }],
            }],
        };
        let sigma = UnresolvedXs::new(&urr).unwrap().cross_sections(5.0e4);
        assert!(sigma.fission > 0.0 && sigma.fission.is_finite());
    }

    #[test]
    fn three_fluctuating_widths_converge() {
        // U-236（ENDF/B-VIII.0）の l=2, J=3/2 の系列。中性子・核分裂・競合の幅がゆらぐ
        let point = |e, gx| UrrPoint {
            e,
            d: 9.18,
            gx,
            gno: 1.836e-3,
            gg: 0.019,
            gf: 3.4e-4,
        };
        let urr = Unresolved {
            lssf: 0,
            spi: 0.0,
            ap: 0.93543,
            ap_energy: None,
            naps: 0,
            l_values: vec![UrrLValue {
                awri: 234.0178,
                l: OrbitalL::D,
                j_values: vec![UrrJValue {
                    aj: 1.5,
                    int: Interpolation::LogLog,
                    amux: WidthDof::One,
                    amun: WidthDof::One,
                    amuf: WidthDof::One,
                    points: vec![
                        point(4.5e4, 0.0),
                        point(5.0e4, 2.03e-5),
                        point(6.0e4, 3.33e-4),
                    ],
                }],
            }],
        };
        let sigma = UnresolvedXs::new(&urr).unwrap().cross_sections(5.0e4);
        assert!(sigma.fission > 0.0 && sigma.fission.is_finite());
    }

    #[test]
    fn shielded_cross_sections_at_infinite_dilution_match_averages() {
        let urr = s_wave(
            vec![point(5.0e3, 1.0e-3), point(2.0e4, 4.0e-3)],
            Interpolation::LinLin,
        );
        let xs = UnresolvedXs::new(&urr).unwrap();
        let e = 1.0e4;
        let average = xs.cross_sections(e);
        let [dilute] = xs.shielded(e, 300.0, &[1.0e10], Sigma::default()).unwrap()[..] else {
            panic!("one σ0");
        };
        for (shielded, infinite) in [
            (dilute.capture, average.capture),
            (dilute.elastic, average.elastic),
            (dilute.total, average.total),
            (dilute.current_total, average.total),
        ] {
            assert!(
                (shielded - infinite).abs() <= 1e-6 * infinite,
                "{shielded} vs {infinite}"
            );
        }
    }

    #[test]
    fn sequence_without_capture_width_is_skipped_as_in_averages() {
        // Γγ = 0 の系列は、平均断面積と同じく 0 として扱う（x → 0 で θ = Γ/Δ が 0 に近づくため）
        let mut zero_capture = point(5.0e3, 1.0e-3);
        zero_capture.gg = 0.0;
        let urr = s_wave(vec![zero_capture], Interpolation::LinLin);
        let xs = UnresolvedXs::new(&urr).unwrap();
        let average = xs.cross_sections(1.0e4);
        let [dilute] = xs
            .shielded(1.0e4, 300.0, &[1.0e10], Sigma::default())
            .unwrap()[..]
        else {
            panic!("one σ0");
        };
        assert!((dilute.elastic - average.elastic).abs() <= 1e-12 * average.elastic);
        assert_eq!(dilute.capture, 0.0);
    }

    #[test]
    fn capture_is_shielded_more_at_small_background_and_low_temperature() {
        let urr = s_wave(
            vec![point(5.0e3, 1.0e-3), point(2.0e4, 4.0e-3)],
            Interpolation::LinLin,
        );
        let xs = UnresolvedXs::new(&urr).unwrap();
        let sigma0 = [1.0e10, 1.0e3, 10.0, 1.0];
        let cold = xs
            .shielded(1.0e4, 300.0, &sigma0, Sigma::default())
            .unwrap();
        let hot = xs
            .shielded(1.0e4, 1200.0, &sigma0, Sigma::default())
            .unwrap();
        for w in cold.windows(2) {
            assert!(
                w[1].capture < w[0].capture,
                "{} then {}",
                w[0].capture,
                w[1].capture
            );
        }
        assert!(hot[2].capture > cold[2].capture);
    }

    #[test]
    fn naps_two_needs_radius_table() {
        // NAPS=2 は、散乱半径の表（NRO=1）とは別に AP をチャネル半径に使う指定
        let mut urr = s_wave(vec![point(1.5e3, 1.0e-3)], Interpolation::LinLin);
        urr.naps = 2;
        assert!(UnresolvedXs::new(&urr).is_err());
    }

    #[test]
    fn interpolation_must_agree_across_sequences() {
        let mut urr = s_wave(
            vec![point(1.5e3, 1.0e-3), point(1.0e5, 1.0e-3)],
            Interpolation::LogLog,
        );
        assert_eq!(interpolation(&urr).unwrap(), Interpolation::LogLog);
        let mut other = urr.l_values[0].j_values[0].clone();
        other.int = Interpolation::LinLin;
        urr.l_values[0].j_values.push(other);
        assert!(interpolation(&urr).is_err());
    }
}
