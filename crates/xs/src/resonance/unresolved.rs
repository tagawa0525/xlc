//! 非分離共鳴の無限希釈の平均断面積（NJOY RECONR csunr2 と同じ単準位 Breit-Wigner の式）と、ラダーのサンプリングによる
//! 自己遮蔽（ADR-0177。`ladder`）。
//!
//! 幅のゆらぎの期待値は、NJOY の求積表ではなく適応積分で求める（`fluctuation`）。

use std::f64::consts::PI;

use super::Sigma;
use super::factors::{channel_radius, penetrability, phase_shift};
use super::fluctuation::{expectation, width_ratios};
use super::ladder::{Bondarenko, LadderError, LadderSequence, Sampling, bondarenko};
use crate::endf::{Interpolation, Unresolved, UrrJValue, UrrPoint, WidthDof, scattering_radius_at};
use crate::numeric::NumericError;
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
        "the sampled total cross section is not positive {offset} eV from {energy} eV{}",
        if *floored { " (with the elastic floored at 1 microbarn)" } else { "" }
    )]
    NonPositiveTotal {
        energy: f64,
        offset: f64,
        floored: bool,
    },
    #[error(
        "the ladder sampling needs at least one ladder and one energy, got {ladders} and {samples}"
    )]
    EmptySampling { ladders: usize, samples: usize },
    #[error("unresolved range has no (l, J) sequences")]
    NoSequences,
    #[error("the scattering radius table does not cover {energy} eV")]
    RadiusTable { energy: f64 },
    #[error(transparent)]
    Numeric(#[from] NumericError),
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
    pub fn cross_sections(&self, e: f64) -> Result<Sigma, UrrError> {
        let urr = self.urr;
        let ay = self.scattering_radius(e)?;
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
                let [f1, f2, f3] = fluctuation_integrals(gn, p.gf, p.gg, p.gx, jv)?;
                sigma.elastic += scale * (gn * gn * f1 - 2.0 * gn * sin2);
                sigma.capture += scale * gn * p.gg * f2;
                sigma.fission += scale * gn * p.gf * f3;
            }
        }
        sigma.total = sigma.elastic + sigma.fission + sigma.capture;
        Ok(sigma)
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
    /// エネルギー e (eV)、温度 `temperature` (K) での、バックグラウンド断面積 `sigma0` ごとの実効断面積を、ラダーのサンプリングの
    /// Bondarenko の平均で求める（ADR-0177。`ladder`）。`background` は MF3 のバックグラウンド断面積（LSSF=1 では MF3 が平均断面積そのもの
    /// なので 0 を渡す）。サンプルの数は `sampling`。乱数は種 `seed` から決まり、同じ種なら温度を変えても同じラダーとエネルギーを使う（温度の間の相関を保つ）。
    pub fn shielded(
        &self,
        e: f64,
        temperature: f64,
        sigma0: &[f64],
        background: Sigma,
        sampling: Sampling,
        seed: u64,
    ) -> Result<Bondarenko, UrrError> {
        if temperature.is_nan() || temperature <= 0.0 {
            return Err(UrrError::Temperature(temperature));
        }
        let urr = self.urr;
        let ay = self.scattering_radius(e)?;
        let mut potential = 0.0;
        let mut sequences = Vec::new();
        for lv in &urr.l_values {
            let aa = self.channel_radius(lv.awri, ay);
            let k = wave_number_factor() * lv.awri / (lv.awri + 1.0) * e.sqrt();
            let (rho, rhoc) = (k * aa, k * ay);
            let phi = phase_shift(lv.l, rhoc);
            let ab = 4.0 * PI / (k * k);
            if !lv.j_values.is_empty() {
                potential += ab * f64::from(2 * lv.l.number() + 1) * phi.sin().powi(2);
            }
            // ドップラー幅 Δ = 2 √(E kT / A)
            let delta = 2.0 * (e * BK * temperature / lv.awri).sqrt();
            for jv in &lv.j_values {
                let p = parameters(&jv.points, e);
                let gj = (2.0 * jv.aj + 1.0) / (4.0 * urr.spi + 2.0);
                // 平均中性子幅（換算幅 × 自由度 × 透過因子/ρ × √E。`cross_sections` と同じ）
                let gn =
                    p.gno * f64::from(jv.amun.number()) * penetrability(lv.l, rho) / rho * e.sqrt();
                sequences.push(LadderSequence {
                    scale: ab * gj,
                    gn,
                    gf: p.gf,
                    gg: p.gg,
                    gx: p.gx,
                    dof_n: jv.amun,
                    dof_f: if p.gf > 0.0 { jv.amuf } else { WidthDof::Fixed },
                    dof_x: if p.gx > 0.0 { jv.amux } else { WidthDof::Fixed },
                    d: p.d,
                    delta,
                    cos2: (2.0 * phi).cos(),
                    sin2p: (2.0 * phi).sin(),
                });
            }
        }
        bondarenko(&sequences, potential, background, sigma0, sampling, seed).map_err(|error| {
            match error {
                LadderError::EmptySampling => UrrError::EmptySampling {
                    ladders: sampling.ladders,
                    samples: sampling.samples,
                },
                LadderError::NonPositiveTotal { offset, floored } => UrrError::NonPositiveTotal {
                    energy: e,
                    offset,
                    floored,
                },
            }
        })
    }

    /// 位相シフトに使う散乱半径（NRO=1 なら表から補間する）。
    fn scattering_radius(&self, e: f64) -> Result<f64, UrrError> {
        scattering_radius_at(self.urr.ap, self.urr.ap_energy.as_ref(), e)
            .ok_or(UrrError::RadiusTable { energy: e })
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
fn fluctuation_integrals(
    gn: f64,
    gf: f64,
    gg: f64,
    gx: f64,
    jv: &UrrJValue,
) -> Result<[f64; 3], NumericError> {
    if gn <= 0.0 || gg <= 0.0 {
        return Ok([0.0; 3]);
    }
    let dof_f = if gf > 0.0 { jv.amuf } else { WidthDof::Fixed };
    let dof_x = if gx > 0.0 { jv.amux } else { WidthDof::Fixed };
    let ratios = |y: f64, z: f64| width_ratios(jv.amun, (gf * y + gg + gx * z) / gn);
    let outer = |g: &dyn Fn(f64, f64) -> Result<f64, NumericError>| {
        let value = expectation(dof_x, TOLERANCE[1], |z| {
            expectation(dof_f, TOLERANCE[0], |y| g(y, z))
        })?;
        Ok::<_, NumericError>(value / gn)
    };
    Ok([
        outer(&|y, z| Ok(ratios(y, z)?[1]))?,
        outer(&|y, z| Ok(ratios(y, z)?[0]))?,
        if gf > 0.0 {
            outer(&|y, z| Ok(y * ratios(y, z)?[0]))?
        } else {
            0.0
        },
    ])
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
mod tests;
