//! MF6（生成粒子のエネルギー・角度分布）。散乱行列に使う中性子だけを取り出す。
//!
//! 中性子の分布は LAW=1（連続エネルギー・角度分布）のエネルギー分布だけを扱い、角度は無視する
//! （M0 の近似）。重心系の分布は、群に振り分けるときに重心系で等方として実験室系に変換する。
//! LAW=6（N 体の位相空間分布）は重心系で等方な解析式の分布として持ち、LAW=7（実験室系の角度と
//! エネルギーの表）は角度で積分したエネルギー分布だけを持つ（角度のデータは使わない。ADR-0082）。どちらも ENDF-102 の
//! 定義のとおりで、NJOY GROUPR の表の作り方は真似ない（ADR-0271）。
//! LAW=2（2体反応の角度分布）と LAW=3（等方の2体反応）は、出射エネルギーが角度と Q 値で決まる
//! 2体反応であることだけを記録し、MF4 の離散準位と同じく重心系で等方の2体反応として扱う
//! （角度分布は使わない）。
//! 中性子以外の生成粒子は読み飛ばす。

use super::Error;
use super::mf5::{OutgoingDistributions, check_unit_base, incident_laws};
use super::tab1::{Interpolation, Tab1};
use super::tape::{Reader, Section};

/// 出射エネルギーの座標系（LCT）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Frame {
    Laboratory,
    CenterOfMass,
}

/// 1反応で放出される中性子の多重度と出射エネルギー分布。
#[derive(Debug, Clone, PartialEq)]
pub struct NeutronEmission {
    /// 多重度 y(E)
    pub multiplicity: Tab1,
    pub spectrum: NeutronSpectrum,
    pub frame: Frame,
}

/// 中性子の出射エネルギー分布。
#[derive(Debug, Clone, PartialEq)]
pub enum NeutronSpectrum {
    /// 入射エネルギーごとの表（LAW=1 の b0 と、LAW=7 を角度で積分したもの）
    Tabulated(OutgoingDistributions),
    /// LAW=6 の N 体の位相空間分布（いつも重心系）
    PhaseSpace(PhaseSpace),
}

/// LAW=6 の N 体の位相空間分布（ENDF-102 の 6.2.8）。重心系で等方で、重心系の出射エネルギー E' の密度は
/// C_n √E' (E_max − E')^(3n/2 − 4)。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhaseSpace {
    /// この分布で扱う n 個の粒子の質量の和（中性子の質量を単位。APSX）
    pub total_mass: f64,
    /// 粒子の数 n（NPSX。3〜5）
    pub particles: u8,
}

impl PhaseSpace {
    /// 入射エネルギー e (eV) での重心系の出射エネルギーの上限 E_max = (A_p − 1)/A_p (A/(A + 1) E + Q)。
    /// `awr` は標的の質量比 A、`q` は反応の Q 値 (eV)。しきい値の下では 0 以下になる。
    pub fn max_energy(&self, e: f64, awr: f64, q: f64) -> f64 {
        (self.total_mass - 1.0) / self.total_mass * (awr / (awr + 1.0) * e + q)
    }

    /// 密度の (1 − E'/E_max) の指数 3n/2 − 4
    pub fn exponent(&self) -> f64 {
        1.5 * f64::from(self.particles) - 4.0
    }

    /// u = E'/E_max の密度 √u (1 − u)^(3n/2 − 4) を 0〜1 で 1 にする係数 1/B(3/2, 3n/2 − 3)（ENDF-102 の C_n に
    /// E_max^(3n/2 − 2) と、余弦の範囲の幅 2 を掛けたもの）
    pub fn normalization(&self) -> f64 {
        use std::f64::consts::PI;
        // B(3/2, 3/2) = π/8、B(3/2, 3) = 16/105、B(3/2, 9/2) = 7π/256。`parse` が n を 3〜5 に限る
        match self.particles {
            3 => 8.0 / PI,
            4 => 105.0 / 16.0,
            _ => 256.0 / (7.0 * PI),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Mf6 {
    /// 標的の質量比（重心系から実験室系への変換に使う）
    pub awr: f64,
    /// 連続分布（LAW=1）で与えられた中性子
    pub neutrons: Vec<NeutronEmission>,
    /// 中性子が2体反応（LAW=2, 3）で与えられている
    pub two_body: bool,
}

impl Mf6 {
    /// 中性子に LAW=6 の位相空間分布がある（群に振り分けるのに反応の Q 値が要る）
    pub fn has_phase_space(&self) -> bool {
        self.neutrons
            .iter()
            .any(|n| matches!(n.spectrum, NeutronSpectrum::PhaseSpace(_)))
    }
}

impl Mf6 {
    pub fn parse(section: &Section<'_>) -> Result<Self, Error> {
        let mut r = section.reader();
        let head = r.cont()?;
        // LCT=3 は軽い粒子（中性子を含む）が重心系、重い粒子が実験室系
        let frame = match head.l2 {
            1 => Frame::Laboratory,
            2 | 3 => Frame::CenterOfMass,
            lct => return Err(Error::new(r.line(), format!("LCT={lct} is not supported"))),
        };
        let nk = usize::try_from(head.n1)
            .map_err(|_| Error::new(r.line(), format!("invalid NK {}", head.n1)))?;
        let mut neutrons = Vec::new();
        let mut two_body = false;
        for _ in 0..nk {
            let start = r.line() + 1;
            let product = r.tab1()?;
            let law = product.l2;
            if product.c1 != NEUTRON {
                skip_law(&mut r, law, start)?;
                continue;
            }
            match law {
                1 => neutrons.push(NeutronEmission {
                    multiplicity: product.table,
                    spectrum: NeutronSpectrum::Tabulated(read_law1(&mut r)?),
                    frame,
                }),
                // LAW=6 は LCT によらず重心系（ENDF-102 の 6.2.8）
                6 => {
                    xlc_log::warn!(
                        "neutron LAW=6 phase space is grouped at each incident energy; NJOY GROUPR \
                         tabulates it 1.2 apart in incident energy and interpolates the laboratory \
                         spectra by unit base (ADR-0271)"
                    );
                    neutrons.push(NeutronEmission {
                        multiplicity: product.table,
                        spectrum: NeutronSpectrum::PhaseSpace(read_law6(&mut r, head.c2)?),
                        frame: Frame::CenterOfMass,
                    });
                }
                // LAW=7 は実験室系の分布なので、重心系の LCT とは食い違う
                7 if frame != Frame::Laboratory => {
                    return Err(Error::new(
                        start,
                        format!("neutron LAW=7 (laboratory) with LCT={}", head.l2),
                    ));
                }
                7 => {
                    xlc_log::warn!(
                        "neutron LAW=7 is integrated over the cosines as tabulated; NJOY GROUPR's \
                         ll2lab thins the outgoing energies with swapped interpolation weights, \
                         which smears histogram tables (ADR-0271)"
                    );
                    neutrons.push(NeutronEmission {
                        multiplicity: product.table,
                        spectrum: NeutronSpectrum::Tabulated(read_law7(&mut r)?),
                        frame,
                    });
                }
                // LAW=2 は角度分布つき、LAW=3 は等方の2体反応（どちらも角度分布は使わない）
                2 | 3 => {
                    skip_law(&mut r, law, start)?;
                    two_body = true;
                }
                _ => {
                    return Err(Error::new(
                        start,
                        format!("neutron LAW={law} is not supported"),
                    ));
                }
            }
        }
        if two_body && !neutrons.is_empty() {
            return Err(Error::new(
                r.line(),
                "neutrons with both continuum (LAW=1, 6, 7) and two-body (LAW=2, 3) laws",
            ));
        }
        Ok(Self {
            awr: head.c2,
            neutrons,
            two_body,
        })
    }
}

/// 生成粒子の ZAP（中性子は 1）
const NEUTRON: f64 = 1.0;

/// LAW=1 のエネルギー分布（b0）を読む。離散線（ND≠0）は扱わない。
fn read_law1(r: &mut Reader<'_>) -> Result<OutgoingDistributions, Error> {
    let tab2 = r.tab2()?;
    let lep = match tab2.l2 {
        1 => Interpolation::Histogram,
        2 => Interpolation::LinLin,
        other => {
            return Err(Error::new(
                r.line(),
                format!("LEP={other} is not supported"),
            ));
        }
    };
    let line = r.line();
    let mut incident = Vec::new();
    let mut spectra = Vec::new();
    for _ in 0..tab2.nz {
        let list = r.list()?;
        let (nd, na, nep) = (list.l1, list.l2, list.n2);
        if nd != 0 {
            return Err(Error::new(
                r.line(),
                format!("ND={nd} discrete lines are not supported"),
            ));
        }
        let stride = usize::try_from(na + 2)
            .map_err(|_| Error::new(r.line(), format!("invalid NA {na}")))?;
        if usize::try_from(nep).ok().map(|n| n * stride) != Some(list.values.len()) {
            return Err(Error::new(r.line(), "NW does not match NEP and NA"));
        }
        let x: Vec<f64> = list.values.iter().step_by(stride).copied().collect();
        let y: Vec<f64> = list
            .values
            .iter()
            .skip(1)
            .step_by(stride)
            .copied()
            .collect();
        let table = Tab1 {
            breaks: vec![x.len()],
            laws: vec![lep],
            x,
            y,
        };
        table
            .check()
            .map_err(|message| Error::new(r.line(), message))?;
        incident.push(list.c2);
        spectra.push(table);
    }
    let counts: Vec<usize> = spectra.iter().map(|t| t.x.len()).collect();
    let laws = incident_laws(&tab2.codes, &tab2.breaks, &counts, line)?;
    check_unit_base(&laws, &tab2.breaks, &spectra, line)?;
    Ok(OutgoingDistributions {
        breaks: tab2.breaks,
        laws,
        incident,
        spectra,
    })
}

/// LAW=6 の CONT（APSX、NPSX）を読む。`awr` は標的の質量比。
fn read_law6(r: &mut Reader<'_>, awr: f64) -> Result<PhaseSpace, Error> {
    let cont = r.cont()?;
    let particles = match cont.n2 {
        3 => 3,
        4 => 4,
        5 => 5,
        n => {
            return Err(Error::new(
                r.line(),
                format!(
                    "phase space of NPSX={n} particles is not supported (ENDF-102 gives 3 to 5)"
                ),
            ));
        }
    };
    if !(cont.c1.is_finite() && cont.c1 > 1.0) {
        return Err(Error::new(
            r.line(),
            format!(
                "phase-space mass APSX={} must exceed the neutron mass",
                cont.c1
            ),
        ));
    }
    // APSX が生成粒子の質量の和（およそ A + 1）より小さいのは 2 段の反応の 2 段目で、E_a は 1 段目の反跳のエネルギー
    // （ENDF-102 の 6.2.8）。式 6.28 の 1 段の E_a は使えない。1 段の H-2 と Li-6 の APSX は A + 1 との差が Q/(m_n c²) ほど
    // （0.002〜0.004）なので、中性子の質量の半分を境にする
    if cont.c1 < awr + 1.0 - 0.5 {
        return Err(Error::new(
            r.line(),
            format!(
                "phase-space mass APSX={} below the products' mass {} is a two-step reaction, \
                 which is not supported",
                cont.c1,
                awr + 1.0
            ),
        ));
    }
    Ok(PhaseSpace {
        total_mass: cont.c1,
        particles,
    })
}

/// LAW=7 の実験室系の f(μ, E') を、入射エネルギーごとに余弦で積分したエネルギー分布 g(E') = ∫ f(μ, E') dμ にする。
///
/// 余弦の間は、同じ E' での値を余弦の補間則（lin-lin かヒストグラム）で補間する（ENDF-102 の 6.2.9）。g は余弦の
/// 重み w_i と f(μ_i, E') の和で、どの表も E' の lin-lin かヒストグラムなので、全部の表の E' の区切りの間で線形になる。
/// 区切りでは左と右の値を並べた lin-lin の表にする（ヒストグラムの段は同じ E' の2点になる）。
fn read_law7(r: &mut Reader<'_>) -> Result<OutgoingDistributions, Error> {
    let tab2 = r.tab2()?;
    let line = r.line();
    let mut incident = Vec::new();
    let mut spectra = Vec::new();
    for _ in 0..tab2.nz {
        let start = r.line() + 1;
        let cosines = r.tab2()?;
        let mut mu = Vec::with_capacity(cosines.nz);
        let mut tables = Vec::with_capacity(cosines.nz);
        for _ in 0..cosines.nz {
            let record = r.tab1()?;
            let table = record.table;
            if let Some(k) = (0..table.x.len().saturating_sub(1)).find(|&k| {
                !matches!(
                    table.law_of_interval(k),
                    Interpolation::LinLin | Interpolation::Histogram
                )
            }) {
                return Err(Error::new(
                    r.line(),
                    format!(
                        "LAW=7 outgoing-energy interpolation {:?} is not supported",
                        table.law_of_interval(k)
                    ),
                ));
            }
            mu.push(record.c2);
            tables.push(table);
        }
        let weights = cosine_weights(&mu, &cosines.breaks, &cosines.codes)
            .map_err(|message| Error::new(start, message))?;
        let table = integrate_cosines(&weights, &tables);
        if table.x.len() < 2 {
            return Err(Error::new(
                start,
                "LAW=7 distribution has fewer than 2 outgoing energies",
            ));
        }
        table
            .check()
            .map_err(|message| Error::new(start, message))?;
        incident.push(cosines.c2);
        spectra.push(table);
    }
    let counts: Vec<usize> = spectra.iter().map(|t| t.x.len()).collect();
    let laws = incident_laws(&tab2.codes, &tab2.breaks, &counts, line)?;
    check_unit_base(&laws, &tab2.breaks, &spectra, line)?;
    Ok(OutgoingDistributions {
        breaks: tab2.breaks,
        laws,
        incident,
        spectra,
    })
}

/// 余弦 μ_i の表を余弦の補間則で補間して −1〜1 で積分するときの、各表の重み（lin-lin は台形、ヒストグラムは左の表）。
/// 余弦は昇順で −1 から 1 までであること。
fn cosine_weights(mu: &[f64], breaks: &[usize], codes: &[i64]) -> Result<Vec<f64>, String> {
    if mu
        .windows(2)
        .any(|w| w[0].partial_cmp(&w[1]).is_none_or(|o| o.is_ge()))
    {
        return Err(format!("LAW=7 cosines {mu:?} must ascend"));
    }
    // −1〜1 を覆わない分布は、積分すると和が 1 に足りない
    if mu.first() != Some(&-1.0) || mu.last() != Some(&1.0) {
        return Err(format!("LAW=7 cosines {mu:?} must run from -1 to 1"));
    }
    let mut weights = vec![0.0; mu.len()];
    for k in 0..mu.len().saturating_sub(1) {
        let width = mu[k + 1] - mu[k];
        match codes[breaks.partition_point(|&b| b < k + 2)] {
            1 => weights[k] += width,
            2 => {
                weights[k] += 0.5 * width;
                weights[k + 1] += 0.5 * width;
            }
            code => {
                return Err(format!(
                    "LAW=7 cosine interpolation INT={code} is not supported"
                ));
            }
        }
    }
    Ok(weights)
}

/// 重み `weights` の表 `tables`（E' の lin-lin かヒストグラム）の和を、lin-lin の表にする。
fn integrate_cosines(weights: &[f64], tables: &[Tab1]) -> Tab1 {
    let mut grid: Vec<f64> = tables.iter().flat_map(|t| t.x.iter().copied()).collect();
    grid.sort_by(f64::total_cmp);
    grid.dedup();
    let last = grid.len().saturating_sub(1);
    let (mut x, mut y) = (Vec::new(), Vec::new());
    for (k, &e) in grid.iter().enumerate() {
        let (mut left, mut right) = (0.0, 0.0);
        for (w, t) in weights.iter().zip(tables) {
            let (l, r) = limits(t, e);
            left += w * l;
            right += w * r;
        }
        // 範囲の始めの左と終わりの右は表の外（値 0）なので、点にしない
        if k > 0 {
            x.push(e);
            y.push(left);
        }
        if k < last && (k == 0 || right != left) {
            x.push(e);
            y.push(right);
        }
    }
    Tab1::linear(x, y)
}

/// E' の lin-lin かヒストグラムの表 `t` の、e での左と右の極限。表の範囲の外は 0。
fn limits(t: &Tab1, e: f64) -> (f64, f64) {
    let n = t.x.len();
    if n == 0 || e < t.x[0] || e > t.x[n - 1] {
        return (0.0, 0.0);
    }
    // e 以上の最初の点と、e 以下の最後の点
    let first = t.x.partition_point(|&x| x < e);
    let last = t.x.partition_point(|&x| x <= e).saturating_sub(1);
    let value = |k: usize| {
        let law = t.law_of_interval(k);
        law.interpolate(e, (t.x[k], t.y[k]), (t.x[k + 1], t.y[k + 1]))
    };
    let left = if first == 0 {
        0.0
    } else if t.x[first] == e {
        // 区間 (first − 1, first) の右端。ヒストグラムは左の点の値
        match t.law_of_interval(first - 1) {
            Interpolation::Histogram => t.y[first - 1],
            _ => t.y[first],
        }
    } else {
        value(first - 1)
    };
    let right = if last == n - 1 { 0.0 } else { value(last) };
    (left, right)
}

/// 中性子以外の生成粒子の分布を、LAW ごとのレコード構成に従って読み飛ばす。
fn skip_law(r: &mut Reader<'_>, law: i64, start: usize) -> Result<(), Error> {
    match law {
        0 | 3 | 4 => {}
        1 | 2 | 5 => {
            let tab2 = r.tab2()?;
            for _ in 0..tab2.nz {
                r.list()?;
            }
        }
        6 => {
            r.cont()?;
        }
        7 => {
            let tab2 = r.tab2()?;
            for _ in 0..tab2.nz {
                let inner = r.tab2()?;
                for _ in 0..inner.nz {
                    r.tab1()?;
                }
            }
        }
        _ => return Err(Error::new(start, format!("invalid LAW={law}"))),
    }
    Ok(())
}

#[cfg(test)]
mod tests;
