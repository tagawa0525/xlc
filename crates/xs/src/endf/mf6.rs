//! MF6（生成粒子のエネルギー・角度分布）。散乱行列に使う中性子だけを取り出す。
//!
//! 中性子の分布は LAW=1（連続エネルギー・角度分布）のエネルギー分布だけを扱い、角度は無視する
//! （M0 の近似）。重心系の分布は、群に振り分けるときに重心系で等方として実験室系に変換する。
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
    pub distributions: OutgoingDistributions,
    pub frame: Frame,
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
                    distributions: read_law1(&mut r)?,
                    frame,
                }),
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
                "neutrons with both continuum (LAW=1) and two-body (LAW=2, 3) laws",
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
