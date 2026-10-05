//! MF6（生成粒子のエネルギー・角度分布）。散乱行列に使う中性子だけを取り出す。
//!
//! 中性子の分布は LAW=1（連続エネルギー・角度分布）のエネルギー分布だけを扱い、角度は無視する
//! （M0 の近似）。重心系の分布は、群に振り分けるときに重心系で等方として実験室系に変換する。
//! 中性子以外の生成粒子は読み飛ばす。

use super::Error;
use super::mf5::{OutgoingDistributions, incident_laws};
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
    pub neutrons: Vec<NeutronEmission>,
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
        for _ in 0..nk {
            let start = r.line() + 1;
            let product = r.tab1()?;
            let law = product.l2;
            if product.c1 == NEUTRON {
                if law != 1 {
                    return Err(Error::new(
                        start,
                        format!("neutron LAW={law} is not supported"),
                    ));
                }
                neutrons.push(NeutronEmission {
                    multiplicity: product.table,
                    distributions: read_law1(&mut r)?,
                    frame,
                });
            } else {
                skip_law(&mut r, law, start)?;
            }
        }
        Ok(Self {
            awr: head.c2,
            neutrons,
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
    let laws = incident_laws(&tab2.codes, r.line())?;
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
mod tests {
    use super::*;
    use crate::endf::{Tape, rec};

    fn r(f: [&str; 6]) -> String {
        rec(f, 825, 6, 16)
    }

    /// 中性子（LAW=1、2つの入射エネルギー）と、LAW=0 の別の生成粒子
    fn section(lct: &str, neutron_law: &str, nd: &str) -> String {
        [
            r(["8016.0", "15.85751", "0", lct, "2", "0"]),
            // 別の生成粒子（LAW=0: 多重度だけ）
            r(["8015.0", "14.87", "0", "0", "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["1.6e7", "1.0", "2.0e7", "1.0", "", ""]),
            // 中性子
            r(["1.0", "1.0", "0", neutron_law, "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["1.6e7", "2.0", "2.0e7", "2.0", "", ""]),
            r(["0.0", "0.0", "1", "2", "1", "2"]),
            r(["2", "22", "", "", "", ""]),
            r(["0.0", "1.6e7", nd, "0", "4", "2"]),
            r(["0.0", "2.0e-6", "1.0e6", "0.0", "", ""]),
            r(["0.0", "2.0e7", "0", "1", "9", "3"]),
            // (E', b0, b1) の3つ組が3点
            r(["0.0", "0.0", "0.0", "1.0e6", "2.0e-6", "0.5"]),
            r(["2.0e6", "0.0", "0.0", "", "", ""]),
        ]
        .concat()
    }

    #[test]
    fn reads_law1_neutrons_and_skips_other_products() {
        let tape = Tape::parse(&section("2", "1", "0")).unwrap();
        let mf6 = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap();
        assert_eq!(mf6.neutrons.len(), 1);
        assert_eq!(mf6.awr, 15.85751);
        let n = &mf6.neutrons[0];
        assert_eq!(n.frame, Frame::CenterOfMass);
        assert_eq!(n.multiplicity.eval(1.8e7), Some(2.0));
        let d = &n.distributions;
        assert_eq!(d.incident, vec![1.6e7, 2.0e7]);
        assert_eq!(d.laws, vec![Interpolation::LinLin]);
        assert_eq!(d.spectra[0].x, vec![0.0, 1.0e6]);
        assert_eq!(d.spectra[0].y, vec![2.0e-6, 0.0]);
        // NA=1 なので1点あたり3個の値。b0 だけを取る
        assert_eq!(d.spectra[1].x, vec![0.0, 1.0e6, 2.0e6]);
        assert_eq!(d.spectra[1].y, vec![0.0, 2.0e-6, 0.0]);
    }

    #[test]
    fn other_neutron_law_is_error() {
        let tape = Tape::parse(&section("2", "7", "0")).unwrap();
        let err = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap_err();
        assert_eq!(err.line, 5);
    }

    #[test]
    fn laboratory_frame() {
        let tape = Tape::parse(&section("1", "1", "0")).unwrap();
        let mf6 = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap();
        assert_eq!(mf6.neutrons[0].frame, Frame::Laboratory);
    }

    #[test]
    fn unsupported_frame_is_error() {
        // LCT=4 は生成粒子ごとに座標系が変わる（M0 では扱わない）
        let tape = Tape::parse(&section("4", "1", "0")).unwrap();
        let err = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap_err();
        assert_eq!(err.line, 1);
    }

    #[test]
    fn discrete_lines_are_error() {
        let tape = Tape::parse(&section("2", "1", "1")).unwrap();
        assert!(Mf6::parse(&tape.section(825, 6, 16).unwrap()).is_err());
    }
}
