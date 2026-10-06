//! MF2 MT151（共鳴パラメータ）セクション。
//!
//! M0 では、分離共鳴の Reich-Moore 形式（LRU=1, LRF=3）と、無限希釈の断面積が MF3 にある
//! 非分離共鳴（LRU=2）だけを扱う。それ以外の形式はエラーにする。

use super::tape::{Reader, Section};
use super::{Error, Interpolation};

/// 1核種（同位体1つ）の共鳴パラメータ。
#[derive(Debug, Clone, PartialEq)]
pub struct Mf2 {
    pub za: f64,
    pub awr: f64,
    pub ranges: Vec<ResonanceRange>,
}

/// エネルギー範囲 [el, eh] (eV) の共鳴データ。
#[derive(Debug, Clone, PartialEq)]
pub struct ResonanceRange {
    pub el: f64,
    pub eh: f64,
    pub kind: RangeKind,
}

#[derive(Debug, Clone, PartialEq)]
pub enum RangeKind {
    /// LRU=0: 散乱半径だけ（共鳴なし）
    ScatteringRadius { ap: f64 },
    /// LRU=1, LRF=3
    ReichMoore(ReichMoore),
    /// LRU=1, LRF=2
    Mlbw(Mlbw),
    /// LRU=2, LRF=2
    Unresolved(Unresolved),
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReichMoore {
    /// 標的核のスピン
    pub spi: f64,
    /// 散乱半径 (10^-12 cm)
    pub ap: f64,
    /// 透過因子に使うチャネル半径の選び方（NAPS）
    pub naps: i64,
    pub l_values: Vec<RmLValue>,
}

/// 多準位 Breit-Wigner（MLBW）の共鳴パラメータ。
#[derive(Debug, Clone, PartialEq)]
pub struct Mlbw {
    pub spi: f64,
    /// 散乱半径 (10^-12 cm)
    pub ap: f64,
    pub naps: i64,
    pub l_values: Vec<BwLValue>,
}

/// MLBW の1つの l に属する共鳴。競合反応の幅（LRX≠0）は扱わない。
#[derive(Debug, Clone, PartialEq)]
pub struct BwLValue {
    pub awri: f64,
    pub l: u32,
    pub resonances: Vec<BwResonance>,
}

/// MLBW の共鳴1つ。幅は共鳴エネルギーでの値 (eV)。全幅は Γn + Γγ + Γf とする（GT は使わない）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BwResonance {
    pub er: f64,
    pub aj: f64,
    pub gn: f64,
    pub gg: f64,
    pub gf: f64,
}

/// 非分離共鳴のエネルギー依存のパラメータ（LRU=2, LRF=2）。
#[derive(Debug, Clone, PartialEq)]
pub struct Unresolved {
    /// 自己遮蔽の扱い（LSSF）。1 なら無限希釈の断面積は MF3 にある
    pub lssf: i64,
    pub spi: f64,
    /// 散乱半径 (10^-12 cm)
    pub ap: f64,
    pub naps: i64,
    pub l_values: Vec<UrrLValue>,
}

/// 非分離共鳴の1つの l。
#[derive(Debug, Clone, PartialEq)]
pub struct UrrLValue {
    pub awri: f64,
    pub l: u32,
    pub j_values: Vec<UrrJValue>,
}

/// 非分離共鳴の1つの (l, J) の系列。`int` はパラメータの、エネルギーについての補間則。
/// 幅の分布の自由度は、競合反応 `amux`・中性子 `amun`・核分裂 `amuf`（γ線の自由度 AMUG は
/// 幅のゆらぎを無視するので読まない）。
#[derive(Debug, Clone, PartialEq)]
pub struct UrrJValue {
    pub aj: f64,
    pub int: Interpolation,
    pub amux: f64,
    pub amun: f64,
    pub amuf: f64,
    pub points: Vec<UrrPoint>,
}

/// エネルギー `e` (eV) での平均準位間隔 `d` と平均幅（競合 `gx`、換算中性子幅 `gno`、
/// 捕獲 `gg`、核分裂 `gf`）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UrrPoint {
    pub e: f64,
    pub d: f64,
    pub gx: f64,
    pub gno: f64,
    pub gg: f64,
    pub gf: f64,
}

/// 1つの軌道角運動量 l に属する共鳴。
#[derive(Debug, Clone, PartialEq)]
pub struct RmLValue {
    pub awri: f64,
    /// l ごとの散乱半径。0 なら `ReichMoore::ap` を使う
    pub apl: f64,
    pub l: u32,
    pub resonances: Vec<RmResonance>,
}

/// Reich-Moore の共鳴1つ。幅はすべて共鳴エネルギーでの値 (eV)。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RmResonance {
    pub er: f64,
    /// 全角運動量 J。符号はチャネルスピンを表す
    pub aj: f64,
    pub gn: f64,
    pub gg: f64,
    pub gfa: f64,
    pub gfb: f64,
}

impl Mf2 {
    pub fn parse(section: &Section<'_>) -> Result<Self, Error> {
        let mut r = section.reader();
        let head = r.cont()?;
        if head.n1 != 1 {
            return Err(Error::new(
                r.line(),
                format!("NIS={} isotopes are not supported", head.n1),
            ));
        }
        let isotope = r.cont()?;
        let ner = usize::try_from(isotope.n1)
            .map_err(|_| Error::new(r.line(), format!("invalid NER {}", isotope.n1)))?;

        // NER はファイルの値なので、それに合わせて先に確保しない
        let mut ranges = Vec::new();
        for _ in 0..ner {
            let range = r.cont()?;
            let range_line = r.line();
            let (lru, lrf, nro) = (range.l1, range.l2, range.n1);
            let unsupported = |what: String| Err(Error::new(range_line, what));
            let kind = match lru {
                0 => {
                    let cont = r.cont()?;
                    RangeKind::ScatteringRadius { ap: cont.c2 }
                }
                1 if nro != 0 => {
                    return unsupported(
                        "energy-dependent radius (NRO≠0) in the resolved range".into(),
                    );
                }
                1 if lrf == 3 => RangeKind::ReichMoore(read_reich_moore(&mut r, range.n2)?),
                1 if lrf == 2 => RangeKind::Mlbw(read_mlbw(&mut r, range.n2)?),
                1 => return unsupported(format!("resolved LRF={lrf} is not supported")),
                2 => RangeKind::Unresolved(read_unresolved(&mut r, range.n2)?),
                _ => return unsupported(format!("invalid LRU={lru}")),
            };
            ranges.push(ResonanceRange {
                el: range.c1,
                eh: range.c2,
                kind,
            });
        }
        Ok(Self {
            za: head.c1,
            awr: head.c2,
            ranges,
        })
    }
}

fn read_unresolved(r: &mut Reader<'_>, naps: i64) -> Result<Unresolved, Error> {
    let _ = (r, naps);
    todo!()
}

fn read_reich_moore(r: &mut Reader<'_>, naps: i64) -> Result<ReichMoore, Error> {
    let cont = r.cont()?;
    let nls = usize::try_from(cont.n1)
        .map_err(|_| Error::new(r.line(), format!("invalid NLS {}", cont.n1)))?;
    let mut l_values = Vec::with_capacity(nls.min(MAX_L + 1));
    for _ in 0..nls {
        let list = r.list()?;
        let l = u32::try_from(list.l1)
            .ok()
            .filter(|&l| l as usize <= MAX_L)
            .ok_or_else(|| Error::new(r.line(), format!("unsupported L={}", list.l1)))?;
        // LIST の値は1行に1共鳴（6個）ずつ並ぶ
        let first_value_line = r.line() + 1 - list.values.len().div_ceil(6);
        if let Some(i) = list.values.iter().step_by(6).position(|&er| er == 0.0) {
            return Err(Error::new(
                first_value_line + i,
                "resonance at zero energy (penetrability is zero)",
            ));
        }
        let resonances = list
            .values
            .as_chunks::<6>()
            .0
            .iter()
            .map(|p| RmResonance {
                er: p[0],
                aj: p[1],
                gn: p[2],
                gg: p[3],
                gfa: p[4],
                gfb: p[5],
            })
            .collect();
        l_values.push(RmLValue {
            awri: list.c1,
            apl: list.c2,
            l,
            resonances,
        });
    }
    Ok(ReichMoore {
        spi: cont.c1,
        ap: cont.c2,
        naps,
        l_values,
    })
}

fn read_mlbw(r: &mut Reader<'_>, naps: i64) -> Result<Mlbw, Error> {
    let cont = r.cont()?;
    let spi = cont.c1;
    let nls = usize::try_from(cont.n1)
        .map_err(|_| Error::new(r.line(), format!("invalid NLS {}", cont.n1)))?;
    let mut l_values = Vec::with_capacity(nls.min(MAX_L + 1));
    for _ in 0..nls {
        let header_line = r.line() + 1;
        let list = r.list()?;
        let l = u32::try_from(list.l1)
            .ok()
            .filter(|&l| l as usize <= MAX_L)
            .ok_or_else(|| Error::new(header_line, format!("unsupported L={}", list.l1)))?;
        if list.l2 != 0 {
            return Err(Error::new(
                header_line,
                format!("competitive width (LRX={}) is not supported", list.l2),
            ));
        }
        // 標的スピンと l から作れる J
        let fl = f64::from(l);
        let (ajmin, ajmax) = (((spi - fl).abs() - 0.5).abs(), spi + fl + 0.5);
        let mut resonances = Vec::new();
        for (i, p) in list.values.as_chunks::<6>().0.iter().enumerate() {
            let line = header_line + 1 + i;
            if p[0] == 0.0 {
                return Err(Error::new(
                    line,
                    "resonance at zero energy (penetrability is zero)",
                ));
            }
            let steps = p[1] - ajmin;
            if p[1] > ajmax + 1e-6 || (steps - steps.round()).abs() > 1e-6 || steps < -1e-6 {
                return Err(Error::new(
                    line,
                    format!("J={} is not allowed for L={l}", p[1]),
                ));
            }
            resonances.push(BwResonance {
                er: p[0],
                aj: p[1],
                gn: p[3],
                gg: p[4],
                gf: p[5],
            });
        }
        l_values.push(BwLValue {
            awri: list.c1,
            l,
            resonances,
        });
    }
    Ok(Mlbw {
        spi,
        ap: cont.c2,
        naps,
        l_values,
    })
}

/// 透過因子・位相シフトを実装している l の上限（NJOY と同じ）。
pub(crate) const MAX_L: usize = 4;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::endf::tape::tests::rec;
    use crate::endf::{Interpolation, Tape};

    fn mf2(lines: &[[&str; 6]]) -> Result<Mf2, Error> {
        let text: String = lines.iter().map(|f| rec(*f, 9237, 2, 151)).collect();
        let tape = Tape::parse(&text).unwrap();
        Mf2::parse(&tape.section(9237, 2, 151).unwrap())
    }

    const HEAD: [&str; 6] = ["9.223800+4", "2.360058+2", "0", "0", "1", "0"];

    #[test]
    fn reads_reich_moore_and_unresolved_ranges() {
        let data = mf2(&[
            HEAD,
            ["9.223800+4", "1.0", "0", "0", "2", "0"],
            // RRR: EL EH LRU LRF NRO NAPS
            ["1.0-5", "2.0+4", "1", "3", "0", "0"],
            ["0.0", "0.0", "0", "0", "1", "1"],
            // l=0: AWRI APL L 0 6*NRS NRS
            ["2.360060+2", "9.480000-1", "0", "0", "12", "2"],
            ["-4.0", "0.5", "2.0-2", "2.3-2", "0.0", "1.0-6"],
            ["6.674", "-0.5", "1.5-3", "2.3-2", "0.0", "0.0"],
            // URR: SPI AP LSSF 0 NLS 0
            ["2.0+4", "1.49+5", "2", "2", "0", "0"],
            ["0.0", "0.948", "1", "0", "1", "0"],
            // l=0: AWRI 0 L 0 NJS 0
            ["2.360060+2", "0.0", "0", "0", "1", "0"],
            // J: AJ 0 INT 0 6*NE+6 NE、続いて 0 0 AMUX AMUN AMUG AMUF と ES D GX GNO GG GF
            ["0.5", "0.0", "5", "0", "18", "2"],
            ["0.0", "0.0", "0.0", "1.0", "0.0", "0.0"],
            ["2.0+4", "20.0", "0.0", "2.0-3", "2.3-2", "0.0"],
            ["1.49+5", "18.0", "0.0", "1.8-3", "2.3-2", "0.0"],
        ])
        .unwrap();
        assert_eq!(data.za, 92238.0);
        assert_eq!(data.ranges.len(), 2);
        assert_eq!((data.ranges[0].el, data.ranges[0].eh), (1.0e-5, 2.0e4));
        let RangeKind::ReichMoore(rm) = &data.ranges[0].kind else {
            panic!("expected Reich-Moore");
        };
        assert_eq!((rm.spi, rm.ap, rm.naps), (0.0, 0.0, 0));
        assert_eq!(rm.l_values.len(), 1);
        assert_eq!(rm.l_values[0].apl, 0.948);
        assert_eq!(
            rm.l_values[0].resonances[1],
            RmResonance {
                er: 6.674,
                aj: -0.5,
                gn: 1.5e-3,
                gg: 2.3e-2,
                gfa: 0.0,
                gfb: 0.0
            }
        );
        let RangeKind::Unresolved(urr) = &data.ranges[1].kind else {
            panic!("expected unresolved");
        };
        assert_eq!((urr.lssf, urr.spi, urr.ap, urr.naps), (1, 0.0, 0.948, 0));
        assert_eq!(urr.l_values.len(), 1);
        assert_eq!((urr.l_values[0].awri, urr.l_values[0].l), (236.006, 0));
        let j = &urr.l_values[0].j_values[0];
        assert_eq!((j.aj, j.int), (0.5, Interpolation::LogLog));
        assert_eq!((j.amux, j.amun, j.amuf), (0.0, 1.0, 0.0));
        assert_eq!(
            j.points[1],
            UrrPoint {
                e: 1.49e5,
                d: 18.0,
                gx: 0.0,
                gno: 1.8e-3,
                gg: 2.3e-2,
                gf: 0.0
            }
        );
    }

    #[test]
    fn rejects_unresolved_lrf1() {
        let result = mf2(&[
            HEAD,
            ["9.223800+4", "1.0", "0", "0", "1", "0"],
            ["2.0+4", "1.49+5", "2", "1", "0", "0"],
            ["0.0", "0.948", "1", "0", "1", "0"],
        ]);
        assert!(result.is_err());
    }

    #[test]
    fn reads_mlbw_range() {
        let data = mf2(&[
            ["4.009000+4", "8.913240+1", "0", "0", "1", "0"],
            ["4.009000+4", "1.0", "0", "0", "1", "0"],
            ["1.0-5", "2.0+5", "1", "2", "0", "0"],
            ["0.0", "7.16-1", "0", "0", "1", "0"],
            // l=1: AWRI QX L LRX 6*NRS NRS、値は ER AJ GT GN GG GF
            ["8.913240+1", "0.0", "1", "0", "6", "1"],
            ["2.4e3", "1.5", "0.25", "0.2", "0.05", "0.0"],
        ])
        .unwrap();
        let RangeKind::Mlbw(bw) = &data.ranges[0].kind else {
            panic!("expected MLBW");
        };
        assert_eq!((bw.spi, bw.ap, bw.naps), (0.0, 0.716, 0));
        assert_eq!(bw.l_values[0].l, 1);
        assert_eq!(
            bw.l_values[0].resonances[0],
            BwResonance {
                er: 2.4e3,
                aj: 1.5,
                gn: 0.2,
                gg: 0.05,
                gf: 0.0
            }
        );
    }

    #[test]
    fn competitive_width_is_error() {
        let err = mf2(&[
            ["4.009000+4", "8.913240+1", "0", "0", "1", "0"],
            ["4.009000+4", "1.0", "0", "0", "1", "0"],
            ["1.0-5", "2.0+5", "1", "2", "0", "0"],
            ["0.0", "7.16-1", "0", "0", "1", "0"],
            ["8.913240+1", "-1.0+5", "1", "1", "6", "1"],
            ["2.4e3", "1.5", "0.25", "0.2", "0.05", "0.0"],
        ])
        .unwrap_err();
        assert_eq!(err.line, 5);
    }

    #[test]
    fn reads_scattering_radius_only_range() {
        let data = mf2(&[
            ["1.001000+3", "9.991673-1", "0", "0", "1", "0"],
            ["1.001000+3", "1.0", "0", "0", "1", "0"],
            ["1.0-5", "1.0+5", "0", "0", "0", "0"],
            ["0.5", "1.276553+0", "0", "0", "0", "0"],
        ])
        .unwrap();
        assert_eq!(
            data.ranges[0].kind,
            RangeKind::ScatteringRadius { ap: 1.276553 }
        );
    }

    #[test]
    fn unsupported_formalism_is_error() {
        // LRF=7 (R-matrix limited)
        let err = mf2(&[
            HEAD,
            ["9.223800+4", "1.0", "0", "0", "1", "0"],
            ["1.0-5", "2.0+4", "1", "7", "0", "0"],
        ])
        .unwrap_err();
        assert_eq!(err.line, 3);
    }

    #[test]
    fn energy_dependent_radius_in_reich_moore_is_error() {
        let err = mf2(&[
            HEAD,
            ["9.223800+4", "1.0", "0", "0", "1", "0"],
            ["1.0-5", "2.0+4", "1", "3", "1", "0"],
        ])
        .unwrap_err();
        assert_eq!(err.line, 3);
    }

    #[test]
    fn zero_resonance_energy_is_error() {
        // Er=0 では透過因子が 0 になり、中性子幅を換算できない
        let err = mf2(&[
            HEAD,
            ["9.223800+4", "1.0", "0", "0", "1", "0"],
            ["1.0-5", "2.0+4", "1", "3", "0", "0"],
            ["0.0", "0.0", "0", "0", "1", "1"],
            ["2.360060+2", "9.480000-1", "0", "0", "6", "1"],
            ["0.0", "0.5", "2.0-2", "2.3-2", "0.0", "0.0"],
        ])
        .unwrap_err();
        assert_eq!(err.line, 6);
    }

    #[test]
    fn huge_range_count_is_error_not_abort() {
        let err = mf2(&[HEAD, ["9.223800+4", "1.0", "0", "0", "99999999999", "0"]]).unwrap_err();
        assert_eq!(err.line, 3);
    }

    #[test]
    fn multiple_isotopes_are_error() {
        let err = mf2(&[["9.223800+4", "2.360058+2", "0", "0", "2", "0"]]).unwrap_err();
        assert_eq!(err.line, 1);
    }
}
