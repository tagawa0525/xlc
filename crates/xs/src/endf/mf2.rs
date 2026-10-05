//! MF2 MT151（共鳴パラメータ）セクション。
//!
//! M0 では、分離共鳴の Reich-Moore 形式（LRU=1, LRF=3）と、無限希釈の断面積が MF3 にある
//! 非分離共鳴（LRU=2）だけを扱う。それ以外の形式はエラーにする。

use super::Error;
use super::tape::{Reader, Section};

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
    /// LRU=2。M0 では LSSF=1（断面積は MF3 にある）のときだけ使える
    Unresolved { lssf: i64 },
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

        let mut ranges = Vec::with_capacity(ner);
        for index in 0..ner {
            let range = r.cont()?;
            let range_line = r.line();
            let (lru, lrf, nro) = (range.l1, range.l2, range.n1);
            let unsupported = |what: String| Err(Error::new(range_line, what));
            let kind = match lru {
                0 => {
                    let cont = r.cont()?;
                    RangeKind::ScatteringRadius { ap: cont.c2 }
                }
                1 if lrf != 3 => {
                    return unsupported(format!("resolved LRF={lrf} is not supported"));
                }
                1 if nro != 0 => {
                    return unsupported("energy-dependent radius (NRO≠0) with Reich-Moore".into());
                }
                1 => RangeKind::ReichMoore(read_reich_moore(&mut r, range.n2)?),
                2 => {
                    // 非分離共鳴は LSSF だけを読む。続くパラメータは読まないので、最後の範囲に限る
                    if index + 1 != ner {
                        return unsupported("unresolved range must be the last range".into());
                    }
                    if nro != 0 {
                        r.tab1()?;
                    }
                    RangeKind::Unresolved { lssf: r.cont()?.l1 }
                }
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

/// 透過因子・位相シフトを実装している l の上限（NJOY と同じ）。
pub(crate) const MAX_L: usize = 4;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::endf::Tape;
    use crate::endf::tape::tests::rec;

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
            // URR
            ["2.0+4", "1.49+5", "2", "2", "0", "0"],
            ["0.0", "0.948", "1", "0", "2", "0"],
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
        assert_eq!(data.ranges[1].kind, RangeKind::Unresolved { lssf: 1 });
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
        // LRF=2 (MLBW)
        let err = mf2(&[
            HEAD,
            ["9.223800+4", "1.0", "0", "0", "1", "0"],
            ["1.0-5", "2.0+4", "1", "2", "0", "0"],
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

    #[test]
    fn unresolved_range_must_be_last() {
        let err = mf2(&[
            HEAD,
            ["9.223800+4", "1.0", "0", "0", "2", "0"],
            ["2.0+4", "1.49+5", "2", "2", "0", "0"],
            ["0.0", "0.948", "1", "0", "2", "0"],
            ["1.49+5", "2.0+5", "0", "0", "0", "0"],
        ])
        .unwrap_err();
        assert_eq!(err.line, 3);
    }
}
