//! MF2 MT151（共鳴パラメータ）セクション。
//!
//! 分離共鳴の Reich-Moore 形式（LRU=1, LRF=3）と MLBW 形式（LRU=1, LRF=2）、
//! エネルギー依存のパラメータによる非分離共鳴（LRU=2, LRF=2）を扱う。それ以外の形式はエラーにする。

use super::tape::{Reader, Section};
use super::{Error, Interpolation, Tab1};

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
    /// LRU=2, LRF=1。パラメータは読み飛ばし、LSSF だけを持つ（無限希釈の断面積が MF3 にあれば使える）
    UnresolvedLrf1 { lssf: i64 },
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
    pub l: OrbitalL,
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
    /// エネルギー依存の散乱半径（NRO=1）。あれば位相シフトにはこちらを使う
    pub ap_energy: Option<Tab1>,
    pub naps: i64,
    pub l_values: Vec<UrrLValue>,
}

/// 非分離共鳴の1つの l。
#[derive(Debug, Clone, PartialEq)]
pub struct UrrLValue {
    pub awri: f64,
    pub l: OrbitalL,
    pub j_values: Vec<UrrJValue>,
}

/// 非分離共鳴の1つの (l, J) の系列。`int` はパラメータの、エネルギーについての補間則。
/// 幅の分布の自由度は、競合反応 `amux`・中性子 `amun`・核分裂 `amuf`（γ線の自由度 AMUG は
/// 幅のゆらぎを無視するので読まない）。
#[derive(Debug, Clone, PartialEq)]
pub struct UrrJValue {
    pub aj: f64,
    pub int: Interpolation,
    pub amux: u32,
    pub amun: u32,
    pub amuf: u32,
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
    pub l: OrbitalL,
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
        // 核分裂幅がエネルギー依存か（非分離共鳴 LRF=1 の形式を決める）
        let lfw = isotope.l2;
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
                2 if lrf == 1 || lrf == 2 => {
                    let ap_energy = match nro {
                        0 => None,
                        1 => Some(r.tab1()?.table),
                        _ => return unsupported(format!("invalid NRO={nro}")),
                    };
                    if lrf == 2 {
                        RangeKind::Unresolved(read_unresolved(&mut r, range.n2, ap_energy)?)
                    } else {
                        RangeKind::UnresolvedLrf1 {
                            lssf: skip_unresolved_lrf1(&mut r, lfw)?,
                        }
                    }
                }
                2 => return unsupported(format!("unresolved LRF={lrf} is not supported")),
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

fn read_unresolved(
    r: &mut Reader<'_>,
    naps: i64,
    ap_energy: Option<Tab1>,
) -> Result<Unresolved, Error> {
    let cont = r.cont()?;
    let count = |n: i64, line: usize, what: &str| {
        usize::try_from(n).map_err(|_| Error::new(line, format!("invalid {what} {n}")))
    };
    let nls = count(cont.n1, r.line(), "NLS")?;
    let mut l_values = Vec::with_capacity(nls.min(OrbitalL::COUNT));
    for _ in 0..nls {
        let l_cont = r.cont()?;
        let l = OrbitalL::from_endf(l_cont.l1)
            .ok_or_else(|| Error::new(r.line(), format!("unsupported L={}", l_cont.l1)))?;
        let njs = count(l_cont.n1, r.line(), "NJS")?;
        let mut j_values = Vec::with_capacity(njs.min(OrbitalL::COUNT));
        for _ in 0..njs {
            let header_line = r.line() + 1;
            let list = r.list()?;
            let int = Interpolation::from_code(list.l1).ok_or_else(|| {
                Error::new(
                    header_line,
                    format!("invalid interpolation INT={}", list.l1),
                )
            })?;
            let (chunks, rest) = list.values.as_chunks::<6>();
            let ([header, points @ ..], []) = (chunks, rest) else {
                return Err(Error::new(
                    header_line,
                    "unresolved parameters must be 6 values per energy after a header",
                ));
            };
            // 幅の分布の自由度は、NJOY と同じく最も近い整数に丸める。χ² 分布の求積があるのは 0〜4
            let dof = |v: f64| {
                let rounded = v.round();
                if (0.0..=4.0).contains(&rounded) {
                    Ok(rounded as u32)
                } else {
                    Err(Error::new(
                        header_line + 1,
                        format!("degrees of freedom {v} must round to an integer from 0 to 4"),
                    ))
                }
            };
            let points: Vec<UrrPoint> = points
                .iter()
                .map(|p| UrrPoint {
                    e: p[0],
                    d: p[1],
                    gx: p[2],
                    gno: p[3],
                    gg: p[4],
                    gf: p[5],
                })
                .collect();
            // 同じエネルギーの繰り返しは、パラメータの段差を表す
            if points.is_empty() || points.windows(2).any(|w| w[1].e < w[0].e) {
                return Err(Error::new(
                    header_line,
                    "unresolved energies must be present and not decreasing",
                ));
            }
            j_values.push(UrrJValue {
                aj: list.c1,
                int,
                amux: dof(header[2])?,
                amun: dof(header[3])?,
                amuf: dof(header[5])?,
                points,
            });
        }
        l_values.push(UrrLValue {
            awri: l_cont.c1,
            l,
            j_values,
        });
    }
    Ok(Unresolved {
        lssf: cont.l1,
        spi: cont.c1,
        ap: cont.c2,
        ap_energy,
        naps,
        l_values,
    })
}

/// 非分離共鳴 LRF=1 のパラメータを読み飛ばし、LSSF を返す。
fn skip_unresolved_lrf1(r: &mut Reader<'_>, lfw: i64) -> Result<i64, Error> {
    if lfw == 0 {
        // SPI AP LSSF 0 NLS 0 と、l ごとに J のパラメータの LIST
        let cont = r.cont()?;
        for _ in 0..cont.n1 {
            r.list()?;
        }
        Ok(cont.l1)
    } else {
        // エネルギーの LIST（SPI AP LSSF 0 NE NLS）と、l ごとに CONT（NJS）と J ごとの LIST
        let energies = r.list()?;
        for _ in 0..energies.n2 {
            let l_cont = r.cont()?;
            for _ in 0..l_cont.n1 {
                r.list()?;
            }
        }
        Ok(energies.l1)
    }
}

fn read_reich_moore(r: &mut Reader<'_>, naps: i64) -> Result<ReichMoore, Error> {
    let cont = r.cont()?;
    let nls = usize::try_from(cont.n1)
        .map_err(|_| Error::new(r.line(), format!("invalid NLS {}", cont.n1)))?;
    let mut l_values = Vec::with_capacity(nls.min(OrbitalL::COUNT));
    for _ in 0..nls {
        let list = r.list()?;
        let l = OrbitalL::from_endf(list.l1)
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
    let mut l_values = Vec::with_capacity(nls.min(OrbitalL::COUNT));
    for _ in 0..nls {
        let header_line = r.line() + 1;
        let list = r.list()?;
        let l = OrbitalL::from_endf(list.l1)
            .ok_or_else(|| Error::new(header_line, format!("unsupported L={}", list.l1)))?;
        if list.l2 != 0 {
            return Err(Error::new(
                header_line,
                format!("competitive width (LRX={}) is not supported", list.l2),
            ));
        }
        // 標的スピンと l から作れる J
        let fl = f64::from(l.number());
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

/// 軌道角運動量の量子数 L。透過因子・位相シフトを実装している 0〜4 だけを扱う（NJOY と同じ）。
/// MF2 を読むときにこの範囲を確かめ、型で持つので、共鳴の計算で範囲外を考えなくてよい。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrbitalL {
    S,
    P,
    D,
    F,
    G,
}

impl OrbitalL {
    /// 扱える L の数
    const COUNT: usize = 5;

    /// ENDF の L の値から作る（範囲外なら `None`）。
    fn from_endf(l: i64) -> Option<Self> {
        match l {
            0 => Some(Self::S),
            1 => Some(Self::P),
            2 => Some(Self::D),
            3 => Some(Self::F),
            4 => Some(Self::G),
            _ => None,
        }
    }

    pub fn number(self) -> u32 {
        match self {
            Self::S => 0,
            Self::P => 1,
            Self::D => 2,
            Self::F => 3,
            Self::G => 4,
        }
    }
}

impl std::fmt::Display for OrbitalL {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.number())
    }
}

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
        assert_eq!(urr.ap_energy, None);
        assert_eq!(urr.l_values.len(), 1);
        assert_eq!(
            (urr.l_values[0].awri, urr.l_values[0].l),
            (236.006, OrbitalL::S)
        );
        let j = &urr.l_values[0].j_values[0];
        assert_eq!((j.aj, j.int), (0.5, Interpolation::LogLog));
        assert_eq!((j.amux, j.amun, j.amuf), (0, 1, 0));
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
    fn reads_energy_dependent_radius_of_unresolved_range() {
        let data = mf2(&[
            HEAD,
            ["9.223800+4", "1.0", "0", "0", "1", "0"],
            ["2.0+4", "1.49+5", "2", "2", "1", "0"],
            // AP(E) の TAB1: 0 0 0 0 NR NP、補間の区切り、(E, AP) の組
            ["0.0", "0.0", "0", "0", "1", "2"],
            ["2", "2", "", "", "", ""],
            ["2.0+4", "0.95", "1.49+5", "0.90", "", ""],
            ["0.0", "0.948", "1", "0", "1", "0"],
            ["2.360060+2", "0.0", "0", "0", "1", "0"],
            ["0.5", "0.0", "2", "0", "12", "1"],
            ["0.0", "0.0", "0.0", "1.0", "0.0", "0.0"],
            ["2.0+4", "20.0", "0.0", "2.0-3", "2.3-2", "0.0"],
        ])
        .unwrap();
        let RangeKind::Unresolved(urr) = &data.ranges[0].kind else {
            panic!("expected unresolved");
        };
        let radius = urr.ap_energy.as_ref().unwrap();
        assert_eq!(
            (radius.x.clone(), radius.y.clone()),
            (vec![2.0e4, 1.49e5], vec![0.95, 0.90])
        );
        assert_eq!(urr.ap, 0.948);
    }

    /// 非分離共鳴 LRF=1 の範囲のあとに LRU=0 の範囲を置いて、パラメータを読み飛ばせたことを確かめる。
    fn assert_skips_lrf1(lfw: &str, parameters: &[[&str; 6]]) {
        assert_skips_lrf1_with_radius(lfw, "0", parameters);
    }

    fn assert_skips_lrf1_with_radius(lfw: &str, nro: &str, parameters: &[[&str; 6]]) {
        let mut lines = vec![
            HEAD,
            ["9.423900+4", "1.0", "0", lfw, "2", "0"],
            ["5.7+3", "4.0+4", "2", "1", nro, "0"],
        ];
        lines.extend_from_slice(parameters);
        lines.push(["4.0+4", "2.0+7", "0", "0", "0", "0"]);
        lines.push(["0.0", "0.5", "0", "0", "0", "0"]);
        let data = mf2(&lines).unwrap();
        assert_eq!(data.ranges[0].kind, RangeKind::UnresolvedLrf1 { lssf: 1 });
        assert_eq!(data.ranges[1].kind, RangeKind::ScatteringRadius { ap: 0.5 });
    }

    #[test]
    fn skips_unresolved_lrf1_without_fission_widths() {
        assert_skips_lrf1(
            "0",
            &[
                // SPI AP LSSF 0 NLS 0、l ごとに AWRI 0 L 0 6*NJS NJS と D AJ AMUN GNO GG GF
                ["0.0", "0.888", "1", "0", "1", "0"],
                ["2.379916+2", "0.0", "0", "0", "6", "1"],
                ["13.1", "0.5", "1.0", "1.5-3", "3.1-2", "0.0"],
            ],
        );
    }

    #[test]
    fn skips_unresolved_lrf1_with_fission_widths() {
        assert_skips_lrf1(
            "1",
            &[
                // SPI AP LSSF 0 NE NLS と ES、l ごとに AWRI 0 L 0 NJS 0、
                // J ごとに 0 0 L MUF NE+6 0 と D AJ AMUN GNO GG 0 GF...
                ["0.0", "0.888", "1", "0", "2", "1"],
                ["5.7+3", "4.0+4", "", "", "", ""],
                ["2.379916+2", "0.0", "0", "0", "1", "0"],
                ["0.0", "0.0", "0", "1", "8", "0"],
                ["13.1", "0.5", "1.0", "1.5-3", "3.1-2", "0.0"],
                ["1.2-3", "1.4-3", "", "", "", ""],
            ],
        );
    }

    #[test]
    fn skips_unresolved_lrf1_with_energy_dependent_radius() {
        assert_skips_lrf1_with_radius(
            "0",
            "1",
            &[
                // AP(E) の TAB1 が先に来る
                ["0.0", "0.0", "0", "0", "1", "2"],
                ["2", "2", "", "", "", ""],
                ["5.7+3", "0.89", "4.0+4", "0.88", "", ""],
                ["0.0", "0.888", "1", "0", "1", "0"],
                ["2.379916+2", "0.0", "0", "0", "6", "1"],
                ["13.1", "0.5", "1.0", "1.5-3", "3.1-2", "0.0"],
            ],
        );
    }

    #[test]
    fn rounds_degrees_of_freedom_like_njoy() {
        // NJOY は自由度を最も近い整数に丸めて使う（Kr-83 の AMUX=1.0123 など）
        let data = mf2(&[
            HEAD,
            ["9.223800+4", "1.0", "0", "0", "1", "0"],
            ["2.0+4", "1.49+5", "2", "2", "0", "0"],
            ["0.0", "0.948", "1", "0", "1", "0"],
            ["2.360060+2", "0.0", "0", "0", "1", "0"],
            ["0.5", "0.0", "2", "0", "18", "2"],
            ["0.0", "0.0", "1.0123", "0.9", "0.0", "2.4"],
            // 同じエネルギーの繰り返しはパラメータの段差を表す（Th-229 など）
            ["2.0+4", "20.0", "0.0", "2.0-3", "2.3-2", "0.0"],
            ["2.0+4", "18.0", "0.0", "1.8-3", "2.3-2", "0.0"],
        ])
        .unwrap();
        let RangeKind::Unresolved(urr) = &data.ranges[0].kind else {
            panic!("expected unresolved");
        };
        let j = &urr.l_values[0].j_values[0];
        assert_eq!((j.amux, j.amun, j.amuf), (1, 1, 2));
        assert_eq!(j.points.len(), 2);
    }

    #[test]
    fn rejects_degrees_of_freedom_beyond_quadrature() {
        let result = mf2(&[
            HEAD,
            ["9.223800+4", "1.0", "0", "0", "1", "0"],
            ["2.0+4", "1.49+5", "2", "2", "0", "0"],
            ["0.0", "0.948", "1", "0", "1", "0"],
            ["2.360060+2", "0.0", "0", "0", "1", "0"],
            ["0.5", "0.0", "2", "0", "12", "1"],
            ["0.0", "0.0", "0.0", "5.0", "0.0", "0.0"],
            ["2.0+4", "20.0", "0.0", "2.0-3", "2.3-2", "0.0"],
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
        assert_eq!(bw.l_values[0].l, OrbitalL::P);
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
