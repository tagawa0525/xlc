//! MF5（出射中性子のエネルギー分布）。表形式（LF=1）と、θ(E) で決まる解析式のうち
//! 単純な核分裂スペクトル（LF=7）と蒸発スペクトル（LF=9）を扱う。それ以外の形式は、
//! 対象核種で必要になるまで追加しない。

use super::Error;
use super::tab1::{ClampedTab1, Interpolation, Tab1};
use super::tape::{Reader, Section};

/// 入射エネルギーごとに与えられた、出射エネルギーの確率密度の表。
#[derive(Debug, Clone, PartialEq)]
pub struct OutgoingDistributions {
    /// 入射エネルギー間の補間則の範囲（NBT）
    pub breaks: Vec<usize>,
    /// 入射エネルギー間の補間則
    pub laws: Vec<IncidentLaw>,
    /// 入射エネルギー (eV)
    pub incident: Vec<f64>,
    /// 各入射エネルギーでの出射エネルギーの確率密度 (1/eV)
    pub spectra: Vec<Tab1>,
}

/// 部分分布の形。解析式の分布は、出射エネルギー E' を 0 ≤ E' ≤ E - U に制限する。
#[derive(Debug, Clone, PartialEq)]
pub enum EnergyLaw {
    /// LF=1: 入射エネルギーごとの表
    Tabulated(OutgoingDistributions),
    /// LF=7: 単純な核分裂スペクトル f(E') ∝ √E' exp(-E'/θ(E))
    Maxwellian {
        restriction: f64,
        theta: ClampedTab1,
    },
    /// LF=9: 蒸発スペクトル f(E') ∝ E' exp(-E'/θ(E))
    Evaporation {
        restriction: f64,
        theta: ClampedTab1,
    },
}

/// MF5 の1反応。部分分布の確率 p_k(E) と、その分布の組。
#[derive(Debug, Clone, PartialEq)]
pub struct Mf5 {
    pub subsections: Vec<(Tab1, EnergyLaw)>,
}

impl Mf5 {
    pub fn parse(section: &Section<'_>) -> Result<Self, Error> {
        let mut r = section.reader();
        let head = r.cont()?;
        let nk = usize::try_from(head.n1)
            .map_err(|_| Error::new(r.line(), format!("invalid NK {}", head.n1)))?;
        let mut subsections = Vec::new();
        for _ in 0..nk {
            let start = r.line() + 1;
            let probability = r.tab1()?;
            // 解析式の分布では、p(E) の TAB1 の C1 が U で、θ(E) の TAB1 が続く
            let restriction = probability.c1;
            let law = match probability.l2 {
                1 => EnergyLaw::Tabulated(parse_outgoing(&mut r)?),
                7 => EnergyLaw::Maxwellian {
                    restriction,
                    theta: positive_theta(&mut r)?,
                },
                9 => EnergyLaw::Evaporation {
                    restriction,
                    theta: positive_theta(&mut r)?,
                },
                lf => return Err(Error::new(start, format!("LF={lf} is not supported"))),
            };
            subsections.push((probability.table, law));
        }
        Ok(Self { subsections })
    }
}

/// 解析式の分布の θ(E) を読む。θ は分布の温度なので正でなければならない。
fn positive_theta(r: &mut Reader<'_>) -> Result<ClampedTab1, Error> {
    let theta = r.tab1()?.table;
    if let Some(bad) = theta.y.iter().find(|&&t| t <= 0.0 || t.is_nan()) {
        return Err(Error::new(
            r.line(),
            format!("theta must be positive, got {bad}"),
        ));
    }
    ClampedTab1::new(theta).ok_or_else(|| Error::new(r.line(), "theta table is empty"))
}

/// TAB2 と、それに続く NE 個の TAB1（C2 が入射エネルギー）を読む。
/// MF5 の LF=1 と MF6 の出射分布で共通の形。
pub fn parse_outgoing(r: &mut Reader<'_>) -> Result<OutgoingDistributions, Error> {
    let tab2 = r.tab2()?;
    let line = r.line();
    let mut incident = Vec::new();
    let mut spectra = Vec::new();
    for _ in 0..tab2.nz {
        let record = r.tab1()?;
        incident.push(record.c2);
        spectra.push(record.table);
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

/// 入射エネルギー間の出射分布の補間則（ENDF の INT）。
///
/// INT=1〜5 はふつうの補間、21〜25 は単位基底（出射エネルギーの範囲を 0〜1 に正規化して補間する）。
/// 11〜15 は対応点の補間で、隣り合う分布の点の数が違うと定義できないので、NJOY と同じく単位基底に
/// 切り替える。点の数が同じ対応点の補間は扱わない。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IncidentLaw {
    /// 補間の式（単位基底では正規化した座標で使う）
    pub law: Interpolation,
    pub unit_base: bool,
}

/// INT の値の方式（十の位）
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Scheme {
    Cartesian,
    CorrespondingPoint,
    UnitBase,
}

/// 入射エネルギー間の補間則を、出射分布の点の数を見て決める。`counts` は各入射エネルギーの分布の点の数。
pub(super) fn incident_laws(
    codes: &[i64],
    breaks: &[usize],
    counts: &[usize],
    line: usize,
) -> Result<Vec<IncidentLaw>, Error> {
    let unsupported = |code: i64| {
        Error::new(
            line,
            format!("incident-energy interpolation INT={code} is not supported"),
        )
    };
    codes
        .iter()
        .enumerate()
        .map(|(range, &code)| {
            let scheme = match code / 10 {
                0 => Scheme::Cartesian,
                1 => Scheme::CorrespondingPoint,
                2 => Scheme::UnitBase,
                _ => return Err(unsupported(code)),
            };
            let law = Interpolation::from_code(code % 10).ok_or_else(|| unsupported(code))?;
            // この補間則の範囲の区間（点 k と k+1 の間）
            let first = if range == 0 { 0 } else { breaks[range - 1] - 1 };
            let intervals = first..breaks[range].saturating_sub(1);
            let unit_base = match scheme {
                Scheme::Cartesian => false,
                Scheme::UnitBase => true,
                Scheme::CorrespondingPoint => {
                    if intervals.clone().all(|k| counts[k] == counts[k + 1]) {
                        return Err(Error::new(
                            line,
                            format!(
                                "corresponding-point interpolation (INT={code}) with equal grids \
                                 is not supported"
                            ),
                        ));
                    }
                    tracing::warn!(
                        code,
                        "corresponding-point interpolation between grids of different sizes; \
                         using unit-base interpolation as NJOY does"
                    );
                    true
                }
            };
            Ok(IncidentLaw { law, unit_base })
        })
        .collect()
}

/// 単位基底で補間する分布は、2点以上で、出射エネルギーの範囲の幅が正であること（範囲で正規化するため）。
pub(super) fn check_unit_base(
    laws: &[IncidentLaw],
    breaks: &[usize],
    spectra: &[Tab1],
    line: usize,
) -> Result<(), Error> {
    for (range, law) in laws.iter().enumerate() {
        if !law.unit_base {
            continue;
        }
        let first = if range == 0 { 0 } else { breaks[range - 1] - 1 };
        for spectrum in &spectra[first..breaks[range].min(spectra.len())] {
            let width = match spectrum.x[..] {
                [low, .., high] => high - low,
                _ => 0.0,
            };
            if width <= 0.0 || width.is_nan() {
                return Err(Error::new(
                    line,
                    "unit-base interpolation needs outgoing spectra of positive width",
                ));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::endf::{Tape, rec};

    /// 入射エネルギー 1e-5 と 2e7 eV で、出射分布が [lo, hi] の一様分布の MF5。
    pub(crate) fn two_point_mf5(
        lf: &str,
        first: (&str, &str, &str),
        second: (&str, &str, &str),
    ) -> String {
        let r = |f: [&str; 6]| rec(f, 1, 5, 18);
        [
            r(["9.223500+4", "2.330248+2", "0", "0", "1", "0"]),
            r(["0.0", "0.0", "0", lf, "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["1.0-5", "1.0", "2.0+7", "1.0", "", ""]),
            r(["0.0", "0.0", "0", "0", "1", "2"]),
            // 入射エネルギーの間はふつうの線形補間
            r(["2", "2", "", "", "", ""]),
            r(["0.0", "1.0-5", "0", "0", "1", "2"]),
            r(["2", "1", "", "", "", ""]),
            r([first.0, first.2, first.1, "0.0", "", ""]),
            r(["0.0", "2.0+7", "0", "0", "1", "2"]),
            r(["2", "1", "", "", "", ""]),
            r([second.0, second.2, second.1, "0.0", "", ""]),
        ]
        .concat()
    }

    #[test]
    fn reads_tabulated_distribution() {
        let text = two_point_mf5(
            "1",
            ("1.0+6", "2.0+7", "5.263158-8"),
            ("1.0-5", "1.0+6", "1.0-6"),
        );
        let tape = Tape::parse(&text).unwrap();
        let mf5 = Mf5::parse(&tape.section(1, 5, 18).unwrap()).unwrap();
        assert_eq!(mf5.subsections.len(), 1);
        let (p, EnergyLaw::Tabulated(d)) = &mf5.subsections[0] else {
            panic!("{:?}", mf5.subsections[0].1);
        };
        assert_eq!(p.eval(1.0), Some(1.0));
        assert_eq!(d.incident, vec![1.0e-5, 2.0e7]);
        assert_eq!(
            d.laws,
            vec![IncidentLaw {
                law: Interpolation::LinLin,
                unit_base: false
            }]
        );
        assert_eq!(d.spectra[0].x, vec![1.0e6, 2.0e7]);
        assert_eq!(d.spectra[1].y, vec![1.0e-6, 0.0]);
    }

    /// p(E) = 1 で、θ(E) が 2 点の解析式の MF5。U は p(E) の TAB1 の C1。
    fn analytic_mf5(lf: &str, restriction: &str) -> String {
        let r = |f: [&str; 6]| rec(f, 1, 5, 16);
        [
            r(["8.017000+3", "1.683600+1", "0", "0", "1", "0"]),
            r([restriction, "0.0", "0", lf, "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["4.0+6", "1.0", "2.0+7", "1.0", "", ""]),
            r(["0.0", "0.0", "0", "0", "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["4.0+6", "2.0+5", "2.0+7", "6.0+5", "", ""]),
        ]
        .concat()
    }

    #[test]
    fn reads_analytic_spectra() {
        for (lf, restriction) in [("7", "-3.0+7"), ("9", "1.0+6")] {
            let tape = Tape::parse(&analytic_mf5(lf, restriction)).unwrap();
            let mf5 = Mf5::parse(&tape.section(1, 5, 16).unwrap()).unwrap();
            let [(p, law)] = &mf5.subsections[..] else {
                panic!("{:?}", mf5.subsections);
            };
            assert_eq!(p.eval(1.0e7), Some(1.0));
            let (EnergyLaw::Maxwellian { restriction, theta }
            | EnergyLaw::Evaporation { restriction, theta }) = law
            else {
                panic!("{law:?}");
            };
            assert_eq!(theta.eval(1.2e7), 4.0e5);
            match lf {
                "7" => {
                    assert!(matches!(law, EnergyLaw::Maxwellian { .. }) && *restriction == -3.0e7)
                }
                _ => assert!(matches!(law, EnergyLaw::Evaporation { .. }) && *restriction == 1.0e6),
            }
        }
    }

    #[test]
    fn theta_must_be_positive() {
        let text = analytic_mf5("9", "1.0+6").replace("2.0+5", "0.0+0");
        let tape = Tape::parse(&text).unwrap();
        assert!(Mf5::parse(&tape.section(1, 5, 16).unwrap()).is_err());
    }

    #[test]
    fn other_lf_is_error() {
        let text = two_point_mf5("11", ("1.0+6", "2.0+7", "1.0"), ("1.0-5", "1.0+6", "1.0"));
        let tape = Tape::parse(&text).unwrap();
        let err = Mf5::parse(&tape.section(1, 5, 18).unwrap()).unwrap_err();
        assert_eq!(err.line, 2);
    }
}
