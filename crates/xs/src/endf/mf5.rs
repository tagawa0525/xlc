//! MF5（出射中性子のエネルギー分布）。M0 では表形式（LF=1）だけを扱う。

use super::Error;
use super::tab1::{Interpolation, Tab1};
use super::tape::{Reader, Section};

/// 入射エネルギーごとに与えられた、出射エネルギーの確率密度の表。
#[derive(Debug, Clone, PartialEq)]
pub struct OutgoingDistributions {
    /// 入射エネルギー間の補間則の範囲（NBT）
    pub breaks: Vec<usize>,
    /// 入射エネルギー間の補間則。unit-base（INT=21, 22）は通常の補間として扱う（M0 の近似）
    pub laws: Vec<Interpolation>,
    /// 入射エネルギー (eV)
    pub incident: Vec<f64>,
    /// 各入射エネルギーでの出射エネルギーの確率密度 (1/eV)
    pub spectra: Vec<Tab1>,
}

/// MF5 の1反応。部分分布の確率 p_k(E) と、その分布の組。
#[derive(Debug, Clone, PartialEq)]
pub struct Mf5 {
    pub subsections: Vec<(Tab1, OutgoingDistributions)>,
}

impl Mf5 {
    pub fn parse(section: &Section<'_>) -> Result<Self, Error> {
        todo!()
    }
}

/// TAB2 と、それに続く NE 個の TAB1（C2 が入射エネルギー）を読む。
/// MF5 の LF=1 と MF6 の出射分布で共通の形。
pub fn parse_outgoing(r: &mut Reader<'_>) -> Result<OutgoingDistributions, Error> {
    todo!()
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
            r(["2", "22", "", "", "", ""]),
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
        let (p, d) = &mf5.subsections[0];
        assert_eq!(p.eval(1.0), Some(1.0));
        assert_eq!(d.incident, vec![1.0e-5, 2.0e7]);
        assert_eq!(d.laws, vec![Interpolation::LinLin]);
        assert_eq!(d.spectra[0].x, vec![1.0e6, 2.0e7]);
        assert_eq!(d.spectra[1].y, vec![1.0e-6, 0.0]);
    }

    #[test]
    fn other_lf_is_error() {
        let text = two_point_mf5("9", ("1.0+6", "2.0+7", "1.0"), ("1.0-5", "1.0+6", "1.0"));
        let tape = Tape::parse(&text).unwrap();
        let err = Mf5::parse(&tape.section(1, 5, 18).unwrap()).unwrap_err();
        assert_eq!(err.line, 2);
    }
}
