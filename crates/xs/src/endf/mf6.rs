//! MF6（生成粒子のエネルギー・角度分布）。散乱行列に使う中性子だけを取り出す。
//!
//! 中性子の分布は LAW=1（連続エネルギー・角度分布）のエネルギー分布だけを扱い、角度は無視する
//! （M0 の近似）。中性子以外の生成粒子は読み飛ばす。

use super::Error;
use super::mf5::OutgoingDistributions;
use super::tab1::Tab1;
use super::tape::Section;

/// 1反応で放出される中性子の多重度と出射エネルギー分布。
#[derive(Debug, Clone, PartialEq)]
pub struct NeutronEmission {
    /// 多重度 y(E)
    pub multiplicity: Tab1,
    pub distributions: OutgoingDistributions,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Mf6 {
    pub neutrons: Vec<NeutronEmission>,
}

impl Mf6 {
    pub fn parse(section: &Section<'_>) -> Result<Self, Error> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::endf::{Interpolation, Tape, rec};

    fn r(f: [&str; 6]) -> String {
        rec(f, 825, 6, 16)
    }

    /// 中性子（LAW=1、2つの入射エネルギー）と、LAW=0 の別の生成粒子
    fn section(neutron_law: &str, nd: &str) -> String {
        [
            r(["8016.0", "15.85751", "0", "2", "2", "0"]),
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
        let tape = Tape::parse(&section("1", "0")).unwrap();
        let mf6 = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap();
        assert_eq!(mf6.neutrons.len(), 1);
        let n = &mf6.neutrons[0];
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
        let tape = Tape::parse(&section("7", "0")).unwrap();
        let err = Mf6::parse(&tape.section(825, 6, 16).unwrap()).unwrap_err();
        assert_eq!(err.line, 5);
    }

    #[test]
    fn discrete_lines_are_error() {
        let tape = Tape::parse(&section("1", "1")).unwrap();
        assert!(Mf6::parse(&tape.section(825, 6, 16).unwrap()).is_err());
    }
}
