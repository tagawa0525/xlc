//! MF4（出射中性子の角度分布）。等方（LTT=0）、Legendre 係数（LTT=1）、表（LTT=2）と、低いエネルギーで
//! Legendre 係数・高いエネルギーで表を使う混合（LTT=3）を扱う。

use super::Error;
use super::mf6::Frame;
use super::tab1::{Interpolation, Tab1};
use super::tape::{Reader, Section};

/// 入射エネルギーごとの Legendre 係数 a_l（l = 1..NL）。f(μ) = Σ (2l+1)/2 a_l P_l(μ)、a_0 = 1。
#[derive(Debug, Clone, PartialEq)]
pub struct LegendreTable {
    /// 入射エネルギー間の補間則の範囲（NBT）と補間則
    pub breaks: Vec<usize>,
    pub laws: Vec<Interpolation>,
    /// 入射エネルギー (eV)
    pub energies: Vec<f64>,
    /// [入射エネルギー][l − 1]
    pub coefficients: Vec<Vec<f64>>,
}

/// 入射エネルギーごとの散乱角の余弦 μ の確率密度 f(μ) の表。
#[derive(Debug, Clone, PartialEq)]
pub struct AngleTables {
    pub breaks: Vec<usize>,
    pub laws: Vec<Interpolation>,
    pub energies: Vec<f64>,
    pub tables: Vec<Tab1>,
}

/// 角度分布の形。
#[derive(Debug, Clone, PartialEq)]
pub enum Angular {
    /// LTT=0（LI=1）: 等方
    Isotropic,
    /// LTT=1
    Legendre(LegendreTable),
    /// LTT=2
    Tabulated(AngleTables),
    /// LTT=3: `legendre` の最後の入射エネルギーまでは Legendre 係数、その上は表
    Mixed {
        legendre: LegendreTable,
        tabulated: AngleTables,
    },
}

/// MF4 の1反応。
#[derive(Debug, Clone, PartialEq)]
pub struct Mf4 {
    pub frame: Frame,
    pub angular: Angular,
}

impl Mf4 {
    pub fn parse(section: &Section<'_>) -> Result<Self, Error> {
        let _ = section;
        Err(Error::new(0, "MF4 is not supported yet"))
    }

    /// 入射エネルギー `e` での角度分布の Legendre 係数 a_0..a_order（a_0 = 1。`frame` の系の μ について）。
    pub fn moments(&self, e: f64, order: usize) -> Vec<f64> {
        let _ = (e, order);
        Vec::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::endf::Tape;
    use crate::endf::tape::tests::rec;

    fn r(f: [&str; 6]) -> String {
        rec(f, 825, 4, 2)
    }

    fn parse(text: &str) -> Result<Mf4, Error> {
        let tape = Tape::parse(text).unwrap();
        Mf4::parse(&tape.section(825, 4, 2).unwrap())
    }

    fn close(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len(), "{actual:?}");
        for (a, e) in actual.iter().zip(expected) {
            assert!((a - e).abs() < 1e-12, "{actual:?} vs {expected:?}");
        }
    }

    /// Legendre 係数（LTT=1）の部分: 入射エネルギー 1e5 と 2e5 eV で a_1 = 0.1 と 0.3、a_2 = 0.02 と 0.04
    fn legendre_part() -> String {
        [
            r(["0.0", "0.0", "0", "0", "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["0.0", "1.0e5", "0", "0", "2", "0"]),
            r(["0.1", "0.02", "", "", "", ""]),
            r(["0.0", "2.0e5", "0", "0", "2", "0"]),
            r(["0.3", "0.04", "", "", "", ""]),
        ]
        .concat()
    }

    /// 表（LTT=2）の部分: 2e5 eV で f(μ) = 0.5 + 0.3 μ、4e5 eV で等方
    fn tabulated_part() -> String {
        [
            r(["0.0", "0.0", "0", "0", "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["0.0", "2.0e5", "0", "0", "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["-1.0", "0.2", "1.0", "0.8", "", ""]),
            r(["0.0", "4.0e5", "0", "0", "1", "2"]),
            r(["2", "2", "", "", "", ""]),
            r(["-1.0", "0.5", "1.0", "0.5", "", ""]),
        ]
        .concat()
    }

    fn head(ltt: &str, li: &str, lct: &str) -> String {
        [
            r(["8016.0", "15.85751", "0", ltt, "0", "0"]),
            r(["0.0", "15.85751", li, lct, "0", "0"]),
        ]
        .concat()
    }

    #[test]
    fn isotropic_distribution() {
        let mf4 = parse(&head("0", "1", "1")).unwrap();
        assert_eq!(mf4.frame, Frame::Laboratory);
        assert_eq!(mf4.angular, Angular::Isotropic);
        close(&mf4.moments(1.0e6, 2), &[1.0, 0.0, 0.0]);
    }

    #[test]
    fn legendre_coefficients_are_interpolated_between_energies() {
        let mf4 = parse(&(head("1", "0", "2") + &legendre_part())).unwrap();
        assert_eq!(mf4.frame, Frame::CenterOfMass);
        close(&mf4.moments(1.5e5, 3), &[1.0, 0.2, 0.03, 0.0]);
        // 次数を打ち切る
        close(&mf4.moments(1.0e5, 1), &[1.0, 0.1]);
    }

    #[test]
    fn tabulated_distribution_gives_its_moments() {
        // f(μ) = 0.5 + 0.3 μ の a_1 = ∫ f μ dμ = 0.2、a_2 = 0
        let mf4 = parse(&(head("2", "0", "2") + &tabulated_part())).unwrap();
        close(&mf4.moments(2.0e5, 2), &[1.0, 0.2, 0.0]);
        close(&mf4.moments(3.0e5, 2), &[1.0, 0.1, 0.0]);
    }

    #[test]
    fn mixed_distribution_switches_at_the_last_legendre_energy() {
        let text = head("3", "0", "2") + &legendre_part() + &tabulated_part();
        let mf4 = parse(&text).unwrap();
        close(&mf4.moments(1.5e5, 2), &[1.0, 0.2, 0.03]);
        close(&mf4.moments(3.0e5, 2), &[1.0, 0.1, 0.0]);
    }

    #[test]
    fn transformation_matrix_is_not_supported() {
        // NM > 0（Legendre 係数の変換行列）は使われていないので扱わない
        let text = [
            r(["8016.0", "15.85751", "0", "1", "0", "1"]),
            r(["0.0", "15.85751", "0", "2", "0", "1"]),
        ]
        .concat();
        assert!(parse(&text).is_err());
    }
}
