//! MF7 MT4（熱中性子の非干渉性非弾性散乱の S(α, β)）。
//!
//! 主散乱体の S(α, β) と有効温度だけを読む。副散乱体は読み飛ばす（NJOY THERMR も、
//! 自由ガス型の副散乱体は使わない）。α の格子は全 β で同じで、α・β それぞれの補間則は
//! 1つであることを求める（LEAPR が作るデータはそうなっている）。

use super::Error;
use super::tab1::{Interpolation, Tab1};
use super::tape::Section;

/// 非干渉性非弾性散乱の S(α, β)。
#[derive(Debug, Clone, PartialEq)]
pub struct Mf7Incoherent {
    /// LAT=1 なら α, β は温度によらず kT = 0.0253 eV で作られている
    pub lat: i64,
    /// LASYM=1 なら S は β について非対称（負の β も持つ）
    pub lasym: i64,
    /// 主散乱体の1原子あたりの自由原子の散乱断面積 (barn)
    pub sigma_free: f64,
    /// 主散乱体の質量比
    pub awr: f64,
    /// S(α, β) を使う入射エネルギーの上限 (eV)
    pub e_max: f64,
    pub alpha: Vec<f64>,
    pub alpha_law: Interpolation,
    pub beta: Vec<f64>,
    pub beta_law: Interpolation,
    pub temperatures: Vec<f64>,
    /// [温度][β の番号 × α の数 + α の番号] の S(α, β)
    pub s: Vec<Vec<f64>>,
    /// 主散乱体の有効温度 (K)（温度の関数）
    pub effective_temperature: Tab1,
}

impl Mf7Incoherent {
    pub fn parse(section: &Section<'_>) -> Result<Self, Error> {
        let _ = section;
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::endf::Tape;
    use crate::endf::tape::tests::rec;

    fn mf7(lln: &str, second_alpha: &str, s: [&str; 3]) -> Result<Mf7Incoherent, Error> {
        let lines: Vec<[&str; 6]> = vec![
            ["1.001000+3", "9.991673-1", "0", "1", "0", "0"],
            // LIST: 0 0 LLN 0 NI NS、B = σfree M0, εmax, A0, Emax, 0, M0
            ["0.0", "0.0", lln, "0", "6", "0"],
            ["4.0", "10.0", "9.9917-1", "5.0", "0.0", "2.0"],
            // TAB2: β は 2 点、ln S が β について線形（INT=4）
            ["0.0", "0.0", "0", "0", "1", "2"],
            ["2", "4", "", "", "", ""],
            // β = 0: TAB1（T0 β LT 0 NR NP）と、2つ目の温度の LIST（T β LI 0 NP 0）
            ["2.936000+2", "0.0", "1", "0", "1", "3"],
            ["3", "4", "", "", "", ""],
            ["0.1", s[0], "1.0", s[1], "10.0", s[2]],
            ["6.000000+2", "0.0", "4", "0", "3", "0"],
            ["0.9", "0.6", "0.02", "", "", ""],
            // β = 2
            ["2.936000+2", "2.0", "1", "0", "1", "3"],
            ["3", "4", "", "", "", ""],
            ["0.1", "0.2", second_alpha, "0.3", "10.0", "0.005"],
            ["6.000000+2", "2.0", "4", "0", "3", "0"],
            ["0.25", "0.35", "0.006", "", "", ""],
            // 有効温度の TAB1
            ["0.0", "0.0", "0", "0", "1", "2"],
            ["2", "2", "", "", "", ""],
            ["2.936000+2", "1400.0", "6.000000+2", "1500.0", "", ""],
        ];
        let text: String = lines.iter().map(|f| rec(*f, 1, 7, 4)).collect();
        let tape = Tape::parse(&text).unwrap();
        Mf7Incoherent::parse(&tape.section(1, 7, 4).unwrap())
    }

    #[test]
    fn reads_principal_scatterer() {
        let data = mf7("0", "1.0", ["1.0", "0.5", "0.01"]).unwrap();
        assert_eq!((data.lat, data.lasym), (1, 0));
        assert_eq!((data.sigma_free, data.awr, data.e_max), (2.0, 0.99917, 5.0));
        assert_eq!(data.alpha, vec![0.1, 1.0, 10.0]);
        assert_eq!(data.beta, vec![0.0, 2.0]);
        assert_eq!(
            (data.alpha_law, data.beta_law),
            (Interpolation::LogLin, Interpolation::LogLin)
        );
        assert_eq!(data.temperatures, vec![293.6, 600.0]);
        assert_eq!(data.s[0], vec![1.0, 0.5, 0.01, 0.2, 0.3, 0.005]);
        assert_eq!(data.s[1], vec![0.9, 0.6, 0.02, 0.25, 0.35, 0.006]);
        assert_eq!(data.effective_temperature.eval(600.0), Some(1500.0));
    }

    #[test]
    fn rejects_logarithm_of_s() {
        // ENDF/B-VIII.0 の熱中性子散乱のデータはすべて LLN=0 なので、ln S での格納は扱わない
        assert!(mf7("1", "1.0", ["1.0", "0.5", "0.01"]).is_err());
    }

    #[test]
    fn alpha_grid_must_be_shared_by_all_beta() {
        assert!(mf7("0", "1.5", ["1.0", "0.5", "0.01"]).is_err());
    }
}
