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
        let mut r = section.reader();
        let head = r.cont()?;
        let (lat, lasym) = (head.l2, head.n1);
        let constants = r.list()?;
        let b = &constants.values;
        let line = r.line();
        if b.len() < 6 || b[5] <= 0.0 {
            return Err(Error::new(
                line,
                "MF7 needs the principal scatterer constants",
            ));
        }
        if constants.l1 != 0 {
            return Err(Error::new(
                line,
                "S(α, β) stored as ln S (LLN=1) is not supported",
            ));
        }
        let (sigma_free, awr, e_max) = (b[0] / b[5], b[2], b[3]);

        let betas = r.tab2()?;
        let beta_law = single_law(&betas.codes, r.line())?;
        let mut alpha: Option<(Vec<f64>, Interpolation)> = None;
        let mut temperatures: Vec<f64> = Vec::new();
        let mut beta = Vec::with_capacity(betas.nz);
        // [β][温度] の S(α) を読み、最後に [温度][β, α] に並べ替える
        let mut by_beta: Vec<Vec<Vec<f64>>> = Vec::with_capacity(betas.nz);
        for _ in 0..betas.nz {
            let first = r.tab1()?;
            let line = r.line();
            let [law] = first.table.laws[..] else {
                return Err(Error::new(
                    line,
                    "MF7 needs a single interpolation law in α",
                ));
            };
            match &alpha {
                None => alpha = Some((first.table.x.clone(), law)),
                Some((grid, grid_law)) if *grid == first.table.x && *grid_law == law => {}
                Some(_) => {
                    return Err(Error::new(line, "α grid differs between β values"));
                }
            }
            beta.push(first.c2);
            let mut rows = vec![first.table.y];
            let mut temps = vec![first.c1];
            for _ in 0..first.l1 {
                let list = r.list()?;
                temps.push(list.c1);
                rows.push(list.values);
            }
            if temperatures.is_empty() {
                temperatures = temps;
            } else if temperatures != temps {
                return Err(Error::new(r.line(), "temperatures differ between β values"));
            }
            if rows.iter().any(|row| row.len() != rows[0].len()) {
                return Err(Error::new(
                    r.line(),
                    "S(α) has a different length at another temperature",
                ));
            }
            by_beta.push(rows);
        }
        let (alpha, alpha_law) =
            alpha.ok_or_else(|| Error::new(r.line(), "MF7 has no β values"))?;
        let s = (0..temperatures.len())
            .map(|t| {
                by_beta
                    .iter()
                    .flat_map(|rows| rows[t].iter().copied())
                    .collect()
            })
            .collect();
        let effective_temperature = r.tab1()?.table;
        Ok(Self {
            lat,
            lasym,
            sigma_free,
            awr,
            e_max,
            alpha,
            alpha_law,
            beta,
            beta_law,
            temperatures,
            s,
            effective_temperature,
        })
    }
}

/// 補間則がただ1つであることを確かめる。
fn single_law(codes: &[i64], line: usize) -> Result<Interpolation, Error> {
    match codes {
        [code] => Interpolation::from_code(*code)
            .ok_or_else(|| Error::new(line, format!("invalid interpolation code {code}"))),
        _ => Err(Error::new(line, "MF7 needs a single interpolation law")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::endf::Tape;
    use crate::endf::tape::tests::rec;

    fn mf7(lln: &str, second_alpha: &str, s: [&str; 3]) -> Result<Mf7Incoherent, Error> {
        mf7_with(lln, second_alpha, s, "0.0", "1")
    }

    /// `list_beta` は β = 0 の2つ目の温度の LIST の β、`lt` は β = 0 の TAB1 の LT。
    fn mf7_with(
        lln: &str,
        second_alpha: &str,
        s: [&str; 3],
        list_beta: &str,
        lt: &str,
    ) -> Result<Mf7Incoherent, Error> {
        let lines: Vec<[&str; 6]> = vec![
            ["1.001000+3", "9.991673-1", "0", "1", "0", "0"],
            // LIST: 0 0 LLN 0 NI NS、B = σfree M0, εmax, A0, Emax, 0, M0
            ["0.0", "0.0", lln, "0", "6", "0"],
            ["4.0", "10.0", "9.9917-1", "5.0", "0.0", "2.0"],
            // TAB2: β は 2 点、ln S が β について線形（INT=4）
            ["0.0", "0.0", "0", "0", "1", "2"],
            ["2", "4", "", "", "", ""],
            // β = 0: TAB1（T0 β LT 0 NR NP）と、2つ目の温度の LIST（T β LI 0 NP 0）
            ["2.936000+2", "0.0", lt, "0", "1", "3"],
            ["3", "4", "", "", "", ""],
            ["0.1", s[0], "1.0", s[1], "10.0", s[2]],
            ["6.000000+2", list_beta, "4", "0", "3", "0"],
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

    #[test]
    fn temperature_records_must_have_the_same_beta() {
        assert!(mf7_with("0", "1.0", ["1.0", "0.5", "0.01"], "0.5", "1").is_err());
    }

    #[test]
    fn number_of_temperatures_must_not_be_negative() {
        assert!(mf7_with("0", "1.0", ["1.0", "0.5", "0.01"], "0.0", "-1").is_err());
    }
}
