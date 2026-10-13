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
            let lt = usize::try_from(first.l1)
                .map_err(|_| Error::new(line, format!("invalid LT {}", first.l1)))?;
            for _ in 0..lt {
                let list = r.list()?;
                if list.c2 != first.c2 {
                    return Err(Error::new(
                        r.line(),
                        format!(
                            "β {} differs from {} of its first temperature",
                            list.c2, first.c2
                        ),
                    ));
                }
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
mod tests;
