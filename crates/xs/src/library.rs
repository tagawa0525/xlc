//! 多群ライブラリ（`xs` から `lattice` に渡すファイル）。
//!
//! 人が読めて diff を取れるテキスト形式。数値は f64 の最短の往復可能な表記で書く。
//! 核種ごとに温度の表を持ち、温度ごとに σ0（背景断面積）ごとの群定数を持つ。σ0 は全核種で共通で、
//! 降順に並べ、最初の値を無限希釈とみなす。
//!
//! ```text
//! xlc-multigroup 2
//! groups <G>
//! bounds <G+1 個の群境界 (eV)、降順>
//! sigma0 <σ0 (barn)、降順>
//! nuclide <名前>
//! za <ZA>
//! awr <AWR>
//! temperature <K>
//! chi <G 個>
//! scatter
//! <G 行 × G 個: 無限希釈の σ(g→g')、行が入射の群>
//! shielded <σ0>          （σ0 の数だけ、sigma0 の順に繰り返す）
//! total <G 個>
//! elastic <G 個>
//! fission <G 個>
//! nu-fission <G 個>
//! temperature <K>        （温度の数だけ繰り返す）
//! ...
//! end
//! ```

use std::fmt;

use crate::multigroup::GroupStructure;

/// 1つの σ0 での群定数 (barn)。
#[derive(Debug, Clone, PartialEq)]
pub struct Shielded {
    pub total: Vec<f64>,
    pub elastic: Vec<f64>,
    pub fission: Vec<f64>,
    pub nu_fission: Vec<f64>,
}

/// 1つの温度での群定数。`shielded` は `Library::sigma0` の順。
#[derive(Debug, Clone, PartialEq)]
pub struct TemperatureData {
    /// 処理温度 (K)
    pub temperature: f64,
    pub chi: Vec<f64>,
    /// 無限希釈の散乱行列 [入射の群][出射の群]
    pub scatter: Vec<Vec<f64>>,
    pub shielded: Vec<Shielded>,
}

/// 1核種の、温度ごとの群定数。
#[derive(Debug, Clone, PartialEq)]
pub struct NuclideData {
    pub name: String,
    pub za: f64,
    pub awr: f64,
    pub temperatures: Vec<TemperatureData>,
}

/// 1核種の、1つの温度での無限希釈の群定数 (barn)。核分裂のない核種では fission, nu-fission, chi は 0。
#[derive(Debug, Clone, PartialEq)]
pub struct Nuclide {
    pub name: String,
    pub za: f64,
    pub awr: f64,
    /// 処理温度 (K)
    pub temperature: f64,
    pub total: Vec<f64>,
    pub fission: Vec<f64>,
    pub nu_fission: Vec<f64>,
    pub chi: Vec<f64>,
    /// [入射の群][出射の群]
    pub scatter: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Library {
    pub groups: GroupStructure,
    /// 背景断面積 (barn)、降順。最初の値を無限希釈とみなす
    pub sigma0: Vec<f64>,
    pub nuclides: Vec<NuclideData>,
}

/// ライブラリの読み込みエラー。`line` は1始まりの行番号。
#[derive(Debug, Clone, PartialEq)]
pub struct ParseError {
    pub line: usize,
    pub message: String,
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "line {}: {}", self.line, self.message)
    }
}

impl std::error::Error for ParseError {}

const FORMAT: &str = "xlc-multigroup";
const VERSION: &str = "2";

impl Library {
    pub fn to_text(&self) -> String {
        todo!()
    }

    pub fn parse(text: &str) -> Result<Self, ParseError> {
        todo!()
    }

    /// 核種 `name` の、温度 `temperature` (K) での無限希釈の群定数。温度は表にあること（補間しない）。
    pub fn nuclide(&self, name: &str, temperature: f64) -> Result<Nuclide, String> {
        todo!()
    }
}

fn numbers_line(keyword: &str, values: &[f64]) -> String {
    let mut out = keyword.to_string();
    for v in values {
        out += &format!(" {v:e}");
    }
    out + "\n"
}

fn error(line: usize, message: impl Into<String>) -> ParseError {
    ParseError {
        line,
        message: message.into(),
    }
}

fn usize_value(line: usize, token: &str) -> Result<usize, ParseError> {
    token
        .parse()
        .map_err(|_| error(line, format!("invalid count `{token}`")))
}

/// 空行を飛ばして、行番号付きでトークンに分けた行を順に返す。
struct Lines<'a> {
    lines: Vec<(usize, Vec<&'a str>)>,
    pos: usize,
}

impl<'a> Lines<'a> {
    fn new(text: &'a str) -> Self {
        let lines = text
            .lines()
            .enumerate()
            .map(|(i, l)| (i + 1, l.split_whitespace().collect::<Vec<_>>()))
            .filter(|(_, tokens)| !tokens.is_empty())
            .collect();
        Self { lines, pos: 0 }
    }

    fn has_more(&self) -> bool {
        self.pos < self.lines.len()
    }

    fn next_line(&mut self, expected: &str) -> Result<(usize, Vec<&'a str>), ParseError> {
        let end = self.lines.last().map_or(1, |(n, _)| n + 1);
        let (line, tokens) =
            self.lines.get(self.pos).cloned().ok_or_else(|| {
                error(end, format!("unexpected end of file, expected {expected}"))
            })?;
        self.pos += 1;
        Ok((line, tokens))
    }

    /// `keyword value` の行（値は1語）
    fn keyword_value(&mut self, keyword: &str) -> Result<(usize, &'a str), ParseError> {
        let (line, tokens) = self.next_line(keyword)?;
        match tokens.as_slice() {
            [k, value] if *k == keyword => Ok((line, value)),
            _ => Err(error(line, format!("expected `{keyword} <value>`"))),
        }
    }

    fn keyword_number(&mut self, keyword: &str) -> Result<f64, ParseError> {
        Ok(self.keyword_numbers(keyword, 1)?.1[0])
    }

    /// `keyword v1 v2 ...` の行（数値はちょうど n 個）
    fn keyword_numbers(
        &mut self,
        keyword: &str,
        n: usize,
    ) -> Result<(usize, Vec<f64>), ParseError> {
        let (line, tokens) = self.next_line(keyword)?;
        if tokens[0] != keyword {
            return Err(error(line, format!("expected `{keyword}`")));
        }
        Ok((line, parse_numbers(line, &tokens[1..], n)?))
    }

    /// `keyword v1 v2 ...` の行（数値の個数は問わない）
    fn keyword_list(&mut self, keyword: &str) -> Result<(usize, Vec<f64>), ParseError> {
        let (line, tokens) = self.next_line(keyword)?;
        if tokens[0] != keyword {
            return Err(error(line, format!("expected `{keyword}`")));
        }
        Ok((line, parse_numbers(line, &tokens[1..], tokens.len() - 1)?))
    }

    /// 数値だけの行（ちょうど n 個）
    fn numbers(&mut self, n: usize) -> Result<Vec<f64>, ParseError> {
        let (line, tokens) = self.next_line("a row of numbers")?;
        parse_numbers(line, &tokens, n)
    }

    /// `temperature` の行に続く、1つの温度の群定数。
    fn temperature_data(
        &mut self,
        temperature: f64,
        ng: usize,
        sigma0: &[f64],
    ) -> Result<TemperatureData, ParseError> {
        let chi = self.keyword_numbers("chi", ng)?.1;
        self.keyword_numbers("scatter", 0)?;
        let scatter = (0..ng)
            .map(|_| self.numbers(ng))
            .collect::<Result<_, _>>()?;
        let mut shielded = Vec::with_capacity(sigma0.len());
        for &s0 in sigma0 {
            let (line, value) = self.keyword_numbers("shielded", 1)?;
            if value[0] != s0 {
                return Err(error(line, format!("expected `shielded {s0:e}`")));
            }
            shielded.push(Shielded {
                total: self.keyword_numbers("total", ng)?.1,
                elastic: self.keyword_numbers("elastic", ng)?.1,
                fission: self.keyword_numbers("fission", ng)?.1,
                nu_fission: self.keyword_numbers("nu-fission", ng)?.1,
            });
        }
        Ok(TemperatureData {
            temperature,
            chi,
            scatter,
            shielded,
        })
    }
}

fn parse_numbers(line: usize, tokens: &[&str], n: usize) -> Result<Vec<f64>, ParseError> {
    if tokens.len() != n {
        return Err(error(
            line,
            format!("expected {n} values, found {}", tokens.len()),
        ));
    }
    tokens
        .iter()
        .map(|t| {
            t.parse::<f64>()
                .ok()
                .filter(|v| v.is_finite())
                .ok_or_else(|| error(line, format!("invalid number `{t}`")))
        })
        .collect()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    fn shielded(scale: f64, fissile: bool) -> Shielded {
        Shielded {
            total: vec![11.5 * scale, 1.0 / 3.0],
            elastic: vec![4.0 * scale, 15.0],
            fission: if fissile {
                vec![1.2 * scale, 580.5]
            } else {
                vec![0.0, 0.0]
            },
            nu_fission: if fissile {
                vec![3.0 * scale, 1414.1]
            } else {
                vec![0.0, 0.0]
            },
        }
    }

    fn temperature(t: f64, fissile: bool) -> TemperatureData {
        TemperatureData {
            temperature: t,
            chi: if fissile {
                vec![1.0, 0.0]
            } else {
                vec![0.0, 0.0]
            },
            scatter: vec![vec![4.0, 0.01], vec![1.0e-7, 15.0]],
            shielded: vec![shielded(1.0, fissile), shielded(0.5, fissile)],
        }
    }

    pub(crate) fn sample() -> Library {
        let nuclide = |name: &str, fissile: bool, temperatures: &[f64]| NuclideData {
            name: name.into(),
            za: 92235.0,
            awr: 233.0248,
            temperatures: temperatures
                .iter()
                .map(|&t| temperature(t, fissile))
                .collect(),
        };
        Library {
            groups: GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]),
            sigma0: vec![1.0e10, 10.0],
            nuclides: vec![
                nuclide("U235", true, &[293.6, 900.0]),
                nuclide("H1", false, &[293.6]),
            ],
        }
    }

    #[test]
    fn round_trips_exactly() {
        let library = sample();
        let text = library.to_text();
        assert!(text.starts_with("xlc-multigroup 2\n"));
        assert_eq!(Library::parse(&text).unwrap(), library);
    }

    #[test]
    fn gives_infinite_dilution_at_a_tabulated_temperature() {
        let library = sample();
        let u235 = library.nuclide("U235", 900.0).unwrap();
        assert_eq!((u235.name.as_str(), u235.temperature), ("U235", 900.0));
        assert_eq!(u235.total, vec![11.5, 1.0 / 3.0]);
        assert_eq!(u235.chi, vec![1.0, 0.0]);
        assert!(library.nuclide("U235", 600.0).is_err());
        assert!(library.nuclide("O16", 293.6).is_err());
    }

    #[test]
    fn wrong_number_of_values_is_error() {
        let text = sample().to_text().replacen("total 1.15e1 ", "total ", 1);
        let err = Library::parse(&text).unwrap_err();
        let line = text.lines().position(|l| l.starts_with("total")).unwrap() + 1;
        assert_eq!(err.line, line);
    }

    #[test]
    fn shielded_blocks_must_follow_sigma0() {
        let text = sample()
            .to_text()
            .replacen("shielded 1e1", "shielded 3e1", 1);
        assert!(Library::parse(&text).is_err());
    }

    #[test]
    fn sigma0_must_be_descending() {
        let text = sample()
            .to_text()
            .replacen("sigma0 1e10 1e1", "sigma0 1e1 1e10", 1);
        assert!(Library::parse(&text).is_err());
    }

    #[test]
    fn unknown_format_is_error() {
        let err = Library::parse("xlc-multigroup 1\n").unwrap_err();
        assert_eq!(err.line, 1);
    }

    #[test]
    fn invalid_group_count_is_error() {
        for groups in ["0", "18446744073709551615"] {
            let text = format!("xlc-multigroup 2\ngroups {groups}\nbounds 1e0\n");
            let err = Library::parse(&text).unwrap_err();
            assert_eq!(err.line, 2, "groups {groups}");
        }
    }

    #[test]
    fn missing_end_is_error() {
        let text = sample().to_text();
        let truncated = text.trim_end().strip_suffix("end").unwrap();
        assert!(Library::parse(truncated).is_err());
    }
}
