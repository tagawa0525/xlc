//! 多群ライブラリ（`xs` から `lattice` に渡すファイル）。
//!
//! 人が読めて diff を取れるテキスト形式。数値は f64 の最短の往復可能な表記で書く。
//!
//! ```text
//! xlc-multigroup 1
//! groups <G>
//! bounds <G+1 個の群境界 (eV)、降順>
//! nuclide <名前>
//! za <ZA>
//! awr <AWR>
//! temperature <K>
//! total <G 個>
//! fission <G 個>
//! nu-fission <G 個>
//! chi <G 個>
//! scatter
//! <G 行 × G 個: σ(g→g')、行が入射の群>
//! end
//! ```

use std::fmt;

use crate::multigroup::GroupStructure;

/// 1核種の多群定数 (barn)。核分裂のない核種では fission, nu-fission, chi は 0。
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
    pub nuclides: Vec<Nuclide>,
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
const VERSION: &str = "1";

impl Library {
    pub fn to_text(&self) -> String {
        let mut out = format!("{FORMAT} {VERSION}\n");
        out += &format!("groups {}\n", self.groups.len());
        out += &numbers_line("bounds", self.groups.bounds());
        for n in &self.nuclides {
            out += &format!("nuclide {}\n", n.name);
            out += &format!("za {:e}\n", n.za);
            out += &format!("awr {:e}\n", n.awr);
            out += &format!("temperature {:e}\n", n.temperature);
            out += &numbers_line("total", &n.total);
            out += &numbers_line("fission", &n.fission);
            out += &numbers_line("nu-fission", &n.nu_fission);
            out += &numbers_line("chi", &n.chi);
            out += "scatter\n";
            for row in &n.scatter {
                let values: Vec<String> = row.iter().map(|v| format!("{v:e}")).collect();
                out += &values.join(" ");
                out += "\n";
            }
            out += "end\n";
        }
        out
    }

    pub fn parse(text: &str) -> Result<Self, ParseError> {
        let mut lines = Lines::new(text);
        let (line, header) = lines.next_line("format header")?;
        if header != [FORMAT, VERSION] {
            return Err(error(line, format!("expected `{FORMAT} {VERSION}`")));
        }
        let (line, ng) = lines.keyword_value("groups")?;
        let ng = usize_value(line, ng)?;
        let (line, bounds) = lines.keyword_numbers("bounds", ng + 1)?;
        if !(bounds.windows(2).all(|w| w[0] > w[1]) && bounds.iter().all(|&b| b > 0.0)) {
            return Err(error(line, "bounds must be positive and descending"));
        }
        let groups = GroupStructure::new(bounds);
        let mut nuclides = Vec::new();
        while lines.has_more() {
            let (line, name) = lines.keyword_value("nuclide")?;
            if name.is_empty() {
                return Err(error(line, "missing nuclide name"));
            }
            let za = lines.keyword_number("za")?;
            let awr = lines.keyword_number("awr")?;
            let temperature = lines.keyword_number("temperature")?;
            let total = lines.keyword_numbers("total", ng)?.1;
            let fission = lines.keyword_numbers("fission", ng)?.1;
            let nu_fission = lines.keyword_numbers("nu-fission", ng)?.1;
            let chi = lines.keyword_numbers("chi", ng)?.1;
            lines.keyword_numbers("scatter", 0)?;
            let scatter = (0..ng)
                .map(|_| lines.numbers(ng))
                .collect::<Result<_, _>>()?;
            lines.keyword_numbers("end", 0)?;
            nuclides.push(Nuclide {
                name: name.to_string(),
                za,
                awr,
                temperature,
                total,
                fission,
                nu_fission,
                chi,
                scatter,
            });
        }
        Ok(Self { groups, nuclides })
    }

    pub fn nuclide(&self, name: &str) -> Option<&Nuclide> {
        self.nuclides.iter().find(|n| n.name == name)
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

    /// 数値だけの行（ちょうど n 個）
    fn numbers(&mut self, n: usize) -> Result<Vec<f64>, ParseError> {
        let (line, tokens) = self.next_line("a row of numbers")?;
        parse_numbers(line, &tokens, n)
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
mod tests {
    use super::*;

    fn sample() -> Library {
        let groups = GroupStructure::new(vec![1.0e7, 0.625, 1.0e-5]);
        let nuclide = |name: &str, fissile: bool| Nuclide {
            name: name.into(),
            za: 92235.0,
            awr: 233.0248,
            temperature: 293.6,
            total: vec![11.5, 1.0 / 3.0],
            fission: if fissile {
                vec![1.2, 580.5]
            } else {
                vec![0.0, 0.0]
            },
            nu_fission: if fissile {
                vec![3.0, 1414.1]
            } else {
                vec![0.0, 0.0]
            },
            chi: if fissile {
                vec![1.0, 0.0]
            } else {
                vec![0.0, 0.0]
            },
            scatter: vec![vec![4.0, 0.01], vec![1.0e-7, 15.0]],
        };
        Library {
            groups,
            nuclides: vec![nuclide("U235", true), nuclide("H1", false)],
        }
    }

    #[test]
    fn round_trips_exactly() {
        let library = sample();
        let text = library.to_text();
        assert!(text.starts_with("xlc-multigroup 1\n"));
        assert_eq!(Library::parse(&text).unwrap(), library);
    }

    #[test]
    fn finds_nuclide_by_name() {
        let library = sample();
        assert_eq!(library.nuclide("H1").unwrap().name, "H1");
        assert!(library.nuclide("O16").is_none());
    }

    #[test]
    fn wrong_number_of_values_is_error() {
        let text = sample().to_text().replace("total 1.15e1 ", "total ");
        let err = Library::parse(&text).unwrap_err();
        let line = text.lines().position(|l| l.starts_with("total")).unwrap() + 1;
        assert_eq!(err.line, line);
    }

    #[test]
    fn unknown_format_is_error() {
        let err = Library::parse("xlc-multigroup 2\n").unwrap_err();
        assert_eq!(err.line, 1);
    }

    #[test]
    fn missing_end_is_error() {
        let text = sample().to_text();
        let truncated = text.trim_end().strip_suffix("end").unwrap();
        assert!(Library::parse(truncated).is_err());
    }
}
