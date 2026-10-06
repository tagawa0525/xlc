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

use serde::{Deserialize, Serialize};

use crate::multigroup::GroupStructure;

/// 1つの σ0 での群定数 (barn)。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Shielded {
    pub total: Vec<f64>,
    pub elastic: Vec<f64>,
    pub fission: Vec<f64>,
    pub nu_fission: Vec<f64>,
}

/// 1つの温度での群定数。`shielded` は `Library::sigma0` の順。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TemperatureData {
    /// 処理温度 (K)
    pub temperature: f64,
    pub chi: Vec<f64>,
    /// 無限希釈の散乱行列 [入射の群][出射の群]
    pub scatter: Vec<Vec<f64>>,
    pub shielded: Vec<Shielded>,
}

/// 1核種の、温度ごとの群定数。
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
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

/// MessagePack のファイルの形式の名前
const MSGPACK_FORMAT: &str = "xlc-multigroup";

/// MessagePack のライブラリの読み込みエラー。
#[derive(Debug, thiserror::Error)]
pub enum LibraryError {
    #[error("invalid MessagePack: {0}")]
    Decode(#[from] rmp_serde::decode::Error),
    #[error("the file must be a MessagePack map with field names")]
    NotMap,
    #[error("expected format {MSGPACK_FORMAT}, found {0}")]
    Format(String),
    #[error("{0}")]
    Invalid(String),
}

/// ファイルの中身（フィールド名付きの map 形式で書く）。
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct LibraryFile {
    format: String,
    group_bounds: Vec<f64>,
    sigma0: Vec<f64>,
    nuclides: Vec<NuclideData>,
}

impl Library {
    /// MessagePack（フィールド名付きの map 形式）で書く。
    pub fn to_msgpack(&self) -> Vec<u8> {
        todo!()
    }

    /// MessagePack から読み、群境界・σ0 の順序と配列の長さを確かめる。
    pub fn from_msgpack(bytes: &[u8]) -> Result<Self, LibraryError> {
        todo!("{bytes:?}")
    }

    pub fn to_text(&self) -> String {
        let mut out = format!("{FORMAT} {VERSION}\n");
        out += &format!("groups {}\n", self.groups.len());
        out += &numbers_line("bounds", self.groups.bounds());
        out += &numbers_line("sigma0", &self.sigma0);
        for n in &self.nuclides {
            out += &format!("nuclide {}\n", n.name);
            out += &format!("za {:e}\n", n.za);
            out += &format!("awr {:e}\n", n.awr);
            for t in &n.temperatures {
                out += &format!("temperature {:e}\n", t.temperature);
                out += &numbers_line("chi", &t.chi);
                out += "scatter\n";
                for row in &t.scatter {
                    let values: Vec<String> = row.iter().map(|v| format!("{v:e}")).collect();
                    out += &values.join(" ");
                    out += "\n";
                }
                for (s0, s) in self.sigma0.iter().zip(&t.shielded) {
                    out += &format!("shielded {s0:e}\n");
                    out += &numbers_line("total", &s.total);
                    out += &numbers_line("elastic", &s.elastic);
                    out += &numbers_line("fission", &s.fission);
                    out += &numbers_line("nu-fission", &s.nu_fission);
                }
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
        let n_bounds = ng
            .checked_add(1)
            .filter(|_| ng > 0)
            .ok_or_else(|| error(line, format!("invalid number of groups {ng}")))?;
        let (line, bounds) = lines.keyword_numbers("bounds", n_bounds)?;
        if !(bounds.windows(2).all(|w| w[0] > w[1]) && bounds.iter().all(|&b| b > 0.0)) {
            return Err(error(line, "bounds must be positive and descending"));
        }
        let groups = GroupStructure::new(bounds);
        let (line, sigma0) = lines.keyword_list("sigma0")?;
        if sigma0.is_empty()
            || !(sigma0.windows(2).all(|w| w[0] > w[1]) && sigma0.iter().all(|&s| s > 0.0))
        {
            return Err(error(line, "sigma0 must be positive and descending"));
        }
        let mut nuclides = Vec::new();
        while lines.has_more() {
            let (line, name) = lines.keyword_value("nuclide")?;
            if name.is_empty() {
                return Err(error(line, "missing nuclide name"));
            }
            let za = lines.keyword_number("za")?;
            let awr = lines.keyword_number("awr")?;
            let mut temperatures = Vec::new();
            loop {
                let (line, tokens) = lines.next_line("`temperature` or `end`")?;
                match tokens.as_slice() {
                    ["end"] if !temperatures.is_empty() => break,
                    ["temperature", value] => {
                        let temperature = parse_numbers(line, &[value], 1)?[0];
                        temperatures.push(lines.temperature_data(temperature, ng, &sigma0)?);
                    }
                    _ => return Err(error(line, "expected `temperature <K>` or `end`")),
                }
            }
            nuclides.push(NuclideData {
                name: name.to_string(),
                za,
                awr,
                temperatures,
            });
        }
        Ok(Self {
            groups,
            sigma0,
            nuclides,
        })
    }

    /// 核種 `name` の、温度 `temperature` (K) での無限希釈の群定数。温度は表にあること（補間しない）。
    pub fn nuclide(&self, name: &str, temperature: f64) -> Result<Nuclide, String> {
        let data = self
            .nuclides
            .iter()
            .find(|n| n.name == name)
            .ok_or_else(|| format!("nuclide {name} is not in the library"))?;
        let t = data
            .temperatures
            .iter()
            .find(|t| (t.temperature - temperature).abs() <= 1.0e-9 * temperature)
            .ok_or_else(|| {
                let available: Vec<f64> = data.temperatures.iter().map(|t| t.temperature).collect();
                format!("{name} has no data at {temperature} K (available: {available:?})")
            })?;
        let dilute = &t.shielded[0];
        Ok(Nuclide {
            name: data.name.clone(),
            za: data.za,
            awr: data.awr,
            temperature: t.temperature,
            total: dilute.total.clone(),
            fission: dilute.fission.clone(),
            nu_fission: dilute.nu_fission.clone(),
            chi: t.chi.clone(),
            scatter: t.scatter.clone(),
        })
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
            // 1.7929089665572825e-9 は、10 進表記を経ると最下位のビットがずれやすい値
            total: vec![11.5 * scale, 1.7929089665572825e-9],
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
                nuclide("U-235", true, &[293.6, 900.0]),
                nuclide("H-1", false, &[293.6]),
            ],
        }
    }

    #[test]
    fn round_trips_bit_exactly() {
        let library = sample();
        assert_eq!(
            Library::from_msgpack(&library.to_msgpack()).unwrap(),
            library
        );
    }

    #[test]
    fn gives_infinite_dilution_at_a_tabulated_temperature() {
        let library = sample();
        let u235 = library.nuclide("U-235", 900.0).unwrap();
        assert_eq!((u235.name.as_str(), u235.temperature), ("U-235", 900.0));
        assert_eq!(u235.total, vec![11.5, 1.7929089665572825e-9]);
        assert_eq!(u235.chi, vec![1.0, 0.0]);
        assert!(library.nuclide("U-235", 600.0).is_err());
        assert!(library.nuclide("O-16", 293.6).is_err());
    }

    fn sample_file() -> LibraryFile {
        let library = sample();
        LibraryFile {
            format: "xlc-multigroup".into(),
            group_bounds: library.groups.bounds().to_vec(),
            sigma0: library.sigma0,
            nuclides: library.nuclides,
        }
    }

    /// サンプルを書き換え、フィールド名付きの map 形式で書いて読む。
    fn read_modified(change: impl FnOnce(&mut LibraryFile)) -> Result<Library, LibraryError> {
        let mut file = sample_file();
        change(&mut file);
        Library::from_msgpack(&rmp_serde::to_vec_named(&file).unwrap())
    }

    #[test]
    fn written_as_a_map_with_the_format_name() {
        let bytes = sample().to_msgpack();
        let value: std::collections::BTreeMap<String, serde::de::IgnoredAny> =
            rmp_serde::from_slice(&bytes).unwrap();
        let keys: Vec<&str> = value.keys().map(String::as_str).collect();
        assert_eq!(keys, vec!["format", "group_bounds", "nuclides", "sigma0"]);
        #[derive(Deserialize)]
        struct Head {
            format: String,
        }
        let head: Head = rmp_serde::from_slice(&bytes).unwrap();
        assert_eq!(head.format, "xlc-multigroup");
    }

    #[test]
    fn array_form_is_rejected() {
        // 既定の配列形式ではフィールド名がなく、順番の変わったファイルを黙って誤読しうる
        let bytes = rmp_serde::to_vec(&sample_file()).unwrap();
        assert!(matches!(
            Library::from_msgpack(&bytes),
            Err(LibraryError::NotMap)
        ));
    }

    #[test]
    fn unknown_and_missing_fields_are_rejected() {
        #[derive(Serialize)]
        struct Extra {
            format: &'static str,
            group_bounds: Vec<f64>,
            sigma0: Vec<f64>,
            nuclides: Vec<NuclideData>,
            comment: &'static str,
        }
        #[derive(Serialize)]
        struct Missing {
            format: &'static str,
            group_bounds: Vec<f64>,
            nuclides: Vec<NuclideData>,
        }
        let file = sample_file();
        let extra = Extra {
            format: "xlc-multigroup",
            group_bounds: file.group_bounds.clone(),
            sigma0: file.sigma0.clone(),
            nuclides: file.nuclides.clone(),
            comment: "x",
        };
        let missing = Missing {
            format: "xlc-multigroup",
            group_bounds: file.group_bounds.clone(),
            nuclides: file.nuclides.clone(),
        };
        for bytes in [
            rmp_serde::to_vec_named(&extra).unwrap(),
            rmp_serde::to_vec_named(&missing).unwrap(),
        ] {
            assert!(matches!(
                Library::from_msgpack(&bytes),
                Err(LibraryError::Decode(_))
            ));
        }
    }

    #[test]
    fn other_formats_are_rejected() {
        let err = read_modified(|f| f.format = "xlc-fewgroup".into()).unwrap_err();
        assert!(matches!(err, LibraryError::Format(_)), "{err}");
    }

    #[test]
    fn wrong_number_of_values_is_error() {
        let err = read_modified(|f| {
            f.nuclides[0].temperatures[0].shielded[0].total.pop();
        })
        .unwrap_err();
        assert!(err.to_string().contains("U-235"), "{err}");
        assert!(
            read_modified(|f| f.nuclides[1].temperatures[0].scatter[1]
                .pop()
                .map(|_| ())
                .unwrap())
            .is_err()
        );
        assert!(read_modified(|f| f.nuclides[1].temperatures[0].chi.clear()).is_err());
    }

    #[test]
    fn shielded_sets_must_match_sigma0() {
        assert!(
            read_modified(|f| {
                f.nuclides[0].temperatures[0].shielded.pop();
            })
            .is_err()
        );
    }

    #[test]
    fn bounds_and_sigma0_must_be_descending() {
        assert!(read_modified(|f| f.sigma0 = vec![10.0, 1.0e10]).is_err());
        assert!(read_modified(|f| f.group_bounds = vec![1.0e-5, 0.625, 1.0e7]).is_err());
        assert!(read_modified(|f| f.group_bounds = vec![1.0e7]).is_err());
    }

    #[test]
    fn nuclides_need_temperatures() {
        assert!(read_modified(|f| f.nuclides[1].temperatures.clear()).is_err());
    }

    #[test]
    fn broken_bytes_are_rejected() {
        assert!(matches!(
            Library::from_msgpack(&[0x81, 0xc1]),
            Err(LibraryError::Decode(_))
        ));
        assert!(Library::from_msgpack(&[]).is_err());
    }
}
