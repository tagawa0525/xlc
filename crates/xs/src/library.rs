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

impl Library {
    pub fn to_text(&self) -> String {
        todo!()
    }

    pub fn parse(text: &str) -> Result<Self, ParseError> {
        todo!()
    }

    pub fn nuclide(&self, name: &str) -> Option<&Nuclide> {
        todo!()
    }
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
        let text = sample().to_text().replace("total 11.5 ", "total ");
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
