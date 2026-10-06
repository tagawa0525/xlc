//! 炉心の入力ファイル（TOML）。M0 では1種類の少数群定数で均質な直方体の炉心を扱う。
//!
//! ```toml
//! few-group = "pincell-fg.msgpack"   # 少数群定数のファイル
//! size = [300, 300, 366]             # x, y, z の寸法 (cm)
//! nodes = [15, 15, 18]               # x, y, z のノード数（等分割）
//! # 各軸の（低い側, 高い側）の境界条件（reflective または zero-flux）。6面とも書く
//! boundary = { x = ["reflective", "zero-flux"], y = ["zero-flux", "zero-flux"], z = ["zero-flux", "zero-flux"] }
//! ```

use std::path::PathBuf;

use crate::nem::Boundary;

/// 入力ファイルの読み込みエラー。
#[derive(Debug, thiserror::Error)]
pub enum InputError {
    #[error("{0}")]
    Toml(#[from] toml::de::Error),
    #[error("{0}")]
    Invalid(String),
}

#[derive(Debug, Clone, PartialEq)]
pub struct CoreInput {
    pub fewgroup: PathBuf,
    pub size: [f64; 3],
    pub nodes: [usize; 3],
    pub boundaries: [Boundary; 6],
}

impl CoreInput {
    pub fn parse(text: &str) -> Result<Self, InputError> {
        todo!("{text}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = r#"
# bare core
few-group = "pincell-fg.msgpack"
size = [300, 300, 366]
nodes = [15, 15, 18]
boundary = { x = ["reflective", "zero-flux"], y = ["zero-flux", "zero-flux"], z = ["zero-flux", "zero-flux"] }  # x=0 is a symmetry plane
"#;

    fn replaced(from: &str, to: &str) -> Result<CoreInput, InputError> {
        assert!(SAMPLE.contains(from), "{from}");
        CoreInput::parse(&SAMPLE.replace(from, to))
    }

    #[test]
    fn parses_sample() {
        let input = CoreInput::parse(SAMPLE).unwrap();
        assert_eq!(input.fewgroup, PathBuf::from("pincell-fg.msgpack"));
        assert_eq!(input.size, [300.0, 300.0, 366.0]);
        assert_eq!(input.nodes, [15, 15, 18]);
        // -x, +x, -y, +y, -z, +z の順
        assert_eq!(input.boundaries[0], Boundary::Reflective);
        assert!(
            input.boundaries[1..]
                .iter()
                .all(|&b| b == Boundary::ZeroFlux)
        );
    }

    #[test]
    fn invalid_inputs_are_errors() {
        assert!(replaced("size = [300, 300, 366]", "size = [300, 300]").is_err());
        assert!(replaced("size = [300, 300, 366]", "size = [300, -1, 366]").is_err());
        assert!(replaced("nodes = [15, 15, 18]", "nodes = [15, 0, 18]").is_err());
        assert!(replaced("\"reflective\"", "\"vacuum\"").is_err());
        assert!(replaced("z = [\"zero-flux\", \"zero-flux\"] ", "").is_err());
        assert!(replaced("y = [\"zero-flux\", \"zero-flux\"]", "y = [\"zero-flux\"]").is_err());
        // 必須のキーがない、知らないキー、TOML の誤りは位置つき
        assert!(replaced("nodes = [15, 15, 18]\n", "").is_err());
        assert!(replaced("nodes = [15, 15, 18]", "nodes = [15, 15, 18]\npower = 1").is_err());
        let err = replaced("size = [300, 300, 366]", "size = [300, 300").unwrap_err();
        assert!(err.to_string().contains("line"), "{err}");
    }
}
