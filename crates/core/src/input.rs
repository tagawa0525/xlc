//! 炉心の入力ファイル。M0 では1種類の少数群定数で均質な直方体の炉心を扱う。
//!
//! 1行に `キーワード 値...` を書く。`#` から行末まではコメント。
//!
//! ```text
//! fewgroup pincell.fg          # 少数群定数のファイル
//! size 300 300 366             # x, y, z の寸法 (cm)
//! nodes 15 15 18               # x, y, z のノード数（等分割）
//! boundary zero-flux zero-flux zero-flux zero-flux zero-flux zero-flux
//!                              # -x +x -y +y -z +z（reflective または zero-flux）
//! ```

use std::path::PathBuf;

use crate::nem::Boundary;

const KEYWORDS: [&str; 4] = ["fewgroup", "size", "nodes", "boundary"];

/// 全要素を変換できたときだけ固定長の配列を返す。
fn parse_each<T, const N: usize>(tokens: &[&str], f: impl Fn(&str) -> Option<T>) -> Option<[T; N]> {
    let values: Vec<T> = tokens.iter().map(|t| f(t)).collect::<Option<_>>()?;
    values.try_into().ok()
}

#[derive(Debug, Clone, PartialEq)]
pub struct CoreInput {
    pub fewgroup: PathBuf,
    pub size: [f64; 3],
    pub nodes: [usize; 3],
    pub boundaries: [Boundary; 6],
}

impl CoreInput {
    /// エラーは「line N: ...」の形の文字列で返す。
    pub fn parse(text: &str) -> Result<Self, String> {
        let mut entries: Vec<(&str, usize, Vec<&str>)> = Vec::new();
        for (index, raw) in text.lines().enumerate() {
            let line = index + 1;
            let content = raw.split('#').next().unwrap_or("");
            let tokens: Vec<&str> = content.split_whitespace().collect();
            let Some((&key, values)) = tokens.split_first() else {
                continue;
            };
            if !KEYWORDS.contains(&key) {
                return Err(format!("line {line}: unknown keyword `{key}`"));
            }
            if entries.iter().any(|(k, ..)| *k == key) {
                return Err(format!("line {line}: `{key}` is given twice"));
            }
            entries.push((key, line, values.to_vec()));
        }
        let get = |key: &str, count: usize| -> Result<(usize, &[&str]), String> {
            let (_, line, values) = entries
                .iter()
                .find(|(k, ..)| *k == key)
                .ok_or_else(|| format!("`{key}` is missing"))?;
            if values.len() != count {
                return Err(format!("line {line}: `{key}` takes {count} values"));
            }
            Ok((*line, values.as_slice()))
        };

        let (line, size) = get("size", 3)?;
        let size = parse_each(size, |t| {
            t.parse::<f64>().ok().filter(|v| *v > 0.0 && v.is_finite())
        })
        .ok_or_else(|| format!("line {line}: sizes must be positive numbers"))?;
        let (line, nodes) = get("nodes", 3)?;
        let nodes = parse_each(nodes, |t| t.parse::<usize>().ok().filter(|&n| n > 0))
            .ok_or_else(|| format!("line {line}: node counts must be positive integers"))?;
        let (line, boundaries) = get("boundary", 6)?;
        let boundaries = parse_each(boundaries, |t| match t {
            "reflective" => Some(Boundary::Reflective),
            "zero-flux" => Some(Boundary::ZeroFlux),
            _ => None,
        })
        .ok_or_else(|| format!("line {line}: boundaries must be `reflective` or `zero-flux`"))?;
        Ok(Self {
            fewgroup: PathBuf::from(get("fewgroup", 1)?.1[0]),
            size,
            nodes,
            boundaries,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "
# bare core
fewgroup pincell.fg
size 300 300 366
nodes 15 15 18
boundary reflective zero-flux zero-flux zero-flux zero-flux zero-flux  # x=0 is a symmetry plane
";

    #[test]
    fn parses_sample() {
        let input = CoreInput::parse(SAMPLE).unwrap();
        assert_eq!(input.fewgroup, PathBuf::from("pincell.fg"));
        assert_eq!(input.size, [300.0, 300.0, 366.0]);
        assert_eq!(input.nodes, [15, 15, 18]);
        assert_eq!(input.boundaries[0], Boundary::Reflective);
        assert_eq!(input.boundaries[5], Boundary::ZeroFlux);
    }

    fn line_of(text: &str, key: &str) -> usize {
        text.lines().position(|l| l.starts_with(key)).unwrap() + 1
    }

    #[test]
    fn errors_report_lines() {
        for (from, to, key) in [
            ("size 300 300 366", "size 300 300", "size"),
            ("size 300 300 366", "size 300 -1 366", "size"),
            ("nodes 15 15 18", "nodes 15 0 18", "nodes"),
            ("boundary reflective", "boundary vacuum", "boundary"),
        ] {
            let text = SAMPLE.replace(from, to);
            let err = CoreInput::parse(&text).unwrap_err();
            assert!(
                err.starts_with(&format!("line {}:", line_of(&text, key))),
                "{err}"
            );
        }
    }

    #[test]
    fn missing_unknown_and_duplicate_keys_are_errors() {
        let without: String = SAMPLE
            .lines()
            .filter(|l| !l.starts_with("nodes"))
            .collect::<Vec<_>>()
            .join("\n");
        assert!(CoreInput::parse(&without).unwrap_err().contains("nodes"));
        assert!(CoreInput::parse(&format!("{SAMPLE}\npower 1")).is_err());
        assert!(CoreInput::parse(&format!("{SAMPLE}\nsize 1 1 1")).is_err());
    }
}
