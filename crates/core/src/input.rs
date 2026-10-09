//! 炉心の入力ファイル（TOML）。M0 では1種類の少数群定数で均質な直方体の炉心を扱う。
//!
//! ```toml
//! few-group = "pincell-fg.msgpack"   # 少数群定数のファイル
//! size = [300, 300, 366]             # x, y, z の寸法 (cm)
//! nodes = [15, 15, 18]               # x, y, z のノード数（等分割）
//! # 各軸の（低い側, 高い側）の境界条件（reflective、zero-flux、vacuum、marshak のいずれか）。6面とも書く。
//! # vacuum は少数群定数の真空の辺の比 J/φ を使い、marshak は P1 の関係（J/φ = 1/2）を置く（ADR-0101）
//! boundary = { x = ["reflective", "zero-flux"], y = ["zero-flux", "zero-flux"], z = ["zero-flux", "zero-flux"] }
//! ```

use std::path::PathBuf;

use serde::Deserialize;

use crate::nem::Boundary;

/// 入力ファイルの読み込みエラー。
#[derive(Debug, thiserror::Error)]
pub enum InputError {
    #[error(transparent)]
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

#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct InputFile {
    few_group: PathBuf,
    size: [f64; 3],
    nodes: [usize; 3],
    boundary: BoundaryEntry,
}

/// 各軸の（低い側, 高い側）
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BoundaryEntry {
    x: [Condition; 2],
    y: [Condition; 2],
    z: [Condition; 2],
}

#[derive(Clone, Copy, Deserialize)]
#[serde(rename_all = "kebab-case")]
enum Condition {
    Reflective,
    ZeroFlux,
    Vacuum,
    Marshak,
}

impl CoreInput {
    pub fn parse(text: &str) -> Result<Self, InputError> {
        let file: InputFile = toml::from_str(text)?;
        if !file.size.iter().all(|&s| s > 0.0 && s.is_finite()) {
            return Err(InputError::Invalid(format!(
                "sizes must be positive, got {:?}",
                file.size
            )));
        }
        if file.nodes.contains(&0) {
            return Err(InputError::Invalid(format!(
                "node counts must be positive, got {:?}",
                file.nodes
            )));
        }
        let b = file.boundary;
        let faces = [b.x[0], b.x[1], b.y[0], b.y[1], b.z[0], b.z[1]].map(|c| match c {
            Condition::Reflective => Boundary::Reflective,
            Condition::ZeroFlux => Boundary::ZeroFlux,
            Condition::Vacuum => Boundary::Vacuum,
            Condition::Marshak => Boundary::Marshak,
        });
        Ok(Self {
            fewgroup: file.few_group,
            size: file.size,
            nodes: file.nodes,
            boundaries: faces,
        })
    }
}

#[cfg(test)]
mod tests;
