//! 少数群定数の分岐のテーブル（`lattice` から `core` に渡すファイル。ADR-0097）。
//!
//! フィールド名付きの map 形式の MessagePack（ADR-0011、ADR-0013）:
//! `format` ("xlc-table"), `bounds`（群境界 (eV)、降順）, `axes`（軸ごとに `name`（`boron`、`moderator`、`fuel`）と
//! `points`（昇順。ppm か K）。`moderator` は格子点ごとの `density` (g/cm³)、`boron` は `b10`（B-10 の原子分率）も持つ）,
//! `split`（集合体を各方向に分けたノードの数。ADR-0108）,
//! `states`（軸の直積の順で、最初の軸がいちばん遅く変わる。状態ごとに split² の部分ノードの定数が続く（部分ノードは
//! -y の行から、各行は -x から）。各状態は少数群定数の `bounds` 以外の項目。`height` は数か nil、`pins` は
//! `dx`、`dy`、`form`、`corners`、`normalized_corners`、`halves` を持つ map か nil で、どちらもキーは省けない）。
//! 軸のないテーブルは状態が1つ。

use serde::{Deserialize, Serialize};

use crate::fewgroup::{FewGroup, FewGroupError, Pins};
use crate::msgpack::{self, MsgpackError};

/// テーブルの軸。点は昇順。
#[derive(Debug, Clone, PartialEq)]
pub enum Axis {
    /// ホウ素濃度 (ppm) と、ホウ素の中の B-10 の原子分率
    Boron { ppm: Vec<f64>, b10: f64 },
    /// 減速材の温度 (K) と、格子点ごとの密度 (g/cm³)
    Moderator {
        temperature: Vec<f64>,
        density: Vec<f64>,
    },
    /// 燃料の温度 (K)
    Fuel { temperature: Vec<f64> },
}

impl Axis {
    pub fn name(&self) -> &'static str {
        match self {
            Self::Boron { .. } => "boron",
            Self::Moderator { .. } => "moderator",
            Self::Fuel { .. } => "fuel",
        }
    }

    pub fn points(&self) -> &[f64] {
        match self {
            Self::Boron { ppm, .. } => ppm,
            Self::Moderator { temperature, .. } | Self::Fuel { temperature } => temperature,
        }
    }
}

/// 軸の並びを確かめる: 同じ軸が2つなく、点が範囲の中で増え、減速材の密度が点ごとに正で、B-10 の割合が 0〜1。
pub fn validate_axes(axes: &[Axis]) -> Result<(), String> {
    let mut names = Vec::new();
    for axis in axes {
        let name = axis.name();
        if names.contains(&name) {
            return Err(format!("axis {name} appears twice"));
        }
        names.push(name);
        let points = axis.points();
        // ホウ素濃度は 0 を含み、温度は正
        let valid = |p: f64| match axis {
            Axis::Boron { .. } => p >= 0.0 && p.is_finite(),
            Axis::Moderator { .. } | Axis::Fuel { .. } => p > 0.0 && p.is_finite(),
        };
        if points.is_empty()
            || !points.iter().all(|&p| valid(p))
            || !points.windows(2).all(|w| w[0] < w[1])
        {
            return Err(format!(
                "axis {name}: points must be valid and increasing, got {points:?}"
            ));
        }
        match axis {
            Axis::Boron { b10, .. } if !(0.0..=1.0).contains(b10) => {
                return Err(format!("axis boron: invalid B-10 atom fraction {b10}"));
            }
            Axis::Moderator {
                temperature,
                density,
            } if density.len() != temperature.len()
                || !density.iter().all(|&d| d > 0.0 && d.is_finite()) =>
            {
                return Err(format!(
                    "axis moderator: {} positive densities are needed, got {density:?}",
                    temperature.len()
                ));
            }
            _ => {}
        }
    }
    Ok(())
}

/// 分岐のテーブル。
#[derive(Debug, Clone, PartialEq)]
pub struct Table {
    pub axes: Vec<Axis>,
    /// 集合体を各方向に分けたノードの数（ADR-0108）。分けなければ 1
    pub split: usize,
    /// 軸の直積の順（最初の軸がいちばん遅く変わる）の状態ごとに、split × split の部分ノードの定数（部分ノードが
    /// いちばん速く変わる。部分ノードは -y の行から、各行は -x から）
    pub states: Vec<FewGroup>,
}

/// テーブルのファイルの読み書きのエラー。
#[derive(Debug, thiserror::Error)]
pub enum TableError {
    #[error(transparent)]
    Msgpack(#[from] MsgpackError),
    #[error("expected format xlc-table, found {0}")]
    Format(String),
    #[error("state {state}")]
    State {
        state: usize,
        #[source]
        source: FewGroupError,
    },
    #[error("{0}")]
    Invalid(String),
}

const FORMAT: &str = "xlc-table";

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct TableFile {
    format: String,
    bounds: Vec<f64>,
    axes: Vec<AxisFile>,
    split: usize,
    states: Vec<StateFile>,
}

/// ファイルの軸（テーブルと媒質で共有する。ADR-0128）
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct AxisFile {
    name: String,
    points: Vec<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    density: Option<Vec<f64>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    b10: Option<f64>,
}

impl From<&Axis> for AxisFile {
    fn from(axis: &Axis) -> Self {
        Self {
            name: axis.name().into(),
            points: axis.points().to_vec(),
            density: match axis {
                Axis::Moderator { density, .. } => Some(density.clone()),
                _ => None,
            },
            b10: match axis {
                Axis::Boron { b10, .. } => Some(*b10),
                _ => None,
            },
        }
    }
}

impl TryFrom<AxisFile> for Axis {
    type Error = String;

    fn try_from(a: AxisFile) -> Result<Self, String> {
        match (a.name.as_str(), a.density, a.b10) {
            ("boron", None, Some(b10)) => Ok(Axis::Boron { ppm: a.points, b10 }),
            ("moderator", Some(density), None) => Ok(Axis::Moderator {
                temperature: a.points,
                density,
            }),
            ("fuel", None, None) => Ok(Axis::Fuel {
                temperature: a.points,
            }),
            (name, ..) => Err(format!(
                "axis {name}: unknown name or wrong fields (boron has b10, moderator has density)"
            )),
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct StateFile {
    k_inf: f64,
    total: Vec<f64>,
    absorption: Vec<f64>,
    nu_fission: Vec<f64>,
    fission: Vec<f64>,
    chi: Vec<f64>,
    diffusion: Vec<f64>,
    scatter: Vec<Vec<f64>>,
    adf: Vec<Vec<f64>>,
    vacuum_ratio: Vec<Option<Vec<f64>>>,
    normalized_adf: Vec<bool>,
    /// nil を書いても、キーは省けない（ADR-0129）
    #[serde(deserialize_with = "nil_or")]
    height: Option<f64>,
    /// nil を書いても、キーは省けない（`deserialize_with` を付けると、serde は省いた `Option` を `None` にしない）
    #[serde(deserialize_with = "nil_or")]
    pins: Option<PinsFile>,
}

/// nil を `None` に読み、キーの欠けはエラーにする（ADR-0129）
fn nil_or<'de, D: serde::Deserializer<'de>, T: Deserialize<'de>>(
    d: D,
) -> Result<Option<T>, D::Error> {
    Option::deserialize(d)
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct PinsFile {
    dx: Vec<f64>,
    dy: Vec<f64>,
    form: Vec<Vec<Vec<f64>>>,
    corners: Vec<Vec<f64>>,
    normalized_corners: Vec<bool>,
    halves: Vec<bool>,
}

impl Table {
    /// 軸のない、状態が1つのテーブル。
    pub fn single(state: FewGroup) -> Self {
        Self {
            axes: Vec::new(),
            split: 1,
            states: vec![state],
        }
    }

    /// 状態ごとの部分ノードの数（split²）
    pub fn sub_nodes(&self) -> usize {
        self.split * self.split
    }

    /// 軸、状態の数、状態の値、状態の間で共通であるべき量（群境界、f* の印、真空の辺、z の厚さ）を確かめる。
    pub fn validate(&self) -> Result<(), TableError> {
        let invalid = |message: String| Err(TableError::Invalid(message));
        validate_axes(&self.axes).map_err(TableError::Invalid)?;
        if self.split == 0 {
            return invalid("split must be positive".into());
        }
        let points: usize = self.axes.iter().map(|a| a.points().len()).product();
        let sub = self.sub_nodes();
        if self.states.len() != points * sub {
            return invalid(format!(
                "{} states for axes of {points} points in all and {sub} sub-nodes per point",
                self.states.len()
            ));
        }
        if self.states.is_empty() {
            return invalid("no states".into());
        }
        for (i, state) in self.states.iter().enumerate() {
            // 状態の間で比べる量は、状態 0 の同じ部分ノードと比べる
            let first = &self.states[i % sub];
            state
                .validate()
                .map_err(|source| TableError::State { state: i, source })?;
            let vacuum = |s: &FewGroup| {
                s.vacuum_ratio
                    .iter()
                    .map(Option::is_some)
                    .collect::<Vec<_>>()
            };
            // 群境界は表全体の量（ファイルには1つだけ書く）なので、状態 0 の最初の定数と比べる
            if state.bounds != self.states[0].bounds
                || state.normalized_adf != first.normalized_adf
                || vacuum(state) != vacuum(first)
                || state.height != first.height
            {
                return invalid(format!(
                    "state {i}: group bounds, f* sides, vacuum sides and the height must equal those of the same sub-node of state 0"
                ));
            }
            // 形状関数は状態ごとに内挿するので、どの状態も持つか持たないかで、セルの境界も同じ
            let cells = |s: &FewGroup| {
                s.pins.as_ref().map(|p| {
                    (
                        p.dx.clone(),
                        p.dy.clone(),
                        p.normalized_corners.clone(),
                        p.halves.clone(),
                    )
                })
            };
            if cells(state) != cells(first) {
                return invalid(format!(
                    "state {i}: pins (form functions) and their cells must equal those of the same sub-node of state 0"
                ));
            }
        }
        Ok(())
    }

    /// MessagePack（フィールド名付きの map 形式）で書く。書く前に `validate` で確かめる。
    pub fn to_msgpack(&self) -> Result<Vec<u8>, TableError> {
        self.validate()?;
        let bounds = self
            .states
            .first()
            .map(|s| s.bounds.clone())
            .unwrap_or_default();
        let axes = self.axes.iter().map(AxisFile::from).collect();
        let states = self
            .states
            .iter()
            .map(|s| StateFile {
                k_inf: s.k_inf,
                total: s.total.clone(),
                absorption: s.absorption.clone(),
                nu_fission: s.nu_fission.clone(),
                fission: s.fission.clone(),
                chi: s.chi.clone(),
                diffusion: s.diffusion.clone(),
                scatter: s.scatter.clone(),
                adf: s.adf.clone(),
                vacuum_ratio: s.vacuum_ratio.clone(),
                normalized_adf: s.normalized_adf.clone(),
                height: s.height,
                pins: s.pins.as_ref().map(|p| PinsFile {
                    dx: p.dx.clone(),
                    dy: p.dy.clone(),
                    form: p.form.clone(),
                    corners: p.corners.clone(),
                    normalized_corners: p.normalized_corners.clone(),
                    halves: p.halves.clone(),
                }),
            })
            .collect();
        Ok(msgpack::to_named(&TableFile {
            format: FORMAT.into(),
            bounds,
            axes,
            split: self.split,
            states,
        })?)
    }

    /// MessagePack から読み、`validate` で値を確かめる。
    pub fn from_msgpack(bytes: &[u8]) -> Result<Self, TableError> {
        let file: TableFile = msgpack::from_value(msgpack::read_value(bytes)?)?;
        if file.format != FORMAT {
            return Err(TableError::Format(file.format));
        }
        let axes = file
            .axes
            .into_iter()
            .map(Axis::try_from)
            .collect::<Result<Vec<_>, _>>()
            .map_err(TableError::Invalid)?;
        let states = file
            .states
            .into_iter()
            .map(|s| FewGroup {
                bounds: file.bounds.clone(),
                k_inf: s.k_inf,
                total: s.total,
                absorption: s.absorption,
                nu_fission: s.nu_fission,
                fission: s.fission,
                chi: s.chi,
                diffusion: s.diffusion,
                scatter: s.scatter,
                adf: s.adf,
                vacuum_ratio: s.vacuum_ratio,
                normalized_adf: s.normalized_adf,
                height: s.height,
                pins: s.pins.map(|p| Pins {
                    dx: p.dx,
                    dy: p.dy,
                    form: p.form,
                    corners: p.corners,
                    normalized_corners: p.normalized_corners,
                    halves: p.halves,
                }),
            })
            .collect();
        let table = Self {
            axes,
            split: file.split,
            states,
        };
        table.validate()?;
        Ok(table)
    }
}

#[cfg(test)]
mod tests;
