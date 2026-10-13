//! 細かい群のまま均質化した媒質（単一集合体の格子計算から、軸方向の反射体の平板の格子計算に渡すファイル。ADR-0128）。
//!
//! フィールド名付きの map 形式の MessagePack（ADR-0011、ADR-0013）:
//! `format` ("xlc-medium"), `bounds`（ライブラリの群境界 (eV)、降順）, `axes`（分岐のテーブルと同じ軸。ADR-0097）,
//! `states`（軸の直積の順で、最初の軸がいちばん遅く変わる。状態ごとに `total`、`scatter`（[モーメント（P0 から）]
//! [入射の群][出射の群]）、`fission`、`nu_fission`、`chi`（核分裂がなければ 0））。軸のない媒質は状態が1つ。

use serde::{Deserialize, Serialize};

use crate::group::GroupStructure;
use crate::msgpack::{self, MsgpackError};
use crate::table::{Axis, AxisFile, validate_axes};

/// χ の和と 1 の差の許容値（生成率で重みづけして和で割った値なので、丸めの程度）
const CHI_TOLERANCE: f64 = 1.0e-12;

/// 分岐の状態ごとの、細かい群の均質な媒質。
#[derive(Debug, Clone, PartialEq)]
pub struct Medium {
    /// 群境界 (eV)、降順
    pub bounds: Vec<f64>,
    pub axes: Vec<Axis>,
    /// 軸の直積の順（最初の軸がいちばん遅く変わる）
    pub states: Vec<MediumState>,
}

/// 1つの状態の媒質の巨視的断面積 (1/cm)。
#[derive(Debug, Clone, PartialEq)]
pub struct MediumState {
    pub total: Vec<f64>,
    /// [モーメント（P0 から）][入射の群][出射の群] の散乱行列
    pub scatter: Vec<Vec<Vec<f64>>>,
    pub fission: Vec<f64>,
    pub nu_fission: Vec<f64>,
    /// 核分裂スペクトル（核種ごとの χ を生成率で重みづけした1つ。和が 1、核分裂がなければ 0）
    pub chi: Vec<f64>,
}

/// 媒質のファイルの読み書きのエラー。
#[derive(Debug, thiserror::Error)]
pub enum MediumError {
    #[error(transparent)]
    Msgpack(#[from] MsgpackError),
    #[error("expected format xlc-medium, found {0}")]
    Format(String),
    #[error("{0}")]
    Invalid(String),
}

const FORMAT: &str = "xlc-medium";

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct MediumFile {
    format: String,
    bounds: Vec<f64>,
    axes: Vec<AxisFile>,
    states: Vec<StateFile>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct StateFile {
    total: Vec<f64>,
    scatter: Vec<Vec<Vec<f64>>>,
    fission: Vec<f64>,
    nu_fission: Vec<f64>,
    chi: Vec<f64>,
}

impl MediumState {
    pub fn groups(&self) -> usize {
        self.total.len()
    }

    /// 配列の長さ、値が有限であること、全断面積が正で、核分裂の量が負でないこと、χ の規格化を確かめる。
    fn validate(&self, ng: usize) -> Result<(), String> {
        let vectors = [
            ("total", &self.total),
            ("fission", &self.fission),
            ("nu_fission", &self.nu_fission),
            ("chi", &self.chi),
        ];
        for (name, v) in vectors {
            if v.len() != ng {
                return Err(format!("`{name}` has {} values for {ng} groups", v.len()));
            }
            if let Some(bad) = v.iter().find(|x| !x.is_finite()) {
                return Err(format!("`{name}` has a non-finite value {bad}"));
            }
        }
        if self.total.iter().any(|&t| t <= 0.0) {
            return Err(format!("`total` must be positive, got {:?}", self.total));
        }
        for (name, v) in [("fission", &self.fission), ("nu_fission", &self.nu_fission)] {
            if v.iter().any(|&x| x < 0.0) {
                return Err(format!("`{name}` must not be negative, got {v:?}"));
            }
        }
        if self.scatter.is_empty() {
            return Err("`scatter` needs at least the P0 moment".into());
        }
        if self
            .scatter
            .iter()
            .any(|m| m.len() != ng || m.iter().any(|row| row.len() != ng))
        {
            return Err(format!("`scatter` moments must be {ng} x {ng}"));
        }
        if let Some(bad) = self
            .scatter
            .iter()
            .flatten()
            .flatten()
            .find(|x| !x.is_finite())
        {
            return Err(format!("`scatter` has a non-finite value {bad}"));
        }
        let sum: f64 = self.chi.iter().sum();
        if self.nu_fission.iter().any(|&x| x > 0.0) {
            if self.chi.iter().any(|&c| c < 0.0) || (sum - 1.0).abs() > CHI_TOLERANCE {
                return Err(format!(
                    "`chi` must be non-negative and sum to 1 with fission, got {:?}",
                    self.chi
                ));
            }
        } else if self.chi.iter().any(|&c| c != 0.0) {
            return Err(format!(
                "`chi` must be 0 without fission, got {:?}",
                self.chi
            ));
        }
        Ok(())
    }
}

impl Medium {
    pub fn groups(&self) -> usize {
        self.bounds.len().saturating_sub(1)
    }

    /// 群境界、軸、状態の数、状態の値、状態の間で散乱のモーメントの数が同じことを確かめる。
    pub fn validate(&self) -> Result<(), MediumError> {
        let invalid = |message: String| Err(MediumError::Invalid(message));
        GroupStructure::new(self.bounds.clone())
            .map_err(|e| MediumError::Invalid(e.to_string()))?;
        validate_axes(&self.axes).map_err(MediumError::Invalid)?;
        // 制御棒の値の定数は、平板の連続の軸の各値で同じ媒質の状態を使って作る（ADR-0128）
        if self.axes.iter().any(|a| matches!(a, Axis::Rod { .. })) {
            return invalid("media have no rod axis (ADR-0128)".into());
        }
        let points: usize = self.axes.iter().map(Axis::len).product();
        if self.states.len() != points {
            return invalid(format!(
                "{} states for axes of {points} points in all",
                self.states.len()
            ));
        }
        let ng = self.groups();
        for (i, state) in self.states.iter().enumerate() {
            state
                .validate(ng)
                .map_err(|e| MediumError::Invalid(format!("state {i}: {e}")))?;
            if state.scatter.len() != self.states[0].scatter.len() {
                return invalid(format!(
                    "state {i} has {} scattering moments, state 0 has {}",
                    state.scatter.len(),
                    self.states[0].scatter.len()
                ));
            }
        }
        Ok(())
    }

    /// MessagePack（フィールド名付きの map 形式）で書く。書く前に `validate` で確かめる。
    pub fn to_msgpack(&self) -> Result<Vec<u8>, MediumError> {
        self.validate()?;
        let states = self
            .states
            .iter()
            .map(|s| StateFile {
                total: s.total.clone(),
                scatter: s.scatter.clone(),
                fission: s.fission.clone(),
                nu_fission: s.nu_fission.clone(),
                chi: s.chi.clone(),
            })
            .collect();
        Ok(msgpack::to_named(&MediumFile {
            format: FORMAT.into(),
            bounds: self.bounds.clone(),
            axes: self.axes.iter().map(AxisFile::from).collect(),
            states,
        })?)
    }

    /// MessagePack から読み、`validate` で値を確かめる。
    pub fn from_msgpack(bytes: &[u8]) -> Result<Self, MediumError> {
        let file: MediumFile = msgpack::from_value(msgpack::read_value(bytes)?)?;
        if file.format != FORMAT {
            return Err(MediumError::Format(file.format));
        }
        let axes = file
            .axes
            .into_iter()
            .map(Axis::try_from)
            .collect::<Result<Vec<_>, _>>()
            .map_err(MediumError::Invalid)?;
        let medium = Self {
            bounds: file.bounds,
            axes,
            states: file
                .states
                .into_iter()
                .map(|s| MediumState {
                    total: s.total,
                    scatter: s.scatter,
                    fission: s.fission,
                    nu_fission: s.nu_fission,
                    chi: s.chi,
                })
                .collect(),
        };
        medium.validate()?;
        Ok(medium)
    }
}

#[cfg(test)]
mod tests;
