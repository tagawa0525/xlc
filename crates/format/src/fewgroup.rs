//! 均質化した少数群定数（`lattice` から `core` に渡すファイル）。
//!
//! フィールド名付きの map 形式の MessagePack（ADR-0011、ADR-0013）:
//! `format` ("xlc-fewgroup"), `bounds`（群境界 (eV)、降順）, `k_inf`, `total`, `absorption`, `nu_fission`,
//! `fission`, `chi`, `diffusion`（G 個ずつ）, `scatter`（[入射の群][出射の群]）, `adf`（[辺（左、右、上、下）][群]）。

use serde::{Deserialize, Serialize};

use crate::msgpack::{self, MsgpackError};

/// セル全体で均質化した少数群定数。断面積は 1/cm、拡散係数は cm。
#[derive(Debug, Clone, PartialEq)]
pub struct FewGroup {
    /// 群境界 (eV)、降順
    pub bounds: Vec<f64>,
    pub k_inf: f64,
    pub total: Vec<f64>,
    pub absorption: Vec<f64>,
    pub nu_fission: Vec<f64>,
    pub fission: Vec<f64>,
    pub chi: Vec<f64>,
    pub diffusion: Vec<f64>,
    /// [入射の群][出射の群]
    pub scatter: Vec<Vec<f64>>,
    /// [辺（左、右、上、下）][群] の不連続因子（辺の上の中性子束 / 均質化した中性子束）
    pub adf: Vec<Vec<f64>>,
}

/// 少数群定数のファイルの読み込みエラー。
#[derive(Debug, thiserror::Error)]
pub enum FewGroupError {
    #[error(transparent)]
    Msgpack(#[from] MsgpackError),
    #[error("expected format xlc-fewgroup, found {0}")]
    Format(String),
    #[error("{0}")]
    Invalid(String),
}

/// ファイルの中身（フィールド名付きの map 形式で書く）。
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct FewGroupFile {
    format: String,
    bounds: Vec<f64>,
    k_inf: f64,
    total: Vec<f64>,
    absorption: Vec<f64>,
    nu_fission: Vec<f64>,
    fission: Vec<f64>,
    chi: Vec<f64>,
    diffusion: Vec<f64>,
    scatter: Vec<Vec<f64>>,
    adf: Vec<Vec<f64>>,
}

impl FewGroup {
    pub fn groups(&self) -> usize {
        self.total.len()
    }

    /// MessagePack（フィールド名付きの map 形式）で書く。
    pub fn to_msgpack(&self) -> Result<Vec<u8>, MsgpackError> {
        msgpack::to_named(&FewGroupFile {
            format: FORMAT.into(),
            bounds: self.bounds.clone(),
            k_inf: self.k_inf,
            total: self.total.clone(),
            absorption: self.absorption.clone(),
            nu_fission: self.nu_fission.clone(),
            fission: self.fission.clone(),
            chi: self.chi.clone(),
            diffusion: self.diffusion.clone(),
            scatter: self.scatter.clone(),
            adf: self.adf.clone(),
        })
    }

    /// MessagePack から読み、`validate` で値を確かめる。
    pub fn from_msgpack(bytes: &[u8]) -> Result<Self, FewGroupError> {
        let file: FewGroupFile = msgpack::from_value(msgpack::read_value(bytes)?)?;
        if file.format != FORMAT {
            return Err(FewGroupError::Format(file.format));
        }
        let fewgroup = Self {
            bounds: file.bounds,
            k_inf: file.k_inf,
            total: file.total,
            absorption: file.absorption,
            nu_fission: file.nu_fission,
            fission: file.fission,
            chi: file.chi,
            diffusion: file.diffusion,
            scatter: file.scatter,
            adf: file.adf,
        };
        fewgroup.validate()?;
        Ok(fewgroup)
    }

    /// 群境界の順序、配列の長さ、値が有限であること、不連続因子が正であることを確かめる。
    /// ファイルから読んだ値にも、ファイルを経ずに作った値にも使う。
    pub fn validate(&self) -> Result<(), FewGroupError> {
        let invalid = |message: String| Err(FewGroupError::Invalid(message));
        let b = &self.bounds;
        if b.len() < 2
            || !b.windows(2).all(|w| w[0] > w[1])
            || !b.iter().all(|&x| x > 0.0 && x.is_finite())
        {
            return invalid(format!("bounds must be positive and descending, got {b:?}"));
        }
        let ng = b.len() - 1;
        let vectors = [
            ("total", &self.total),
            ("absorption", &self.absorption),
            ("nu_fission", &self.nu_fission),
            ("fission", &self.fission),
            ("chi", &self.chi),
            ("diffusion", &self.diffusion),
        ];
        for (name, v) in vectors
            .into_iter()
            .chain(self.scatter.iter().map(|row| ("scatter", row)))
        {
            if v.len() != ng {
                return invalid(format!("`{name}` has {} values for {ng} groups", v.len()));
            }
            if let Some(bad) = v.iter().find(|x| !x.is_finite()) {
                return invalid(format!("`{name}` has a non-finite value {bad}"));
            }
        }
        if self.scatter.len() != ng {
            return invalid(format!(
                "`scatter` has {} rows for {ng} groups",
                self.scatter.len()
            ));
        }
        if self.adf.len() != 4 {
            return invalid(format!(
                "`adf` (discontinuity factors) has {} sides, expected 4",
                self.adf.len()
            ));
        }
        for side in &self.adf {
            if side.len() != ng || !side.iter().all(|v| *v > 0.0 && v.is_finite()) {
                return invalid(format!(
                    "`adf` (discontinuity factors) needs {ng} positive values per side"
                ));
            }
        }
        if !(self.k_inf > 0.0 && self.k_inf.is_finite()) {
            return invalid(format!("invalid k-inf {}", self.k_inf));
        }
        Ok(())
    }
}

const FORMAT: &str = "xlc-fewgroup";

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> FewGroup {
        // 10 進で短く書けない値にして、ビット単位の往復を確かめる
        let x = |v: f64| v / 3.0;
        FewGroup {
            bounds: vec![1.0e7, 0.625, 1.0e-5],
            k_inf: x(3.5),
            total: vec![x(1.6), x(2.2)],
            absorption: vec![x(0.03), x(0.3)],
            nu_fission: vec![x(0.02), x(0.4)],
            fission: vec![x(0.008), x(0.16)],
            chi: vec![1.0, 0.0],
            diffusion: vec![x(4.1), x(1.2)],
            scatter: vec![vec![x(1.5), x(0.05)], vec![x(0.002), x(1.9)]],
            adf: vec![vec![x(3.0), x(3.3)]; 4],
            vacuum_ratio: vec![None, Some(vec![x(1.6), x(1.7)]), None, None],
        }
    }

    fn sample_file() -> FewGroupFile {
        let fg = sample();
        FewGroupFile {
            format: "xlc-fewgroup".into(),
            bounds: fg.bounds,
            k_inf: fg.k_inf,
            total: fg.total,
            absorption: fg.absorption,
            nu_fission: fg.nu_fission,
            fission: fg.fission,
            chi: fg.chi,
            diffusion: fg.diffusion,
            scatter: fg.scatter,
            adf: fg.adf,
            vacuum_ratio: fg.vacuum_ratio,
        }
    }

    fn read_modified(change: impl FnOnce(&mut FewGroupFile)) -> Result<FewGroup, FewGroupError> {
        let mut file = sample_file();
        change(&mut file);
        FewGroup::from_msgpack(&rmp_serde::to_vec_named(&file).unwrap())
    }

    #[test]
    fn round_trips_bit_exactly() {
        let fg = sample();
        assert_eq!(
            FewGroup::from_msgpack(&fg.to_msgpack().unwrap()).unwrap(),
            fg
        );
    }

    #[test]
    fn invalid_files_are_rejected() {
        assert!(read_modified(|f| f.format = "xlc-multigroup".into()).is_err());
        assert!(
            read_modified(|f| {
                f.diffusion.pop();
            })
            .is_err()
        );
        assert!(
            read_modified(|f| {
                f.scatter[0].pop();
            })
            .is_err()
        );
        assert!(read_modified(|f| f.bounds.reverse()).is_err());
        assert!(read_modified(|f| f.total[0] = f64::NAN).is_err());
        assert!(read_modified(|f| f.k_inf = f64::INFINITY).is_err());
        assert!(read_modified(|f| f.vacuum_ratio[1] = Some(vec![0.5])).is_err());
        // 配列形式、末尾の余り
        let array = rmp_serde::to_vec(&sample_file()).unwrap();
        assert!(matches!(
            FewGroup::from_msgpack(&array),
            Err(FewGroupError::Msgpack(MsgpackError::NotMap))
        ));
        let mut bytes = sample().to_msgpack().unwrap();
        bytes.push(0xc0);
        assert!(FewGroup::from_msgpack(&bytes).is_err());
    }

    #[test]
    fn validate_checks_values_in_memory() {
        // ファイルを経ずに作った値（core の Problem など）も同じ検査を通す
        assert!(sample().validate().is_ok());
        let mut fg = sample();
        fg.adf.pop();
        assert!(fg.validate().is_err());
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut fg = sample();
            fg.adf[1][0] = bad;
            let err = fg.validate().unwrap_err().to_string();
            assert!(err.contains("discontinuity"), "{err}");
        }
        let mut fg = sample();
        fg.diffusion.pop();
        assert!(fg.validate().is_err());
        let mut fg = sample();
        fg.k_inf = 0.0;
        assert!(fg.validate().is_err());
    }

    #[test]
    fn vacuum_ratios_are_checked() {
        // 真空の辺の J/φ は、出ていく角度中性子束の μ の平均なので 0 < r ≤ 1（ADR-0101）
        let mut fg = sample();
        fg.vacuum_ratio[1] = Some(vec![1.0, 0.5]);
        assert!(fg.validate().is_ok());
        for bad in [0.0, -0.1, 1.0 + 1e-12, f64::NAN, f64::INFINITY] {
            let mut fg = sample();
            fg.vacuum_ratio[1] = Some(vec![bad, 0.5]);
            let err = fg.validate().unwrap_err().to_string();
            assert!(err.contains("vacuum_ratio"), "{err}");
        }
        let mut fg = sample();
        fg.vacuum_ratio.pop();
        assert!(fg.validate().is_err());
        let mut fg = sample();
        fg.vacuum_ratio[0] = Some(vec![0.5]);
        assert!(fg.validate().is_err());
        // 真空の辺がなければ、どの辺も持たない
        let mut fg = sample();
        fg.vacuum_ratio = vec![None; 4];
        assert!(fg.validate().is_ok());
        assert_eq!(
            FewGroup::from_msgpack(&fg.to_msgpack().unwrap()).unwrap(),
            fg
        );
    }
}
