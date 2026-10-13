//! 巨視的断面積のファイル（C5G7 のように物質の多群断面積が直接与えられる問題のため）。
//!
//! TOML。物質ごとの表に、断面積 (1/cm) を群の順に書く。`scatter` は入射の群ごとの行。核分裂のない物質は
//! `fission`, `nu-fission`, `chi` を省略できるが、書くなら3つともそろえる。χ は和が 1 になるよう規格化する。
//!
//! ```toml
//! groups = 2
//!
//! [material.fuel]
//! total = [0.5, 1.2]
//! fission = [0.003, 0.1]
//! nu-fission = [0.008, 0.25]
//! chi = [1.0, 0.0]
//! scatter = [[0.45, 0.02], [0.0, 0.9]]   # 群 1 から群 1, 2 へ、群 2 から群 1, 2 へ
//!
//! [material.water]
//! total = [0.6, 2.0]
//! scatter = [[0.5, 0.098], [0.002, 1.97]]
//! ```

use std::collections::BTreeMap;

use serde::Deserialize;

use crate::material::{FissionTerm, Material};

/// 巨視的断面積のファイルの読み込みエラー。
#[derive(Debug, thiserror::Error)]
pub enum MacroscopicError {
    #[error(transparent)]
    Toml(#[from] toml::de::Error),
    #[error("material `{material}`: {message}")]
    Invalid { material: String, message: String },
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct MacroscopicFile {
    groups: usize,
    material: BTreeMap<String, MaterialEntry>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct MaterialEntry {
    total: Vec<f64>,
    /// 核分裂のない物質では省略する（3つともそろえる）
    fission: Option<Vec<f64>>,
    nu_fission: Option<Vec<f64>>,
    chi: Option<Vec<f64>>,
    /// [入射の群][出射の群]
    scatter: Vec<Vec<f64>>,
}

/// ファイル中の物質を、名前の順で返す。
pub fn parse(text: &str) -> Result<Vec<(String, Material)>, MacroscopicError> {
    let file: MacroscopicFile = toml::from_str(text)?;
    file.material
        .into_iter()
        .map(|(name, entry)| {
            entry
                .into_material(file.groups)
                .map(|m| (name.clone(), m))
                .map_err(|message| MacroscopicError::Invalid {
                    material: name,
                    message,
                })
        })
        .collect()
}

impl MaterialEntry {
    fn into_material(self, ng: usize) -> Result<Material, String> {
        let check = |what: &str, v: &[f64]| {
            if v.len() != ng {
                Err(format!("`{what}` has {} values for {ng} groups", v.len()))
            } else if let Some(bad) = v.iter().find(|x| !x.is_finite()) {
                Err(format!("`{what}` has a non-finite value {bad}"))
            } else {
                Ok(())
            }
        };
        check("total", &self.total)?;
        if self.total.iter().any(|&t| t <= 0.0) {
            return Err("total cross sections must be positive".into());
        }
        if self.scatter.len() != ng {
            return Err(format!(
                "`scatter` has {} rows for {ng} groups",
                self.scatter.len()
            ));
        }
        for row in &self.scatter {
            check("scatter", row)?;
        }
        let (fission, nu_fission, fission_terms) = match (self.fission, self.nu_fission, self.chi) {
            (None, None, None) => (vec![0.0; ng], vec![0.0; ng], Vec::new()),
            (Some(fission), Some(nu_fission), Some(chi)) => {
                check("fission", &fission)?;
                check("nu-fission", &nu_fission)?;
                check("chi", &chi)?;
                // χ は確率分布なので和で割る（公表値は丸めてあり、和がちょうど 1 にならない）
                let sum: f64 = chi.iter().sum();
                if sum <= 0.0 {
                    return Err(format!("chi sums to {sum}"));
                }
                let term = FissionTerm {
                    chi: chi.iter().map(|c| c / sum).collect(),
                    nu_fission: nu_fission.clone(),
                };
                (fission, nu_fission, vec![term])
            }
            _ => return Err("give `fission`, `nu-fission` and `chi` together".into()),
        };
        let material = Material {
            total: self.total,
            scatter: self.scatter,
            anisotropic: Vec::new(),
            fission,
            nu_fission,
            fission_terms,
        };
        material.validate().map_err(|e| e.to_string())?;
        Ok(material)
    }
}

#[cfg(test)]
mod tests;
