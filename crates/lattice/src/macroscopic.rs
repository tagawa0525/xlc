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
mod tests {
    use super::*;

    const SAMPLE: &str = r#"
groups = 2

[material.fuel]   # 燃料
total = [0.5, 1.2]
fission = [0.003, 0.1]
nu-fission = [0.008, 0.25]
chi = [1.0, 0.0]
scatter = [[0.45, 0.02], [0.0, 0.9]]

[material.water]
total = [0.6, 2.0]
scatter = [[0.5, 0.098], [0.002, 1.97]]
"#;

    #[test]
    fn reads_materials() {
        let materials = parse(SAMPLE).unwrap();
        let names: Vec<&str> = materials.iter().map(|(n, _)| n.as_str()).collect();
        assert_eq!(names, vec!["fuel", "water"]);
        let fuel = &materials[0].1;
        assert_eq!(fuel.total, vec![0.5, 1.2]);
        assert_eq!(fuel.scatter, vec![vec![0.45, 0.02], vec![0.0, 0.9]]);
        assert_eq!(fuel.fission, vec![0.003, 0.1]);
        assert_eq!(fuel.fission_terms.len(), 1);
        assert_eq!(fuel.fission_terms[0].chi, vec![1.0, 0.0]);
        assert_eq!(fuel.fission_terms[0].nu_fission, vec![0.008, 0.25]);
        let water = &materials[1].1;
        assert_eq!(water.nu_fission, vec![0.0, 0.0]);
        assert!(water.fission_terms.is_empty());
    }

    #[test]
    fn chi_is_normalized() {
        // χ は確率分布。C5G7 の χ は有効数字6桁に丸めてあり、和が 1 + 9e-6 になる
        let text = SAMPLE.replace("chi = [1.0, 0.0]", "chi = [0.58791, 0.41176]");
        let chi = &parse(&text).unwrap()[0].1.fission_terms[0].chi;
        let sum = 0.58791 + 0.41176;
        assert_eq!(chi, &vec![0.58791 / sum, 0.41176 / sum]);
        assert!(parse(&SAMPLE.replace("chi = [1.0, 0.0]", "chi = [0.0, 0.0]")).is_err());
    }

    #[test]
    fn invalid_files_are_errors() {
        let replaced = |from: &str, to: &str| {
            assert!(SAMPLE.contains(from), "{from}");
            parse(&SAMPLE.replace(from, to))
        };
        // 群の数と値の数が合わない、散乱の行が足りない
        let err = replaced("total = [0.6, 2.0]", "total = [0.6]").unwrap_err();
        assert!(err.to_string().contains("water"), "{err}");
        assert!(
            replaced(
                "scatter = [[0.5, 0.098], [0.002, 1.97]]",
                "scatter = [[0.5, 0.098]]"
            )
            .is_err()
        );
        // 核分裂の3つの項目はそろえる（Σf がないとピン出力が黙って 0 になる）
        assert!(replaced("chi = [1.0, 0.0]\n", "").is_err());
        assert!(replaced("fission = [0.003, 0.1]\n", "").is_err());
        assert!(replaced("nu-fission = [0.008, 0.25]\n", "").is_err());
        // 同じ名前の物質（TOML の重複キー）、知らないキー、群の数がない
        assert!(replaced("[material.water]", "[material.fuel]").is_err());
        assert!(
            replaced(
                "total = [0.6, 2.0]",
                "total = [0.6, 2.0]\nabsorption = [0.1, 0.1]"
            )
            .is_err()
        );
        assert!(replaced("groups = 2\n", "").is_err());
        // 全断面積は正、値は有限
        assert!(replaced("total = [0.6, 2.0]", "total = [0.6, 0.0]").is_err());
        assert!(replaced("total = [0.6, 2.0]", "total = [0.6, nan]").is_err());
        assert!(replaced("total = [0.6, 2.0]", "total = [0.6, inf]").is_err());
    }
}
