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

use crate::material::{FissionTerm, Material};

/// 巨視的断面積のファイルの読み込みエラー。
#[derive(Debug, thiserror::Error)]
pub enum MacroscopicError {
    #[error("{0}")]
    Toml(#[from] toml::de::Error),
    #[error("material `{material}`: {message}")]
    Invalid { material: String, message: String },
}

/// ファイル中の物質を、名前の順で返す。
pub fn parse(text: &str) -> Result<Vec<(String, Material)>, MacroscopicError> {
    todo!("{text}")
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
