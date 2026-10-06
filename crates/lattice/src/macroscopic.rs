//! 巨視的断面積のファイル（C5G7 のように物質の多群断面積が直接与えられる問題のため）。
//!
//! 1行に `キーワード 値...` を書く。`#` から行末まではコメント。`material` から次の `material` までが
//! 1つの物質で、`scatter` は入射の群の順に群の数だけ書く。核分裂のない物質は `fission`, `nu-fission`,
//! `chi` を省略できる。
//!
//! ```text
//! xlc-macroscopic 1
//! groups 2
//! material fuel
//! total 0.5 1.2                 # 1/cm
//! fission 0.003 0.1
//! nu-fission 0.008 0.25
//! chi 1.0 0.0
//! scatter 0.45 0.02             # 群 1 から群 1, 2 へ
//! scatter 0.0 0.9
//! ```

use crate::material::Material;

/// ファイル中の物質を、名前とともにファイルの順で返す。エラーは「line N: ...」の形の文字列。
pub fn parse(text: &str) -> Result<Vec<(String, Material)>, String> {
    todo!("{text}")
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "
xlc-macroscopic 1
groups 2
# 燃料
material fuel
total 0.5 1.2
fission 0.003 0.1
nu-fission 0.008 0.25
chi 1.0 0.0
scatter 0.45 0.02
scatter 0.0 0.9
material water
total 0.6 2.0
scatter 0.5 0.098
scatter 0.002 1.97
";

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

    fn line_of(text: &str, prefix: &str) -> usize {
        text.lines().position(|l| l.starts_with(prefix)).unwrap() + 1
    }

    #[test]
    fn errors_report_lines() {
        // 群の数と値の数が合わない
        let text = SAMPLE.replace("total 0.6 2.0", "total 0.6");
        let err = parse(&text).unwrap_err();
        assert!(
            err.starts_with(&format!("line {}:", line_of(&text, "total 0.6"))),
            "{err}"
        );
        // 散乱の行が足りない、χ がないのに νΣf がある、同じ名前の物質、知らないキーワード
        assert!(parse(&SAMPLE.replace("scatter 0.002 1.97\n", "")).is_err());
        assert!(parse(&SAMPLE.replace("chi 1.0 0.0\n", "")).is_err());
        assert!(parse(&SAMPLE.replace("material water", "material fuel")).is_err());
        assert!(parse(&format!("{SAMPLE}\nalbedo 1")).is_err());
        // 形式の版と群の数がない
        assert!(parse(&SAMPLE.replace("xlc-macroscopic 1", "")).is_err());
        assert!(parse(&SAMPLE.replace("groups 2", "")).is_err());
        // 全断面積は正
        assert!(parse(&SAMPLE.replace("total 0.6 2.0", "total 0.6 0.0")).is_err());
    }
}
