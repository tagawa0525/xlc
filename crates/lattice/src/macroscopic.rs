//! 巨視的断面積のファイル（C5G7 のように物質の多群断面積が直接与えられる問題のため）。
//!
//! 1行に `キーワード 値...` を書く。χ は和が 1 になるよう規格化する。`#` から行末まではコメント。`material` から次の `material` までが
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

use crate::material::{FissionTerm, Material};

/// ファイル中の物質を、名前とともにファイルの順で返す。エラーは「line N: ...」の形の文字列。
pub fn parse(text: &str) -> Result<Vec<(String, Material)>, String> {
    let mut lines = text.lines().enumerate().filter_map(|(index, raw)| {
        let content = raw.split('#').next().unwrap_or("");
        let tokens: Vec<&str> = content.split_whitespace().collect();
        (!tokens.is_empty()).then_some((index + 1, tokens))
    });
    match lines.next() {
        Some((_, tokens)) if tokens == ["xlc-macroscopic", "1"] => {}
        Some((line, _)) => return Err(format!("line {line}: expected `xlc-macroscopic 1`")),
        None => return Err("the file is empty".into()),
    }
    let ng = match lines.next() {
        Some((line, tokens)) => match tokens[..] {
            ["groups", n] => n
                .parse::<usize>()
                .ok()
                .filter(|&n| n > 0)
                .ok_or_else(|| format!("line {line}: invalid number of groups `{n}`"))?,
            _ => return Err(format!("line {line}: expected `groups N`")),
        },
        None => return Err("`groups` is missing".into()),
    };

    let mut materials: Vec<(String, Draft)> = Vec::new();
    for (line, tokens) in lines {
        let (key, values) = (tokens[0], &tokens[1..]);
        if key == "material" {
            let [name] = values else {
                return Err(format!("line {line}: `material` takes one name"));
            };
            if materials.iter().any(|(n, _)| n == name) {
                return Err(format!("line {line}: material `{name}` is given twice"));
            }
            materials.push((name.to_string(), Draft::new(line)));
            continue;
        }
        let (_, draft) = materials
            .last_mut()
            .ok_or_else(|| format!("line {line}: `{key}` comes before any `material`"))?;
        let numbers = values
            .iter()
            .map(|v| v.parse::<f64>().ok().filter(|x| x.is_finite()))
            .collect::<Option<Vec<f64>>>()
            .filter(|v| v.len() == ng)
            .ok_or_else(|| format!("line {line}: `{key}` takes {ng} numbers"))?;
        let slot = match key {
            "total" => &mut draft.total,
            "fission" => &mut draft.fission,
            "nu-fission" => &mut draft.nu_fission,
            "chi" => &mut draft.chi,
            "scatter" => {
                draft.scatter.push(numbers);
                if draft.scatter.len() > ng {
                    return Err(format!("line {line}: more than {ng} `scatter` rows"));
                }
                continue;
            }
            _ => return Err(format!("line {line}: unknown keyword `{key}`")),
        };
        if slot.replace(numbers).is_some() {
            return Err(format!("line {line}: `{key}` is given twice"));
        }
    }
    materials
        .into_iter()
        .map(|(name, draft)| {
            let line = draft.line;
            draft
                .finish(ng)
                .map(|m| (name.clone(), m))
                .map_err(|e| format!("line {line}: material `{name}`: {e}"))
        })
        .collect()
}

/// 読みかけの1つの物質。
struct Draft {
    /// `material` の行
    line: usize,
    total: Option<Vec<f64>>,
    fission: Option<Vec<f64>>,
    nu_fission: Option<Vec<f64>>,
    chi: Option<Vec<f64>>,
    scatter: Vec<Vec<f64>>,
}

impl Draft {
    fn new(line: usize) -> Self {
        Self {
            line,
            total: None,
            fission: None,
            nu_fission: None,
            chi: None,
            scatter: Vec::new(),
        }
    }

    fn finish(self, ng: usize) -> Result<Material, String> {
        let total = self.total.ok_or("`total` is missing")?;
        if total.iter().any(|&t| t <= 0.0) {
            return Err("total cross sections must be positive".into());
        }
        if self.scatter.len() != ng {
            return Err(format!(
                "{} `scatter` rows given, expected {ng}",
                self.scatter.len()
            ));
        }
        let zeros = vec![0.0; ng];
        let nu_fission = self.nu_fission.unwrap_or_else(|| zeros.clone());
        let fission_terms = if nu_fission.iter().any(|&v| v != 0.0) {
            // χ は確率分布なので和で割る（公表値は丸めてあり、和がちょうど 1 にならない）
            let chi = self.chi.ok_or("`nu-fission` needs `chi`")?;
            let sum: f64 = chi.iter().sum();
            if sum <= 0.0 {
                return Err(format!("chi sums to {sum}"));
            }
            let chi = chi.iter().map(|c| c / sum).collect();
            vec![FissionTerm {
                chi,
                nu_fission: nu_fission.clone(),
            }]
        } else {
            Vec::new()
        };
        let material = Material {
            total,
            scatter: self.scatter,
            fission: self.fission.unwrap_or(zeros),
            nu_fission,
            fission_terms,
        };
        material.validate()?;
        Ok(material)
    }
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

    #[test]
    fn chi_is_normalized() {
        // χ は確率分布。C5G7 の χ は有効数字6桁に丸めてあり、和が 1 + 9e-6 になる
        let text = SAMPLE.replace("chi 1.0 0.0", "chi 0.58791 0.41176");
        let chi = &parse(&text).unwrap()[0].1.fission_terms[0].chi;
        let sum = 0.58791 + 0.41176;
        assert_eq!(chi, &vec![0.58791 / sum, 0.41176 / sum]);
        assert!(parse(&SAMPLE.replace("chi 1.0 0.0", "chi 0.0 0.0")).is_err());
    }
}
