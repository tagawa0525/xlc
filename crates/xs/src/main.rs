//! ENDF ファイルから多群ライブラリを作るコマンド。
//!
//! 使い方: xs <入力ファイル>
//!
//! 入力ファイルの書式（`#` から行末まではコメント。パスは入力ファイルのディレクトリからの相対パス）:
//!
//! ```text
//! output m1.lib                                   # 出力するライブラリ
//! endf-directory ../endf/neutrons                 # 核種の ENDF ファイルのディレクトリ
//! thermal-directory ../endf/thermal_scatt         # 熱中性子散乱のファイルのディレクトリ
//! sigma0 1e10 1e5 1e4 1e3 300 100 30 10 1         # 背景断面積 (barn)、降順（最初が無限希釈）
//! nuclide H1 n-001_H_001.endf temperatures 293.6 600 thermal tsl-HinH2O.endf
//! nuclide U238 n-092_U_238.endf temperatures 293.6 600 900 1200
//! ```
//!
//! `thermal` を指定した核種は、4 eV 未満の弾性散乱を S(α, β) で、ほかは自由ガス模型で扱う。
//! 群構造は WIMS 69群、重み関数は軽水炉向けの iwt=4 とする。

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use xs::endf::Tape;
use xs::library::{Library, NuclideData};
use xs::multigroup::{GroupStructure, Weight};
use xs::pointwise::DEFAULT_TOLERANCE;
use xs::process::{Settings, process_nuclide, read_thermal_scattering};

/// 入力ファイルを解釈した結果（パスは入力ファイルのディレクトリで解決済み）。
#[derive(Debug, Clone, PartialEq)]
struct Input {
    output: PathBuf,
    sigma0: Vec<f64>,
    nuclides: Vec<NuclideInput>,
}

#[derive(Debug, Clone, PartialEq)]
struct NuclideInput {
    name: String,
    endf: PathBuf,
    temperatures: Vec<f64>,
    thermal: Option<PathBuf>,
}

/// 入力ファイルの `nuclide` の行（パスは未解決）。
struct NuclideLine {
    line: usize,
    name: String,
    file: String,
    temperatures: Vec<f64>,
    thermal: Option<String>,
}

/// 入力ファイルを読む。エラーは「line N: ...」の形の文字列で返す。
fn parse_input(text: &str, base: &Path) -> Result<Input, String> {
    let (mut output, mut endf_directory, mut thermal_directory, mut sigma0) =
        (None, None, None, None);
    let mut nuclides: Vec<NuclideLine> = Vec::new();
    let set_once = |slot: &mut Option<String>, line: usize, key: &str, values: &[&str]| {
        if slot.is_some() {
            return Err(format!("line {line}: `{key}` is given twice"));
        }
        match values {
            [value] => {
                *slot = Some(value.to_string());
                Ok(())
            }
            _ => Err(format!("line {line}: `{key}` takes one value")),
        }
    };
    let positive = |line: usize, token: &str| {
        token
            .parse::<f64>()
            .ok()
            .filter(|v| *v > 0.0 && v.is_finite())
            .ok_or_else(|| format!("line {line}: `{token}` must be a positive number"))
    };
    for (index, raw) in text.lines().enumerate() {
        let line = index + 1;
        let tokens: Vec<&str> = raw
            .split('#')
            .next()
            .unwrap_or("")
            .split_whitespace()
            .collect();
        let Some((&key, values)) = tokens.split_first() else {
            continue;
        };
        match key {
            "output" => set_once(&mut output, line, key, values)?,
            "endf-directory" => set_once(&mut endf_directory, line, key, values)?,
            "thermal-directory" => set_once(&mut thermal_directory, line, key, values)?,
            "sigma0" => {
                if sigma0.is_some() {
                    return Err(format!("line {line}: `sigma0` is given twice"));
                }
                let list = values
                    .iter()
                    .map(|v| positive(line, v))
                    .collect::<Result<Vec<_>, _>>()?;
                if list.is_empty() || list.windows(2).any(|w| w[1] >= w[0]) {
                    return Err(format!("line {line}: `sigma0` must be descending"));
                }
                sigma0 = Some(list);
            }
            "nuclide" => {
                let [name, file, "temperatures", rest @ ..] = values else {
                    return Err(format!(
                        "line {line}: expected `nuclide <name> <file> temperatures <K>...`"
                    ));
                };
                let (temperatures, thermal) = match rest.iter().position(|t| *t == "thermal") {
                    Some(at) => match &rest[at + 1..] {
                        [file] => (&rest[..at], Some(file.to_string())),
                        _ => return Err(format!("line {line}: `thermal` takes one file")),
                    },
                    None => (rest, None),
                };
                let temperatures = temperatures
                    .iter()
                    .map(|t| positive(line, t))
                    .collect::<Result<Vec<_>, _>>()?;
                if temperatures.is_empty() {
                    return Err(format!(
                        "line {line}: {name} needs at least one temperature"
                    ));
                }
                if temperatures
                    .iter()
                    .enumerate()
                    .any(|(i, t)| temperatures[..i].contains(t))
                {
                    return Err(format!("line {line}: {name} has a repeated temperature"));
                }
                if nuclides.iter().any(|n| n.name == *name) {
                    return Err(format!("line {line}: nuclide `{name}` is given twice"));
                }
                nuclides.push(NuclideLine {
                    line,
                    name: name.to_string(),
                    file: file.to_string(),
                    temperatures,
                    thermal,
                });
            }
            _ => return Err(format!("line {line}: unknown keyword `{key}`")),
        }
    }
    let missing = |key: &str| format!("`{key}` is missing");
    let output = output.ok_or_else(|| missing("output"))?;
    let endf_directory = base.join(endf_directory.ok_or_else(|| missing("endf-directory"))?);
    let sigma0 = sigma0.ok_or_else(|| missing("sigma0"))?;
    if nuclides.is_empty() {
        return Err(missing("nuclide"));
    }
    let nuclides = nuclides
        .into_iter()
        .map(|n| {
            let NuclideLine {
                line,
                name,
                file,
                temperatures,
                thermal,
            } = n;
            let thermal = match thermal {
                Some(file) => {
                    let directory = thermal_directory.as_ref().ok_or_else(|| {
                        format!("line {line}: `thermal` needs `thermal-directory`")
                    })?;
                    Some(base.join(directory).join(file))
                }
                None => None,
            };
            Ok(NuclideInput {
                name,
                endf: endf_directory.join(file),
                temperatures,
                thermal,
            })
        })
        .collect::<Result<_, String>>()?;
    Ok(Input {
        output: base.join(output),
        sigma0,
        nuclides,
    })
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.as_slice() {
        [path] => read_input(Path::new(path)).and_then(|input| run(&input)),
        _ => Err("usage: xs <input file>".into()),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("error: {message}");
            ExitCode::FAILURE
        }
    }
}

fn read_input(path: &Path) -> Result<Input, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let base = path.parent().unwrap_or(Path::new("."));
    parse_input(&text, base).map_err(|e| format!("{}: {e}", path.display()))
}

fn run(input: &Input) -> Result<(), String> {
    let settings = Settings {
        groups: GroupStructure::wims69(),
        weight: Weight::lwr(),
        sigma0: input.sigma0.clone(),
        thermal_cutoff: THERMAL_CUTOFF,
        tolerance: DEFAULT_TOLERANCE,
    };
    // 核種ごとの処理は独立なので並列に行う
    let results: Vec<Result<NuclideData, String>> = std::thread::scope(|scope| {
        let handles: Vec<_> = input
            .nuclides
            .iter()
            .map(|nuclide| scope.spawn(|| process_file(nuclide, &settings)))
            .collect();
        handles
            .into_iter()
            .map(|h| {
                h.join()
                    .unwrap_or_else(|_| Err("processing panicked".into()))
            })
            .collect()
    });
    let nuclides = results.into_iter().collect::<Result<Vec<_>, _>>()?;
    let library = Library {
        groups: settings.groups,
        sigma0: settings.sigma0,
        nuclides,
    };
    std::fs::write(&input.output, library.to_text())
        .map_err(|e| format!("{}: {e}", input.output.display()))
}

fn process_file(nuclide: &NuclideInput, settings: &Settings) -> Result<NuclideData, String> {
    let start = Instant::now();
    let read = |path: &Path| -> Result<Tape, String> {
        let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        Tape::parse(&text).map_err(|e| format!("{}: {e}", path.display()))
    };
    let name = &nuclide.name;
    let tape = read(&nuclide.endf)?;
    let thermal = match &nuclide.thermal {
        Some(path) => {
            Some(read_thermal_scattering(&read(path)?).map_err(|e| format!("{name}: {e}"))?)
        }
        None => None,
    };
    let data = process_nuclide(
        &tape,
        name,
        &nuclide.temperatures,
        thermal.as_ref(),
        settings,
    )
    .map_err(|e| format!("{name}: {e}"))?;
    eprintln!("{name}: done in {:.1} s", start.elapsed().as_secs_f64());
    Ok(data)
}

/// この入射エネルギー (eV) 未満の弾性散乱を熱領域の扱いにする（WIMS 69群の熱群の上限）
const THERMAL_CUTOFF: f64 = 4.0;

#[cfg(test)]
mod tests {
    use super::*;

    const INPUT: &str = "\
# M1 のライブラリ
output m1.lib
endf-directory endf
thermal-directory tsl
sigma0 1e10 10 1
nuclide H1 h.endf temperatures 293.6 600 thermal h2o.endf   # 水素
nuclide U238 u.endf temperatures 293.6 900
";

    #[test]
    fn parses_input_relative_to_its_directory() {
        let input = parse_input(INPUT, Path::new("work")).unwrap();
        assert_eq!(input.output, PathBuf::from("work/m1.lib"));
        assert_eq!(input.sigma0, vec![1.0e10, 10.0, 1.0]);
        assert_eq!(
            input.nuclides,
            vec![
                NuclideInput {
                    name: "H1".into(),
                    endf: PathBuf::from("work/endf/h.endf"),
                    temperatures: vec![293.6, 600.0],
                    thermal: Some(PathBuf::from("work/tsl/h2o.endf")),
                },
                NuclideInput {
                    name: "U238".into(),
                    endf: PathBuf::from("work/endf/u.endf"),
                    temperatures: vec![293.6, 900.0],
                    thermal: None,
                },
            ]
        );
    }

    #[test]
    fn rejects_bad_input() {
        let base = Path::new(".");
        let without = |key: &str| {
            INPUT
                .lines()
                .filter(|l| !l.starts_with(key))
                .collect::<Vec<_>>()
                .join("\n")
        };
        // 必須の行がない
        for key in ["output", "endf-directory", "sigma0", "nuclide"] {
            assert!(parse_input(&without(key), base).is_err(), "without {key}");
        }
        // thermal を使うのに thermal-directory がない
        assert!(parse_input(&without("thermal-directory"), base).is_err());
        let replaced = |from: &str, to: &str| parse_input(&INPUT.replace(from, to), base);
        // σ0 が降順でない、正でない
        assert!(replaced("sigma0 1e10 10 1", "sigma0 1 10").is_err());
        assert!(replaced("sigma0 1e10 10 1", "sigma0 1e10 0").is_err());
        // 温度がない、正でない、重複している
        assert!(replaced("temperatures 293.6 900", "temperatures").is_err());
        assert!(replaced("temperatures 293.6 900", "temperatures -1").is_err());
        assert!(replaced("temperatures 293.6 900", "temperatures 900 900").is_err());
        // 核種名の重複、知らないキー
        assert!(replaced("nuclide U238", "nuclide H1").is_err());
        assert!(replaced("output m1.lib", "out m1.lib").is_err());
    }
}
