//! ENDF ファイルから多群ライブラリを作るコマンド。
//!
//! 使い方: xs <入力ファイル>
//!
//! 入力ファイルは TOML（パスは入力ファイルのディレクトリからの相対パス）:
//!
//! ```toml
//! output = "m1.msgpack"                       # 出力するライブラリ
//! endf-directory = "../endf/neutrons"         # 核種の ENDF ファイルのディレクトリ
//! thermal-directory = "../endf/thermal_scatt" # 熱中性子散乱のファイルのディレクトリ（thermal を使うとき）
//! sigma0 = [1e10, 1e5, 1e4, 1e3, 300, 100, 30, 10, 1]   # 背景断面積 (barn)、降順（最初が無限希釈）
//!
//! [nuclide.H-1]
//! endf = "n-001_H_001.endf"
//! temperatures = [293.6, 600]
//! thermal = "tsl-HinH2O.endf"
//!
//! [nuclide.U-238]
//! endf = "n-092_U_238.endf"
//! temperatures = [293.6, 600, 900, 1200]
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

/// 入力ファイルの読み込みエラー。
#[derive(Debug, thiserror::Error)]
enum InputError {
    #[error("{0}")]
    Toml(#[from] toml::de::Error),
    #[error("{0}")]
    Invalid(String),
}

/// 入力ファイルを読む（パスは `base` で解決する）。
fn parse_input(text: &str, base: &Path) -> Result<Input, InputError> {
    todo!("{text} {base:?}")
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
    std::fs::write(&input.output, library.to_msgpack())
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

    const INPUT: &str = r#"
# M1 のライブラリ
output = "m1.msgpack"
endf-directory = "endf"
thermal-directory = "tsl"
sigma0 = [1e10, 10, 1]

[nuclide.H-1]          # 水素
endf = "h.endf"
temperatures = [293.6, 600]
thermal = "h2o.endf"

[nuclide.U-238]
endf = "u.endf"
temperatures = [293.6, 900]
"#;

    #[test]
    fn parses_input_relative_to_its_directory() {
        let input = parse_input(INPUT, Path::new("work")).unwrap();
        assert_eq!(input.output, PathBuf::from("work/m1.msgpack"));
        assert_eq!(input.sigma0, vec![1.0e10, 10.0, 1.0]);
        assert_eq!(
            input.nuclides,
            vec![
                NuclideInput {
                    name: "H-1".into(),
                    endf: PathBuf::from("work/endf/h.endf"),
                    temperatures: vec![293.6, 600.0],
                    thermal: Some(PathBuf::from("work/tsl/h2o.endf")),
                },
                NuclideInput {
                    name: "U-238".into(),
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
        let replaced = |from: &str, to: &str| {
            assert!(INPUT.contains(from), "{from}");
            parse_input(&INPUT.replace(from, to), base)
        };
        // 必須のキーがない、知らないキー
        for line in [
            "output = \"m1.msgpack\"\n",
            "endf-directory = \"endf\"\n",
            "sigma0 = [1e10, 10, 1]\n",
        ] {
            assert!(replaced(line, "").is_err(), "without {line}");
        }
        assert!(replaced("output", "out").is_err());
        assert!(replaced("endf = \"u.endf\"", "endf = \"u.endf\"\nlibrary = \"x\"").is_err());
        // 核種がない、同じ名前の核種（TOML の重複キー）
        let header = INPUT.split("[nuclide.H-1]").next().unwrap();
        assert!(parse_input(header, base).is_err());
        assert!(replaced("[nuclide.U-238]", "[nuclide.H-1]").is_err());
        // thermal を使うのに thermal-directory がない
        assert!(replaced("thermal-directory = \"tsl\"\n", "").is_err());
        // σ0 が降順でない、正でない
        assert!(replaced("sigma0 = [1e10, 10, 1]", "sigma0 = [1, 10]").is_err());
        assert!(replaced("sigma0 = [1e10, 10, 1]", "sigma0 = [1e10, 0]").is_err());
        // 温度がない、正でない、重複している
        assert!(replaced("temperatures = [293.6, 900]", "temperatures = []").is_err());
        assert!(replaced("temperatures = [293.6, 900]", "temperatures = [-1]").is_err());
        assert!(replaced("temperatures = [293.6, 900]", "temperatures = [900, 900]").is_err());
        // TOML の文法の誤りは位置つきで報告する
        let err = replaced("sigma0 = [1e10, 10, 1]", "sigma0 = [1e10, 10").unwrap_err();
        assert!(err.to_string().contains("line"), "{err}");
    }
}
