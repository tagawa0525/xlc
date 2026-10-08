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
//! 核種ごとの処理時間はログとして標準エラーに出す（`RUST_LOG=debug` で温度・段階ごとの内訳）。
//!
//! `thermal` を指定した核種は、4 eV 未満の弾性散乱を S(α, β) で、ほかは自由ガス模型で扱う。
//! 群構造は WIMS 69群、重み関数は軽水炉向けの iwt=4 とする。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

use std::collections::BTreeMap;
use std::env::VarError;
use std::io::IsTerminal;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{Context, anyhow};
use serde::Deserialize;
use tracing_subscriber::EnvFilter;
use tracing_subscriber::filter::LevelFilter;
use tracing_subscriber::fmt::format::FmtSpan;
use xlc_format::group::GroupStructure;
use xlc_format::library::{Library, NuclideData, is_nuclide_name};
use xlc_xs::endf::Tape;
use xlc_xs::multigroup::Weight;
use xlc_xs::pointwise::DEFAULT_TOLERANCE;
use xlc_xs::process::{Settings, process_nuclide, read_thermal_scattering, resolved};

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
    #[error(transparent)]
    Toml(#[from] toml::de::Error),
    #[error("{0}")]
    Invalid(String),
}

/// 入力ファイルの中身（パスは未解決）。
#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct InputFile {
    output: PathBuf,
    endf_directory: PathBuf,
    /// 熱中性子散乱を使う核種がなければ省略できる
    thermal_directory: Option<PathBuf>,
    sigma0: Vec<f64>,
    /// 名前ごとの核種
    nuclide: BTreeMap<String, NuclideEntry>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct NuclideEntry {
    endf: PathBuf,
    temperatures: Vec<f64>,
    /// 省略すると、熱領域の弾性散乱を自由ガス模型で扱う
    thermal: Option<PathBuf>,
}

/// 入力ファイルを読む（パスは `base` で解決する）。
fn parse_input(text: &str, base: &Path) -> Result<Input, InputError> {
    let file: InputFile = toml::from_str(text)?;
    let invalid = |message: String| Err(InputError::Invalid(message));
    let positive = |v: &f64| *v > 0.0 && v.is_finite();
    if file.sigma0.is_empty()
        || !file.sigma0.iter().all(positive)
        || file.sigma0.windows(2).any(|w| w[1] >= w[0])
    {
        return invalid(format!(
            "sigma0 must be positive and descending, got {:?}",
            file.sigma0
        ));
    }
    if file.nuclide.is_empty() {
        return invalid("no nuclide is given".into());
    }
    let endf_directory = base.join(&file.endf_directory);
    let mut nuclides = Vec::with_capacity(file.nuclide.len());
    for (name, entry) in file.nuclide {
        if !is_nuclide_name(&name) {
            return invalid(format!(
                "nuclide name `{name}` is not of the form U-235, Am-242m"
            ));
        }
        let t = &entry.temperatures;
        if t.is_empty() || !t.iter().all(positive) {
            return invalid(format!("{name}: temperatures must be positive, got {t:?}"));
        }
        if t.iter().enumerate().any(|(i, x)| t[..i].contains(x)) {
            return invalid(format!("{name}: a temperature is repeated in {t:?}"));
        }
        let thermal = match (&entry.thermal, &file.thermal_directory) {
            (None, _) => None,
            (Some(thermal), Some(directory)) => Some(base.join(directory).join(thermal)),
            (Some(_), None) => {
                return invalid(format!("{name}: `thermal` needs `thermal-directory`"));
            }
        };
        nuclides.push(NuclideInput {
            endf: endf_directory.join(&entry.endf),
            temperatures: entry.temperatures,
            thermal,
            name,
        });
    }
    Ok(Input {
        output: base.join(file.output),
        sigma0: file.sigma0,
        nuclides,
    })
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.as_slice() {
        [path] => init_logging()
            .and_then(|()| read_input(Path::new(path)))
            .and_then(|input| run(&input)),
        _ => Err(anyhow!("usage: xs <input file>")),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error:#}");
            ExitCode::FAILURE
        }
    }
}

/// ログ（tracing）を標準エラーに出す。表示する量は環境変数 `RUST_LOG` で選ぶ（なければ info）。
/// 区間（span）が閉じるときに、その経過時間も出す。
fn init_logging() -> anyhow::Result<()> {
    // EnvFilter::from_env は UTF-8 でない値を黙って空とみなすので、自分で読む
    let directives = match std::env::var("RUST_LOG") {
        Ok(value) => value,
        Err(VarError::NotPresent) => String::new(),
        Err(error @ VarError::NotUnicode(_)) => return Err(anyhow!("RUST_LOG: {error}")),
    };
    // ParseError の表示は原因の文を含むので、原因の連鎖は付けない
    let filter = EnvFilter::builder()
        .with_default_directive(LevelFilter::INFO.into())
        .parse(directives)
        .map_err(|error| anyhow!("RUST_LOG: {error}"))?;
    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_writer(std::io::stderr)
        .with_ansi(std::io::stderr().is_terminal())
        .with_span_events(FmtSpan::CLOSE)
        .init();
    Ok(())
}

fn read_input(path: &Path) -> anyhow::Result<Input> {
    let text = std::fs::read_to_string(path).with_context(|| format!("{}", path.display()))?;
    let base = path.parent().unwrap_or(Path::new("."));
    parse_input(&text, base).with_context(|| format!("{}", path.display()))
}

fn run(input: &Input) -> anyhow::Result<()> {
    let groups = GroupStructure::wims69();
    // テープの読み込みと、核種ごとの処理は独立なので並列に行う
    let tapes: Vec<anyhow::Result<Tape>> = std::thread::scope(|scope| {
        let handles: Vec<_> = input
            .nuclides
            .iter()
            .map(|nuclide| scope.spawn(|| read_tape(&nuclide.endf)))
            .collect();
        handles
            .into_iter()
            .map(|h| {
                h.join()
                    .unwrap_or_else(|_| Err(anyhow!("reading panicked")))
            })
            .collect()
    });
    let tapes = tapes.into_iter().collect::<anyhow::Result<Vec<_>>>()?;
    let pointwise_range = pointwise_range(&input.nuclides, &tapes, &groups)?;
    let settings = Settings {
        groups,
        weight: Weight::lwr(),
        sigma0: input.sigma0.clone(),
        thermal_cutoff: THERMAL_CUTOFF,
        tolerance: DEFAULT_TOLERANCE,
        pointwise_range,
        scattering_order: SCATTERING_ORDER,
    };
    let results: Vec<anyhow::Result<NuclideData>> = std::thread::scope(|scope| {
        let handles: Vec<_> = input
            .nuclides
            .iter()
            .zip(&tapes)
            .map(|(nuclide, tape)| scope.spawn(|| process_file(nuclide, tape, &settings)))
            .collect();
        handles
            .into_iter()
            .map(|h| {
                h.join()
                    .unwrap_or_else(|_| Err(anyhow!("processing panicked")))
            })
            .collect()
    });
    let nuclides = results.into_iter().collect::<anyhow::Result<Vec<_>>>()?;
    let library = Library {
        groups: settings.groups,
        sigma0: settings.sigma0,
        pointwise_range: pointwise_range.map(|(low, high)| [low, high]),
        nuclides,
    };
    std::fs::write(&input.output, library.to_msgpack()?)
        .with_context(|| format!("{}", input.output.display()))
}

fn read_tape(path: &Path) -> anyhow::Result<Tape> {
    let context = || format!("{}", path.display());
    let text = std::fs::read_to_string(path).with_context(context)?;
    Tape::parse(&text).with_context(context)
}

/// ライブラリに持たせるポイントワイズの断面積の範囲 (eV): 熱中性子散乱の打ち切り以上で最も低い群の境界から、全核種の
/// 分離共鳴の上端のうち最も高いもの以下で最も高い群の境界まで（どのセルの超微細群の範囲もこれに収まる）。
/// 分離共鳴を持つ核種がなければ、超微細群の計算をしないので範囲はない（`None`）。
fn pointwise_range(
    nuclides: &[NuclideInput],
    tapes: &[Tape],
    groups: &GroupStructure,
) -> anyhow::Result<Option<(f64, f64)>> {
    let mut top: Option<f64> = None;
    for (nuclide, tape) in nuclides.iter().zip(tapes) {
        let materials = tape.materials();
        let [mat] = materials[..] else {
            return Err(anyhow!(
                "{}: expected one material, found {materials:?}",
                nuclide.name
            ));
        };
        if let Some(r) = resolved(tape, mat).with_context(|| nuclide.name.clone())? {
            top = Some(top.map_or(r.top, |u| u.max(r.top)));
        }
    }
    let Some(top) = top else {
        return Ok(None);
    };
    let bounds = groups.bounds();
    let low = bounds.iter().rev().find(|&&b| b >= THERMAL_CUTOFF);
    let high = bounds.iter().find(|&&b| b <= top);
    match (low, high) {
        (Some(&low), Some(&high)) if high > low => Ok(Some((low, high))),
        _ => Err(anyhow!(
            "no group bounds between the thermal cutoff {THERMAL_CUTOFF} eV and {top} eV"
        )),
    }
}

/// 核種を処理する（処理の時間をログに出す）。
#[tracing::instrument(name = "nuclide", skip_all, fields(name = %nuclide.name))]
fn process_file(
    nuclide: &NuclideInput,
    tape: &Tape,
    settings: &Settings,
) -> anyhow::Result<NuclideData> {
    let name = &nuclide.name;
    let thermal = match &nuclide.thermal {
        Some(path) => {
            Some(read_thermal_scattering(&read_tape(path)?).with_context(|| name.clone())?)
        }
        None => None,
    };
    let data = process_nuclide(
        tape,
        name,
        &nuclide.temperatures,
        thermal.as_ref(),
        settings,
    )
    .with_context(|| name.clone())?;
    Ok(data)
}

/// この入射エネルギー (eV) 未満の弾性散乱を熱領域の扱いにする（WIMS 69群の熱群の上限）
const THERMAL_CUTOFF: f64 = 4.0;

/// 散乱行列のルジャンドルのモーメントの次数（MOC の非等方散乱は格子の入力でこの次数まで選ぶ）
const SCATTERING_ORDER: usize = 3;

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
    fn no_resolved_resonances_means_no_pointwise_range() {
        // 分離共鳴を持つ核種がなければ、超微細群の計算に使うポイントワイズの断面積は要らない（エラーにしない）
        assert!(matches!(
            pointwise_range(&[], &[], &GroupStructure::wims69()),
            Ok(None)
        ));
    }

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
        // 核種名が約束事の形でない
        assert!(replaced("[nuclide.U-238]", "[nuclide.U238]").is_err());
        assert!(replaced("[nuclide.U-238]", "[nuclide.\"\"]").is_err());
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
