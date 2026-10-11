//! ENDF ファイルから多群ライブラリを作るコマンド。
//!
//! 使い方: xs <入力ファイル>
//!         xs range <入力ファイル> <範囲のファイル>
//!         xs nuclide <入力ファイル> <範囲のファイル> <核種> <出力>
//!         xs assemble <入力ファイル> <範囲のファイル> <核種ごとの結果のディレクトリ>
//!
//! 引数が入力ファイルだけのときは、すべての核種を処理して入力の `output` にライブラリを書く。ほかの 3 つは、同じ結果を
//! 段階に分けて作る（ADR-0138）。`range` は、すべての核種の分離共鳴の上端から、ポイントワイズの範囲を求めて書く。
//! `nuclide` は、その範囲で核種 1 つだけを処理し、1 核種のライブラリを書く。`assemble` は、ディレクトリの
//! `<核種>.msgpack` を入力の核種の順に、入力と範囲のファイルと同じ条件で作られたことを確かめて、1 つのライブラリにして、
//! 入力の `output` に書く。
//!
//! 入力ファイルは TOML（パスは入力ファイルのディレクトリからの相対パス）:
//!
//! ```toml
//! output = "m1.msgpack"                       # 出力するライブラリ
//! endf-directory = "../endf/neutrons"         # 核種の ENDF ファイルのディレクトリ
//! thermal-directory = "../endf/thermal_scatt" # 熱中性子散乱のファイルのディレクトリ（thermal を使うとき）
//! sigma0 = [1e10, 1e5, 1e4, 1e3, 300, 100, 30, 10, 1]   # 背景断面積 (barn)、降順（最初が無限希釈）
//! scattering-order = 3                        # 散乱行列の Legendre の次数 L（省略すると 3。ADR-0076）
//!
//! [nuclide.H-1]
//! endf = "n-001_H_001.endf"
//! temperatures = [293.6, 600]
//! thermal = "tsl-HinH2O.endf"
//! sab-upper = "evaluated"     # S(α, β) を使う上限: "evaluated"（MF7 MT4 の B(4)）か eV（省略すると "evaluated"）
//! sab-outside = "sct"         # 表の外: "sct"（THERMR と同じ短時間衝突近似）か "zero"（省略すると "sct"）
//! free-gas-upper = 9118       # S(α, β) の上限からこの eV まで弾性散乱を自由ガスで（省略すると `DEFAULT_FREE_GAS_UPPER`）
//!
//! [nuclide.U-238]
//! endf = "n-092_U_238.endf"
//! temperatures = [293.6, 600, 900, 1200]
//! ```
//!
//! 核種ごとの処理時間はログとして標準エラーに出す（`RUST_LOG=debug` で温度・段階ごとの内訳）。
//!
//! `thermal` を指定した核種は、`sab-upper` 未満の弾性散乱を S(α, β) で、その上の `free-gas-upper` 未満を自由ガス
//! 模型で扱う（ADR-0133）。ほかの核種は 4 eV 未満を自由ガス模型で扱う。
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
use xlc_xs::multigroup::{SabOutside, Weight};
use xlc_xs::pointwise::DEFAULT_TOLERANCE;
use xlc_xs::process::{
    SabUpper, Settings, ThermalScattering, process_nuclide, read_thermal_scattering, resolved,
};

/// 入力ファイルを解釈した結果（パスは入力ファイルのディレクトリで解決済み）。
#[derive(Debug, Clone, PartialEq)]
struct Input {
    output: PathBuf,
    sigma0: Vec<f64>,
    scattering_order: usize,
    nuclides: Vec<NuclideInput>,
}

#[derive(Debug, Clone, PartialEq)]
struct NuclideInput {
    name: String,
    endf: PathBuf,
    temperatures: Vec<f64>,
    thermal: Option<ThermalInput>,
}

/// S(α, β) を持つ核種の熱領域の扱い（ADR-0133）。
#[derive(Debug, Clone, PartialEq)]
struct ThermalInput {
    path: PathBuf,
    upper: SabUpper,
    outside: SabOutside,
    free_gas_upper: f64,
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
    /// 省略すると `DEFAULT_SCATTERING_ORDER`
    scattering_order: Option<usize>,
    /// 名前ごとの核種
    nuclide: BTreeMap<String, NuclideEntry>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct NuclideEntry {
    endf: PathBuf,
    temperatures: Vec<f64>,
    /// 省略すると、熱領域の弾性散乱を自由ガス模型で扱う
    thermal: Option<PathBuf>,
    /// `thermal` を使う核種だけ。省略すると "evaluated"
    sab_upper: Option<SabUpperEntry>,
    /// `thermal` を使う核種だけ。省略すると "sct"
    sab_outside: Option<SabOutsideEntry>,
    /// `thermal` を使う核種だけ。省略すると `DEFAULT_FREE_GAS_UPPER`
    free_gas_upper: Option<f64>,
}

#[derive(Deserialize)]
#[serde(untagged)]
enum SabUpperEntry {
    Energy(f64),
    Name(String),
}

#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
enum SabOutsideEntry {
    Sct,
    Zero,
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
            (None, _) => {
                if entry.sab_upper.is_some()
                    || entry.sab_outside.is_some()
                    || entry.free_gas_upper.is_some()
                {
                    return invalid(format!(
                        "{name}: `sab-upper`, `sab-outside` and `free-gas-upper` need `thermal`"
                    ));
                }
                None
            }
            (Some(thermal), Some(directory)) => {
                let upper = match &entry.sab_upper {
                    None => SabUpper::Evaluated,
                    Some(SabUpperEntry::Name(n)) if n == "evaluated" => SabUpper::Evaluated,
                    Some(SabUpperEntry::Energy(e)) if positive(e) => SabUpper::Energy(*e),
                    Some(_) => {
                        return invalid(format!(
                            "{name}: `sab-upper` must be \"evaluated\" or a positive energy (eV)"
                        ));
                    }
                };
                let outside = match entry.sab_outside {
                    None | Some(SabOutsideEntry::Sct) => SabOutside::ShortCollisionTime,
                    Some(SabOutsideEntry::Zero) => SabOutside::Zero,
                };
                let free_gas_upper = entry.free_gas_upper.unwrap_or(DEFAULT_FREE_GAS_UPPER);
                if !positive(&free_gas_upper) {
                    return invalid(format!(
                        "{name}: `free-gas-upper` must be a positive energy (eV), got {free_gas_upper}"
                    ));
                }
                Some(ThermalInput {
                    path: base.join(directory).join(thermal),
                    upper,
                    outside,
                    free_gas_upper,
                })
            }
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
        scattering_order: file.scattering_order.unwrap_or(DEFAULT_SCATTERING_ORDER),
        nuclides,
    })
}

const USAGE: &str = "usage: xs <input file>
       xs range <input file> <range file>
       xs nuclide <input file> <range file> <nuclide> <output>
       xs assemble <input file> <range file> <nuclide directory>";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result =
        init_logging().and_then(
            |()| match args.iter().map(String::as_str).collect::<Vec<_>>()[..] {
                [path] => read_input(Path::new(path)).and_then(|input| run(&input)),
                ["range", path, out] => {
                    read_input(Path::new(path)).and_then(|input| run_range(&input, Path::new(out)))
                }
                ["nuclide", path, range, name, out] => read_input(Path::new(path))
                    .and_then(|input| run_nuclide(&input, Path::new(range), name, Path::new(out))),
                ["assemble", path, range, directory] => read_input(Path::new(path))
                    .and_then(|input| run_assemble(&input, Path::new(range), Path::new(directory))),
                _ => Err(anyhow!(USAGE)),
            },
        );
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
#[allow(
    clippy::disallowed_methods,
    reason = "RUST_LOG は表示する量だけを選び、計算の結果を変えない（ADR-0137）"
)]
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

/// 入力の核種の ENDF のテープを、並列に読む。
fn read_tapes(nuclides: &[NuclideInput]) -> anyhow::Result<Vec<Tape>> {
    let tapes: Vec<anyhow::Result<Tape>> = std::thread::scope(|scope| {
        let handles: Vec<_> = nuclides
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
    tapes.into_iter().collect()
}

fn settings_of(input: &Input, pointwise_range: Option<(f64, f64)>) -> Settings {
    Settings {
        groups: GroupStructure::wims69(),
        weight: Weight::lwr(),
        sigma0: input.sigma0.clone(),
        thermal_cutoff: THERMAL_CUTOFF,
        tolerance: DEFAULT_TOLERANCE,
        pointwise_range,
        scattering_order: input.scattering_order,
    }
}

fn library_of(settings: Settings, nuclides: Vec<NuclideData>) -> Library {
    Library {
        groups: settings.groups,
        sigma0: settings.sigma0,
        pointwise_range: settings.pointwise_range.map(|(low, high)| [low, high]),
        nuclides,
    }
}

/// ファイルを書く。DVC は段階の前に出力を消すとき、空になったディレクトリも消すので、ディレクトリを作る。
fn write_file(path: &Path, contents: &[u8]) -> anyhow::Result<()> {
    let context = || format!("{}", path.display());
    if let Some(directory) = path.parent() {
        std::fs::create_dir_all(directory).with_context(context)?;
    }
    std::fs::write(path, contents).with_context(context)
}

fn write_library(library: &Library, path: &Path) -> anyhow::Result<()> {
    write_file(path, &library.to_msgpack()?)
}

/// 入力のすべての核種を処理して、ライブラリを書く。
fn run(input: &Input) -> anyhow::Result<()> {
    // テープの読み込みと、核種ごとの処理は独立なので並列に行う
    let tapes = read_tapes(&input.nuclides)?;
    let pointwise_range = pointwise_range(&input.nuclides, &tapes, &GroupStructure::wims69())?;
    let settings = settings_of(input, pointwise_range);
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
    write_library(&library_of(settings, nuclides), &input.output)
}

/// すべての核種の分離共鳴の上端から、ポイントワイズの範囲を求めて書く（`nuclide` の段階の入力になる）。
fn run_range(input: &Input, out: &Path) -> anyhow::Result<()> {
    let tapes = read_tapes(&input.nuclides)?;
    let range = pointwise_range(&input.nuclides, &tapes, &GroupStructure::wims69())?;
    write_file(out, range_to_toml(range).as_bytes())
}

fn read_range(path: &Path) -> anyhow::Result<Option<(f64, f64)>> {
    let context = || format!("{}", path.display());
    let text = std::fs::read_to_string(path).with_context(context)?;
    range_from_toml(&text).with_context(context)
}

/// 核種 `name` だけを、`range` のファイルのポイントワイズの範囲で処理して、1 核種のライブラリを書く。
fn run_nuclide(input: &Input, range: &Path, name: &str, out: &Path) -> anyhow::Result<()> {
    let nuclide = input
        .nuclides
        .iter()
        .find(|n| n.name == name)
        .ok_or_else(|| anyhow!("{name} is not in the input"))?;
    let pointwise_range = read_range(range)?;
    let settings = settings_of(input, pointwise_range);
    let data = process_file(nuclide, &read_tape(&nuclide.endf)?, &settings)?;
    write_library(&library_of(settings, vec![data]), out)
}

/// `directory` の `<核種>.msgpack`（`nuclide` の結果）を、入力の核種の順に 1 つのライブラリにして書く。
/// 各結果が入力と `range` のファイルの条件で作られたことを確かめる。
fn run_assemble(input: &Input, range: &Path, directory: &Path) -> anyhow::Result<()> {
    let range = read_range(range)?;
    let mut parts = Vec::with_capacity(input.nuclides.len());
    for nuclide in &input.nuclides {
        let path = directory.join(format!("{}.msgpack", nuclide.name));
        let bytes = std::fs::read(&path).with_context(|| format!("{}", path.display()))?;
        let part = Library::from_msgpack(&bytes).with_context(|| format!("{}", path.display()))?;
        check_part(nuclide, &part)
            .and_then(|()| check_shared(input, range, &part))
            .with_context(|| format!("{}", path.display()))?;
        parts.push(part);
    }
    write_library(&Library::merge(parts)?, &input.output)
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

/// `range_to_toml` と `range_from_toml` が書く、ポイントワイズの範囲のファイルの中身。
#[derive(Deserialize)]
#[serde(deny_unknown_fields, rename_all = "kebab-case")]
struct RangeFile {
    /// 範囲がなければ空の配列
    pointwise_range: Vec<f64>,
}

/// ポイントワイズの範囲 (eV) を TOML に書く。範囲がなければ空の配列を書く（項目は欠かさない。ADR-0129）。
/// 浮動小数点数は、読み直して同じ値になる最短の表記（`{:?}`）で書く。
fn range_to_toml(range: Option<(f64, f64)>) -> String {
    match range {
        Some((low, high)) => format!("pointwise-range = [{low:?}, {high:?}]\n"),
        None => "pointwise-range = []\n".into(),
    }
}

/// `range_to_toml` の形を読む。
fn range_from_toml(text: &str) -> anyhow::Result<Option<(f64, f64)>> {
    let file: RangeFile = toml::from_str(text)?;
    match file.pointwise_range[..] {
        [] => Ok(None),
        [low, high] if low > 0.0 && high > low && high.is_finite() => Ok(Some((low, high))),
        ref other => Err(anyhow!(
            "pointwise-range must be empty or [low, high] with 0 < low < high, got {other:?}"
        )),
    }
}

/// 核種ごとの結果が、入力のその核種の結果であること（1 核種だけで、名前と温度が入力と同じ）を確かめる。
fn check_part(nuclide: &NuclideInput, part: &Library) -> anyhow::Result<()> {
    let [data] = &part.nuclides[..] else {
        return Err(anyhow!(
            "{}: expected one nuclide, found {}",
            nuclide.name,
            part.nuclides.len()
        ));
    };
    if data.name != nuclide.name {
        return Err(anyhow!("{}: the result is for {}", nuclide.name, data.name));
    }
    let temperatures: Vec<f64> = data.temperatures.iter().map(|t| t.temperature).collect();
    if temperatures != nuclide.temperatures {
        return Err(anyhow!(
            "{}: the result has temperatures {temperatures:?}, the input has {:?}",
            nuclide.name,
            nuclide.temperatures
        ));
    }
    Ok(())
}

/// 核種ごとの結果が、入力と範囲のファイルと同じ条件で作られたこと（群構造、σ0、ポイントワイズの範囲、散乱行列の次数）を
/// 確かめる。古い結果が残っていても、組み立てが黙って食い違うライブラリを作らないため。
fn check_shared(input: &Input, range: Option<(f64, f64)>, part: &Library) -> anyhow::Result<()> {
    if part.groups != GroupStructure::wims69() {
        return Err(anyhow!("the group structure is not WIMS 69"));
    }
    if part.sigma0 != input.sigma0 {
        return Err(anyhow!(
            "sigma0 is {:?}, the input has {:?}",
            part.sigma0,
            input.sigma0
        ));
    }
    let range = range.map(|(low, high)| [low, high]);
    if part.pointwise_range != range {
        return Err(anyhow!(
            "the pointwise range is {:?}, the range file has {range:?}",
            part.pointwise_range
        ));
    }
    let orders = input.scattering_order + 1;
    for data in &part.nuclides {
        for t in &data.temperatures {
            if t.scatter.len() != orders {
                return Err(anyhow!(
                    "{} at {} K has {} scattering orders, the input has {orders}",
                    data.name,
                    t.temperature,
                    t.scatter.len()
                ));
            }
        }
    }
    Ok(())
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
        Some(input) => Some(ThermalScattering {
            data: read_thermal_scattering(&read_tape(&input.path)?)
                .with_context(|| name.clone())?,
            upper: input.upper,
            outside: input.outside,
            free_gas_upper: input.free_gas_upper,
        }),
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

/// S(α, β) を持たない核種で、この入射エネルギー (eV) 未満の弾性散乱を自由ガス模型で扱う（WIMS 69群の熱群の上限）。
/// ポイントワイズの断面積（UFG）の範囲の下限もこれで決める
const THERMAL_CUTOFF: f64 = 4.0;

/// 入力で省略したときの、S(α, β) を持つ核種の弾性散乱を自由ガス模型で扱う上限 (eV)（ADR-0133）
const DEFAULT_FREE_GAS_UPPER: f64 = 9118.0;

/// 入力で省略したときの、散乱行列のルジャンドルのモーメントの次数（ADR-0076）。MOC の非等方散乱は、格子の入力で
/// ライブラリにあるこの次数までを選ぶ
const DEFAULT_SCATTERING_ORDER: usize = 3;

#[cfg(test)]
mod tests {
    use xlc_format::library::TemperatureData;

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
        assert_eq!(input.scattering_order, 3);
        assert_eq!(
            input.nuclides,
            vec![
                NuclideInput {
                    name: "H-1".into(),
                    endf: PathBuf::from("work/endf/h.endf"),
                    temperatures: vec![293.6, 600.0],
                    thermal: Some(ThermalInput {
                        path: PathBuf::from("work/tsl/h2o.endf"),
                        upper: SabUpper::Evaluated,
                        outside: SabOutside::ShortCollisionTime,
                        free_gas_upper: DEFAULT_FREE_GAS_UPPER,
                    }),
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
    fn accepts_a_scattering_order() {
        let text = INPUT.replace("sigma0 =", "scattering-order = 1\nsigma0 =");
        let input = parse_input(&text, Path::new(".")).unwrap();
        assert_eq!(input.scattering_order, 1);
        let text = INPUT.replace("sigma0 =", "scattering-order = 0\nsigma0 =");
        assert_eq!(
            parse_input(&text, Path::new(".")).unwrap().scattering_order,
            0
        );
    }

    #[test]
    fn accepts_the_thermal_scattering_limits() {
        let text = INPUT.replace(
            "thermal = \"h2o.endf\"",
            "thermal = \"h2o.endf\"\nsab-upper = 4\nsab-outside = \"zero\"\nfree-gas-upper = 48.052",
        );
        let input = parse_input(&text, Path::new(".")).unwrap();
        let thermal = input.nuclides[0].thermal.as_ref().unwrap();
        assert_eq!(thermal.upper, SabUpper::Energy(4.0));
        assert_eq!(thermal.outside, SabOutside::Zero);
        assert_eq!(thermal.free_gas_upper, 48.052);
        let text = INPUT.replace(
            "thermal = \"h2o.endf\"",
            "thermal = \"h2o.endf\"\nsab-upper = \"evaluated\"\nsab-outside = \"sct\"",
        );
        let input = parse_input(&text, Path::new(".")).unwrap();
        let thermal = input.nuclides[0].thermal.as_ref().unwrap();
        assert_eq!(thermal.upper, SabUpper::Evaluated);
        assert_eq!(thermal.outside, SabOutside::ShortCollisionTime);
    }

    #[test]
    fn rejects_bad_thermal_scattering_limits() {
        let with = |lines: &str| {
            let text = INPUT.replace(
                "thermal = \"h2o.endf\"",
                &format!("thermal = \"h2o.endf\"\n{lines}"),
            );
            parse_input(&text, Path::new("."))
        };
        for lines in [
            "sab-upper = \"top\"",
            "sab-upper = 0",
            "sab-upper = -4",
            "sab-outside = \"free-gas\"",
            "free-gas-upper = 0",
            "free-gas-upper = \"all\"",
        ] {
            assert!(with(lines).is_err(), "{lines}");
        }
        // S(α, β) を使わない核種には書けない
        for key in [
            "sab-upper = 4",
            "sab-outside = \"sct\"",
            "free-gas-upper = 100",
        ] {
            let text = INPUT.replace(
                "temperatures = [293.6, 900]",
                &format!("temperatures = [293.6, 900]\n{key}"),
            );
            assert!(parse_input(&text, Path::new(".")).is_err(), "{key}");
        }
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
        // 散乱行列の次数が負、整数でない、文字列
        for value in ["-1", "1.5", "\"3\""] {
            let text = format!("scattering-order = {value}\n{INPUT}");
            assert!(parse_input(&text, base).is_err(), "{value}");
        }
        // 温度がない、正でない、重複している
        assert!(replaced("temperatures = [293.6, 900]", "temperatures = []").is_err());
        assert!(replaced("temperatures = [293.6, 900]", "temperatures = [-1]").is_err());
        assert!(replaced("temperatures = [293.6, 900]", "temperatures = [900, 900]").is_err());
        // TOML の文法の誤りは位置つきで報告する
        let err = replaced("sigma0 = [1e10, 10, 1]", "sigma0 = [1e10, 10").unwrap_err();
        assert!(err.to_string().contains("line"), "{err}");
    }

    #[test]
    fn the_pointwise_range_file_round_trips() {
        for range in [Some((4.0, 9118.0)), Some((0.1 + 0.2, 1.0e5)), None] {
            let text = range_to_toml(range);
            assert!(text.contains("pointwise-range"), "{text}");
            assert_eq!(range_from_toml(&text).unwrap(), range);
        }
    }

    #[test]
    fn rejects_a_bad_pointwise_range_file() {
        for text in [
            "",
            "pointwise-range = [4.0]",
            "pointwise-range = [4.0, 5.0, 6.0]",
            "pointwise-range = [5.0, 4.0]",
            "pointwise-range = [4.0, 4.0]",
            "pointwise-range = [0.0, 4.0]",
            "pointwise-range = [4.0, 1e5]\nextra = 1",
            "pointwise-range = \"4 1e5\"",
        ] {
            assert!(range_from_toml(text).is_err(), "{text}");
        }
    }

    fn part_of(name: &str, temperatures: &[f64]) -> Library {
        Library {
            groups: GroupStructure::wims69(),
            sigma0: vec![1.0e10],
            pointwise_range: None,
            nuclides: vec![NuclideData {
                name: name.into(),
                za: 0.0,
                awr: 1.0,
                resolved: None,
                thermal_scattering: false,
                temperatures: temperatures
                    .iter()
                    .map(|&temperature| TemperatureData {
                        temperature,
                        chi: Vec::new(),
                        scatter: Vec::new(),
                        shielded: Vec::new(),
                        pointwise: Default::default(),
                    })
                    .collect(),
            }],
        }
    }

    #[test]
    fn a_nuclide_result_must_match_its_input() {
        let input = parse_input(INPUT, Path::new(".")).unwrap();
        let uranium = &input.nuclides[1];
        assert!(check_part(uranium, &part_of("U-238", &[293.6, 900.0])).is_ok());
        // 名前、温度の値と順、核種の数が入力と違う
        for part in [
            part_of("U-235", &[293.6, 900.0]),
            part_of("U-238", &[293.6]),
            part_of("U-238", &[900.0, 293.6]),
            part_of("U-238", &[293.6, 600.0]),
        ] {
            assert!(check_part(uranium, &part).is_err());
        }
        let mut two = part_of("U-238", &[293.6, 900.0]);
        two.nuclides.extend(part_of("H-1", &[293.6]).nuclides);
        assert!(check_part(uranium, &two).is_err());
        let mut none = part_of("U-238", &[293.6, 900.0]);
        none.nuclides.clear();
        assert!(check_part(uranium, &none).is_err());
    }

    #[test]
    fn a_nuclide_result_must_share_the_conditions_of_the_input() {
        let input = parse_input(INPUT, Path::new(".")).unwrap();
        let range = Some((4.0, 9118.0));
        let mut part = part_of("U-238", &[293.6, 900.0]);
        part.sigma0 = input.sigma0.clone();
        part.pointwise_range = range.map(|(low, high)| [low, high]);
        for t in &mut part.nuclides[0].temperatures {
            t.scatter = vec![Vec::new(); input.scattering_order + 1];
        }
        assert!(check_shared(&input, range, &part).is_ok());
        let mut sigma0 = part.clone();
        sigma0.sigma0 = vec![1.0e10, 5.0, 1.0];
        let mut groups = part.clone();
        groups.groups = GroupStructure::new(vec![1.0e7, 1.0, 1.0e-5]).unwrap();
        let mut order = part.clone();
        order.nuclides[0].temperatures[1].scatter.pop();
        let mut other_range = part.clone();
        other_range.pointwise_range = Some([4.0, 1.0e3]);
        let mut no_range = part.clone();
        no_range.pointwise_range = None;
        for stale in [sigma0, groups, order, other_range, no_range] {
            assert!(check_shared(&input, range, &stale).is_err());
        }
        assert!(check_shared(&input, None, &part).is_err());
    }
}
