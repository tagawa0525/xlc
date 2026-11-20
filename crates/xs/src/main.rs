//! ENDF ファイルから多群ライブラリを作るコマンド。
//!
//! 使い方: xs [--threads <数>] <入力ファイル>
//!         xs [--threads <数>] range <入力ファイル> <範囲のファイル>
//!         xs [--threads <数>] nuclide <入力ファイル> <範囲のファイル> <核種> <出力>
//!         xs [--threads <数>] assemble <入力ファイル> <範囲のファイル> <核種ごとの結果のディレクトリ>
//!         xs [--threads <数>] check <報告> <ENDF のファイル>...
//!
//! `--threads` は同時に使うスレッドの数で、省略すると物理コアの数をこのプロセスに許された CPU の数で抑えたもの
//! （ADR-0067）。結果はスレッドの数によらない。
//!
//! 引数が入力ファイルだけのときは、すべての核種を処理して入力の `output` にライブラリを書く。ほかの 3 つは、同じ結果を
//! ステージに分けて作る（ADR-0138）。`range` は、すべての核種の分離共鳴の上端から、ポイントワイズの範囲を求めて書く。
//! `nuclide` は、その範囲で核種 1 つだけを処理し、1 核種のライブラリを書く。`assemble` は、ディレクトリの
//! `<核種>.msgpack` を入力の核種の順に、入力と範囲のファイルと同じ条件で作られたことを確かめて、1 つのライブラリにして、
//! 入力の `output` に書く。
//!
//! `check` は入力ファイルを使わず、ENDF のファイルごとに、輸送する核種として `nuclide` が読むものの読み込みと形式の検査だけを
//! 行い、止まった理由をファイルごとに報告（TOML の `[[evaluation]]` をパスの順に並べ、`endf` に与えたパス、止まれば `error` に
//! 理由）に書く（ADR-0266）。止まったファイルがあっても成功で終わる。
//!
//! 入力ファイルは TOML（パスは入力ファイルのディレクトリからの相対パス）:
//!
//! ```toml
//! output = "m1.msgpack"                       # 出力するライブラリ
//! endf-directory = "../endf/neutrons"         # 核種の ENDF ファイルのディレクトリ
//! thermal-directory = "../endf/thermal_scatt" # 熱中性子散乱のファイルのディレクトリ（thermal を使うとき）
//! sigma0 = [1e10, 1e5, 1e4, 1e3, 300, 100, 30, 10, 1]   # バックグラウンド断面積 (barn)、降順（最初が無限希釈）
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
//!
//! [nuclide.Cs-134]
//! endf = "n-055_Cs_134.endf"
//! temperatures = [600, 900]
//! activation = true           # 放射化の核種（輸送しない。省略すると false）
//!
//! [nuclide.Nd-144]
//! endf = "n-060_Nd_144.endf"
//! temperatures = [600, 900]
//! pointwise = false           # ポイントワイズの断面積を持たせない（UFG で解かない核種。省略すると true）
//!
//! [nuclide.Rb-87]
//! endf = "n-037_Rb_087.endf"
//! temperatures = [600]
//! # 評価の誤りを直す宣言（省略すると直さない）: ファイルの line 行目の MF/MT の行の field 番目（1〜6）の整数の欄を from から to に
//! corrections = [{ line = 1194, mf = 5, mt = 17, field = 1, from = 1, to = 3 }]
//! ```
//!
//! `activation = true` の核種は、チェーンの反応と核分裂の無限希釈の群の断面積だけを作り、ライブラリの `activation` に置く
//! （ADR-0238。輸送しないので `thermal` を持てず、ポイントワイズの範囲を決めない）。ライブラリは群速度の逆数も持つ。
//!
//! `pointwise = false` の核種は、ライブラリにポイントワイズの断面積を持たせない（ADR-0269。格子計算の `ufg-nuclides` に挙げない
//! 中心の物質の核種にだけ使える）。`activation` と `thermal` とは組めない。ポイントワイズの範囲は、持たせない核種も含めて決める。
//!
//! 核種ごとの処理時間はログとして標準エラーに出す（`RUST_LOG=debug` で温度・段階ごとの内訳）。
//!
//! `corrections` の宣言は、テープを読んだ直後、MF ごとの読み込みの前に、行の MF と MT と欄の元の値を確かめて当て、当てるたびに
//! warn を出す（ADR-0272）。どれかが宣言と違えば止まる。
//!
//! `thermal` を指定した核種は、`sab-upper` 未満の弾性散乱を S(α, β) で、その上の `free-gas-upper` 未満を自由ガス
//! 模型で扱う（ADR-0133）。ほかの核種は 4 eV 未満を自由ガス模型で扱う。
//! 群構造は WIMS 69群、重み関数は軽水炉向けの iwt=4 とする。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

use std::collections::BTreeMap;
use std::env::VarError;
use std::io::IsTerminal;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{Context, anyhow};
use serde::Deserialize;
use tracing_subscriber::EnvFilter;
use tracing_subscriber::filter::LevelFilter;
use tracing_subscriber::fmt::format::FmtSpan;
use xlc_format::group::GroupStructure;
use xlc_format::library::{ActivationData, Library, NuclideData, is_nuclide_name};
use xlc_xs::check::{Checked, check_contents, report_toml};
use xlc_xs::endf::{Correction, Tape};
use xlc_xs::multigroup::{SabOutside, Weight, inverse_velocity};
use xlc_xs::parallel::{in_parallel, set_threads};
use xlc_xs::pointwise::DEFAULT_TOLERANCE;
use xlc_xs::process::{
    SabUpper, Settings, ThermalScattering, process_activation, process_nuclide,
    read_thermal_scattering, resolved,
};
use xlc_xs::range;

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
    /// 放射化の核種（輸送しない。ADR-0238）
    activation: bool,
    /// ライブラリにポイントワイズの断面積を持たせる（false は UFG で解かない核種。ADR-0269）
    pointwise: bool,
    /// 読む前に当てる評価の誤りの修正（ADR-0272）
    corrections: Vec<Correction>,
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
    /// 放射化の核種（ADR-0238）。省略すると輸送する核種
    #[serde(default)]
    activation: bool,
    /// ポイントワイズの断面積を持たせるか（ADR-0269）。省略すると持たせる
    #[serde(default = "keeps_pointwise")]
    pointwise: bool,
    /// 評価の誤りの修正（ADR-0272）。省略すると直さない
    #[serde(default)]
    corrections: Vec<Correction>,
}

/// `pointwise` を省略したときの値（ADR-0269 の既定。今のとおり持たせる）
fn keeps_pointwise() -> bool {
    true
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
        if entry.activation && thermal.is_some() {
            return invalid(format!(
                "{name}: activation nuclides are not transported and take no `thermal`"
            ));
        }
        // 放射化の核種はもともと持たず、S(α, β) の核種は減速材の核種で UFG で解く（ADR-0268、ADR-0269）
        if !entry.pointwise && (entry.activation || thermal.is_some()) {
            return invalid(format!(
                "{name}: `pointwise = false` is for transported nuclides without `thermal` (activation nuclides have no \
                 pointwise data, and S(α, β) nuclides are moderator nuclides the UFG solves)"
            ));
        }
        nuclides.push(NuclideInput {
            endf: endf_directory.join(&entry.endf),
            temperatures: entry.temperatures,
            activation: entry.activation,
            pointwise: entry.pointwise,
            corrections: entry.corrections,
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

const USAGE: &str = "usage: xs [--threads <count>] <input file>
       xs [--threads <count>] range <input file> <range file>
       xs [--threads <count>] nuclide <input file> <range file> <nuclide> <output>
       xs [--threads <count>] assemble <input file> <range file> <nuclide directory>
       xs [--threads <count>] check <report> <ENDF file>...";

/// 引数の先頭の `--threads <数>` を読み、残りの引数と分ける。
fn thread_option(args: &[String]) -> anyhow::Result<(Option<NonZeroUsize>, &[String])> {
    match args {
        [flag, rest @ ..] if flag == "--threads" => match rest {
            [count, rest @ ..] => {
                let threads = count.parse().map_err(|error| {
                    anyhow!("--threads must be a positive integer, got `{count}`: {error}")
                })?;
                Ok((Some(threads), rest))
            }
            [] => Err(anyhow!("--threads needs the number of threads")),
        },
        _ => Ok((None, args)),
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = init_logging()
        .and_then(|()| thread_option(&args))
        .and_then(|(threads, args)| {
            if let Some(threads) = threads {
                set_threads(threads)?;
            }
            Ok(args)
        })
        .and_then(
            |args| match args.iter().map(String::as_str).collect::<Vec<_>>()[..] {
                [path] => read_input(Path::new(path)).and_then(|input| run(&input)),
                ["range", path, out] => {
                    read_input(Path::new(path)).and_then(|input| run_range(&input, Path::new(out)))
                }
                ["nuclide", path, range, name, out] => read_input(Path::new(path))
                    .and_then(|input| run_nuclide(&input, Path::new(range), name, Path::new(out))),
                ["assemble", path, range, directory] => read_input(Path::new(path))
                    .and_then(|input| run_assemble(&input, Path::new(range), Path::new(directory))),
                ["check", out, ref files @ ..] if !files.is_empty() => {
                    run_check(files, Path::new(out))
                }
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
    in_parallel(nuclides.len(), |i| read_nuclide_tape(&nuclides[i]))
}

/// 核種のテープを読み、入力の宣言した修正を当てる（ADR-0272）。
fn read_nuclide_tape(nuclide: &NuclideInput) -> anyhow::Result<Tape> {
    let mut tape = read_tape(&nuclide.endf)?;
    correct_tape(&nuclide.name, &mut tape, &nuclide.corrections)
        .with_context(|| format!("{}", nuclide.endf.display()))?;
    Ok(tape)
}

/// 宣言の順に修正を当て、当てるたびに warn を出す。
fn correct_tape(name: &str, tape: &mut Tape, corrections: &[Correction]) -> anyhow::Result<()> {
    for c in corrections {
        tape.correct(c)?;
        xlc_log::warn!(
            "{name}: corrected line {} (MF{} MT{}) field {} from {} to {}, as the input declares",
            c.line,
            c.mf,
            c.mt,
            c.field,
            c.from,
            c.to
        );
    }
    Ok(())
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

fn library_of(
    settings: Settings,
    nuclides: Vec<NuclideData>,
    activation: Vec<ActivationData>,
) -> Library {
    Library {
        inverse_velocity: inverse_velocity(&settings.groups, &settings.weight),
        groups: settings.groups,
        sigma0: settings.sigma0,
        pointwise_range: settings.pointwise_range.map(|(low, high)| [low, high]),
        nuclides,
        activation,
    }
}

/// 輸送する核種（放射化の核種を除く。ポイントワイズの範囲を決める核種）
fn transported(nuclides: &[NuclideInput]) -> impl Iterator<Item = &NuclideInput> {
    nuclides.iter().filter(|n| !n.activation)
}

/// 核種 1 つの処理の結果（輸送する核種か放射化の核種）
enum Processed {
    Transported(NuclideData),
    Activation(ActivationData),
}

/// 処理の結果を、輸送する核種と放射化の核種に分けて 1 つのライブラリにする
fn library_from(settings: Settings, results: Vec<Processed>) -> Library {
    let (mut nuclides, mut activation) = (Vec::new(), Vec::new());
    for result in results {
        match result {
            Processed::Transported(data) => nuclides.push(data),
            Processed::Activation(data) => activation.push(data),
        }
    }
    library_of(settings, nuclides, activation)
}

/// 輸送する核種のテープだけを読む（ポイントワイズの範囲を求める）
fn transported_tapes(nuclides: &[NuclideInput]) -> anyhow::Result<(Vec<NuclideInput>, Vec<Tape>)> {
    let transported: Vec<NuclideInput> = transported(nuclides).cloned().collect();
    let tapes = read_tapes(&transported)?;
    Ok((transported, tapes))
}

/// ファイルを書く。DVC はステージの前に出力を消すとき、空になったディレクトリも消すので、ディレクトリを作る。
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
    let (kept, kept_tapes): (Vec<NuclideInput>, Vec<Tape>) = input
        .nuclides
        .iter()
        .zip(&tapes)
        .filter(|(n, _)| !n.activation)
        .map(|(n, t)| (n.clone(), t.clone()))
        .unzip();
    let pointwise_range = pointwise_range(&kept, &kept_tapes, &GroupStructure::wims69())?;
    let settings = settings_of(input, pointwise_range);
    let results = in_parallel(input.nuclides.len(), |i| {
        process_file(&input.nuclides[i], &tapes[i], &settings)
    })?;
    write_library(&library_from(settings, results), &input.output)
}

/// すべての核種の分離共鳴の上端から、ポイントワイズの範囲を求めて書く（`nuclide` のステージの入力になる）。
fn run_range(input: &Input, out: &Path) -> anyhow::Result<()> {
    let (nuclides, tapes) = transported_tapes(&input.nuclides)?;
    let range = pointwise_range(&nuclides, &tapes, &GroupStructure::wims69())?;
    write_file(out, range::to_toml(range).as_bytes())
}

fn read_range(path: &Path) -> anyhow::Result<Option<(f64, f64)>> {
    let context = || format!("{}", path.display());
    let text = std::fs::read_to_string(path).with_context(context)?;
    range::from_toml(&text).with_context(context)
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
    let data = process_file(nuclide, &read_nuclide_tape(nuclide)?, &settings)?;
    write_library(&library_from(settings, vec![data]), out)
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

/// ENDF のファイルごとに読み込みと形式の検査を行い（ADR-0266）、止まった理由をファイルごとに報告に書く。
fn run_check(paths: &[&str], out: &Path) -> anyhow::Result<()> {
    let checked = in_parallel(paths.len(), |i| check_file(paths[i]))?;
    let stopped = checked.iter().filter(|c| c.error.is_some()).count();
    xlc_log::info!(
        "{stopped} of {} evaluations stop the reading checks",
        checked.len()
    );
    write_file(out, report_toml(&checked)?.as_bytes())
}

/// ファイル `path` の検査の結果。ファイルが読めなければエラー（中身がテープとして読めないのは、検査で止まった理由）。
fn check_file(path: &str) -> anyhow::Result<Checked> {
    let bytes = std::fs::read(path).with_context(|| path.to_string())?;
    Ok(Checked {
        endf: path.to_string(),
        error: check_contents(&bytes).err(),
    })
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

/// 核種ごとの結果が、入力のその核種の結果であること（1 核種だけで、名前と温度が入力と同じ）を確かめる。
fn check_part(nuclide: &NuclideInput, part: &Library) -> anyhow::Result<()> {
    let (name, temperatures): (&str, Vec<f64>) =
        match (&part.nuclides[..], &part.activation[..], nuclide.activation) {
            ([data], [], false) => (
                &data.name,
                data.temperatures.iter().map(|t| t.temperature).collect(),
            ),
            ([], [data], true) => (
                &data.name,
                data.temperatures.iter().map(|t| t.temperature).collect(),
            ),
            (nuclides, activation, _) => {
                return Err(anyhow!(
                    "{}: expected one {} nuclide, found {} transported and {} activation nuclides",
                    nuclide.name,
                    if nuclide.activation {
                        "activation"
                    } else {
                        "transported"
                    },
                    nuclides.len(),
                    activation.len()
                ));
            }
        };
    if name != nuclide.name {
        return Err(anyhow!("{}: the result is for {name}", nuclide.name));
    }
    if temperatures != nuclide.temperatures {
        return Err(anyhow!(
            "{}: the result has temperatures {temperatures:?}, the input has {:?}",
            nuclide.name,
            nuclide.temperatures
        ));
    }
    if let Some(data) = part.nuclides.first()
        && data
            .temperatures
            .iter()
            .any(|t| t.pointwise.is_some() != nuclide.pointwise)
    {
        return Err(anyhow!(
            "{}: the result {} pointwise data, the input has `pointwise = {}`",
            nuclide.name,
            if nuclide.pointwise { "lacks" } else { "has" },
            nuclide.pointwise
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
    let velocity = inverse_velocity(&part.groups, &Weight::lwr());
    if part.inverse_velocity != velocity {
        return Err(anyhow!(
            "the inverse velocities are {:?}, the weight function gives {velocity:?}",
            part.inverse_velocity
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

/// 輸送する核種の結果 `data` を、入力 `nuclide` の `pointwise` に合わせる（false ならポイントワイズの断面積を外す）。
/// 処理の中で取り出さないようにしないのは、取り出しの時間が測ったばらつきの中だから（ADR-0269 の案 P1）
fn keep_pointwise(nuclide: &NuclideInput, mut data: NuclideData) -> NuclideData {
    if !nuclide.pointwise {
        for t in &mut data.temperatures {
            t.pointwise = None;
        }
    }
    data
}

/// 核種を処理する（処理の時間をログに出す）。
#[tracing::instrument(name = "nuclide", skip_all, fields(name = %nuclide.name))]
fn process_file(
    nuclide: &NuclideInput,
    tape: &Tape,
    settings: &Settings,
) -> anyhow::Result<Processed> {
    let name = &nuclide.name;
    if nuclide.activation {
        let data = process_activation(tape, name, &nuclide.temperatures, settings)
            .with_context(|| name.clone())?;
        return Ok(Processed::Activation(data));
    }
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
    Ok(Processed::Transported(keep_pointwise(nuclide, data)))
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
#[path = "main/tests.rs"]
mod tests;
