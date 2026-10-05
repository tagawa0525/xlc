//! ENDF ファイルから多群ライブラリを作るコマンド。
//!
//! 使い方: xs [--temperature <K>] <出力ファイル> <名前>=<ENDF ファイル>...

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use xs::endf::Tape;
use xs::library::{Library, Nuclide};
use xs::multigroup::{GroupStructure, Weight};
use xs::pointwise::DEFAULT_TOLERANCE;
use xs::process::{Settings, process_nuclide};

/// コマンドライン引数を解釈した結果。
#[derive(Debug, Clone, PartialEq)]
struct Command {
    output: PathBuf,
    /// (核種名, ENDF ファイル)
    nuclides: Vec<(String, PathBuf)>,
    /// 自由ガス模型の温度 (K)
    temperature: f64,
}

const USAGE: &str = "usage: xs [--temperature <K>] <output> <NAME>=<ENDF>...";

/// 既定の温度 (K)。kT = 0.0253 eV（重み関数の熱 Maxwell 分布）に合わせる。
const DEFAULT_TEMPERATURE: f64 = 293.6;

fn parse_args(args: &[String]) -> Result<Command, String> {
    let mut temperature = DEFAULT_TEMPERATURE;
    let mut rest = args;
    if let [flag, value, tail @ ..] = rest
        && flag == "--temperature"
    {
        temperature = value
            .parse::<f64>()
            .ok()
            .filter(|t| *t > 0.0 && t.is_finite())
            .ok_or_else(|| format!("invalid temperature `{value}`"))?;
        rest = tail;
    }
    let [output, specs @ ..] = rest else {
        return Err(USAGE.into());
    };
    if specs.is_empty() {
        return Err(USAGE.into());
    }
    let mut nuclides: Vec<(String, PathBuf)> = Vec::new();
    for spec in specs {
        let (name, path) = spec
            .split_once('=')
            .filter(|(name, path)| {
                !name.is_empty() && !path.is_empty() && !name.contains(char::is_whitespace)
            })
            .ok_or_else(|| {
                format!("expected <NAME>=<ENDF> without spaces in NAME, found `{spec}`")
            })?;
        if nuclides.iter().any(|(n, _)| n == name) {
            return Err(format!("nuclide `{name}` is given twice"));
        }
        nuclides.push((name.to_string(), PathBuf::from(path)));
    }
    Ok(Command {
        output: PathBuf::from(output),
        nuclides,
        temperature,
    })
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match parse_args(&args).and_then(|command| run(&command)) {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("error: {message}");
            ExitCode::FAILURE
        }
    }
}

fn run(command: &Command) -> Result<(), String> {
    let settings = Settings {
        groups: GroupStructure::wims69(),
        weight: Weight::lwr(),
        temperature: command.temperature,
        thermal_cutoff: THERMAL_CUTOFF,
        tolerance: DEFAULT_TOLERANCE,
    };
    // 核種ごとの処理は独立なので並列に行う
    let results: Vec<Result<Nuclide, String>> = std::thread::scope(|scope| {
        let handles: Vec<_> = command
            .nuclides
            .iter()
            .map(|(name, path)| scope.spawn(|| process_file(name, path, &settings)))
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
        nuclides,
    };
    std::fs::write(&command.output, library.to_text())
        .map_err(|e| format!("{}: {e}", command.output.display()))
}

fn process_file(name: &str, path: &Path, settings: &Settings) -> Result<Nuclide, String> {
    let start = Instant::now();
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let tape = Tape::parse(&text).map_err(|e| format!("{}: {e}", path.display()))?;
    let nuclide = process_nuclide(&tape, name, settings).map_err(|e| format!("{name}: {e}"))?;
    eprintln!("{name}: done in {:.1} s", start.elapsed().as_secs_f64());
    Ok(nuclide)
}

/// この入射エネルギー (eV) 未満の弾性散乱を自由ガス模型で扱う（WIMS 69群の熱群の上限）
const THERMAL_CUTOFF: f64 = 4.0;

#[cfg(test)]
mod tests {
    use super::*;

    fn args(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn parses_nuclides_and_default_temperature() {
        let c = parse_args(&args(&["lib.txt", "H1=h.endf", "U235=u.endf"])).unwrap();
        assert_eq!(c.output, PathBuf::from("lib.txt"));
        assert_eq!(
            c.nuclides,
            vec![
                ("H1".to_string(), PathBuf::from("h.endf")),
                ("U235".to_string(), PathBuf::from("u.endf"))
            ]
        );
        assert_eq!(c.temperature, DEFAULT_TEMPERATURE);
    }

    #[test]
    fn parses_temperature() {
        let c = parse_args(&args(&["--temperature", "565", "lib.txt", "H1=h.endf"])).unwrap();
        assert_eq!(c.temperature, 565.0);
    }

    #[test]
    fn rejects_bad_arguments() {
        assert!(parse_args(&args(&[])).is_err());
        assert!(parse_args(&args(&["lib.txt"])).is_err());
        assert!(parse_args(&args(&["lib.txt", "h.endf"])).is_err());
        assert!(parse_args(&args(&["lib.txt", "=h.endf"])).is_err());
        assert!(parse_args(&args(&["--temperature", "-1", "lib.txt", "H1=h"])).is_err());
        assert!(parse_args(&args(&["lib.txt", "H1=a", "H1=b"])).is_err());
        // ライブラリでは名前を空白で区切るので、空白を含む名前は書けても読めない
        assert!(parse_args(&args(&["lib.txt", "U 235=u.endf"])).is_err());
    }
}
