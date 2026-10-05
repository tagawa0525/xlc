//! ENDF ファイルから多群ライブラリを作るコマンド。
//!
//! 使い方: xs [--temperature <K>] <出力ファイル> <名前>=<ENDF ファイル>...

use std::path::PathBuf;
use std::process::ExitCode;

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
    todo!()
}

fn main() -> ExitCode {
    todo!()
}

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
    }
}
