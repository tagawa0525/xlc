//! core のスレッドの数: ケースの `threads` があればその数、無ければ `RAYON_NUM_THREADS`、どちらも無ければ物理コアの数で解き、始めにログ
//! （標準エラー）に出す。両方を与えるとエラー（ADR-0234）。
#![allow(
    clippy::unwrap_used,
    clippy::disallowed_macros,
    clippy::disallowed_methods,
    clippy::panic,
    reason = "a test: a failed run should stop it with the message"
)]

use std::path::PathBuf;
use std::process::Command;

/// サイクルの無いケース（`threads` の行を足せる）を一時ディレクトリに書き、core を実行して（サイクルを読むところで止まる）、標準エラーを返す
fn stderr_of(name: &str, threads: Option<&str>, rayon_num_threads: Option<&str>) -> String {
    let dir: PathBuf =
        std::env::temp_dir().join(format!("xlc-core-threads-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut case = String::from(
        "cycle = \"missing-cycle.toml\"\naxial-mesh = [0.0, 10.0]\nradial-nodes = 1\n\
         boundary = { x = [\"reflective\", \"reflective\"], y = [\"reflective\", \"reflective\"], \
         z = [\"reflective\", \"reflective\"] }\n",
    );
    if let Some(count) = threads {
        case.push_str(&format!("threads = {count}\n"));
    }
    let path = dir.join("case.toml");
    std::fs::write(&path, case).unwrap();
    let mut command = Command::new(env!("CARGO_BIN_EXE_core"));
    command
        .arg(&path)
        .env_remove("RUST_LOG")
        .env_remove("RAYON_NUM_THREADS");
    if let Some(count) = rayon_num_threads {
        command.env("RAYON_NUM_THREADS", count);
    }
    let output = command.output().unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
    assert!(!output.status.success(), "the missing cycle should fail");
    String::from_utf8(output.stderr).unwrap()
}

/// ログの `<key>=` の値
fn field<'a>(stderr: &'a str, key: &str) -> &'a str {
    stderr
        .split_whitespace()
        .find_map(|word| word.strip_prefix(key)?.strip_prefix('='))
        .unwrap_or_else(|| panic!("no `{key}=` in the log:\n{stderr}"))
}

#[test]
fn uses_the_physical_cores_without_an_input_or_rayon_num_threads() {
    // 既定は物理コアの数を、このプロセスに許された論理 CPU の数で抑えたもの（ADR-0234、ADR-0067）
    let allowed = std::thread::available_parallelism().unwrap().get();
    let stderr = stderr_of("default", None, None);
    let threads: usize = field(&stderr, "threads").parse().unwrap();
    assert_eq!(threads, num_cpus::get_physical().min(allowed));
    assert_eq!(field(&stderr, "from"), "\"default\"");
    assert!(stderr.contains("missing-cycle.toml"), "{stderr}");
}

#[test]
fn uses_the_case_threads_when_given() {
    let stderr = stderr_of("input", Some("3"), None);
    assert_eq!(field(&stderr, "threads"), "3");
    assert_eq!(field(&stderr, "from"), "\"input\"");
}

#[test]
fn uses_rayon_num_threads_when_given() {
    let stderr = stderr_of("rayon", None, Some("2"));
    assert_eq!(field(&stderr, "threads"), "2");
    assert_eq!(field(&stderr, "from"), "\"RAYON_NUM_THREADS\"");
}

#[test]
fn case_threads_and_rayon_num_threads_together_are_an_error() {
    let stderr = stderr_of("both", Some("3"), Some("2"));
    assert!(
        stderr.contains("error:") && stderr.contains("RAYON_NUM_THREADS"),
        "{stderr}"
    );
    assert!(!stderr.contains("threads="), "no pool:\n{stderr}");
}

#[test]
fn rayon_num_threads_that_is_not_a_positive_count_is_an_error() {
    for value in ["abc", "0", ""] {
        let stderr = stderr_of("bad", None, Some(value));
        assert!(
            stderr.contains("error: RAYON_NUM_THREADS"),
            "RAYON_NUM_THREADS={value:?} should stop core:\n{stderr}"
        );
        assert!(
            !stderr.contains("threads="),
            "no pool for {value:?}:\n{stderr}"
        );
    }
}
