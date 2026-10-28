//! lattice のスレッドの数: `RAYON_NUM_THREADS` が無ければ、このプロセスに許された論理 CPU（SMT のハードウェアスレッドを
//! 含む）の数で解き、あればその数で解く。どちらも始めにログ（標準エラー）に出す（ADR-0214）。
#![allow(
    clippy::unwrap_used,
    clippy::disallowed_macros,
    clippy::disallowed_methods,
    clippy::panic,
    reason = "a test: a failed run should stop it with the message"
)]

use std::process::Command;

/// 入力の無いパスで lattice を実行し（入力を読む前に止まる）、標準エラーを返す
fn stderr_of(rayon_num_threads: Option<&str>) -> String {
    let mut command = Command::new(env!("CARGO_BIN_EXE_lattice"));
    command
        .arg("/nonexistent/xlc-threads-test.toml")
        .env_remove("RUST_LOG")
        .env_remove("RAYON_NUM_THREADS");
    if let Some(count) = rayon_num_threads {
        command.env("RAYON_NUM_THREADS", count);
    }
    let output = command.output().unwrap();
    assert!(!output.status.success(), "the missing input should fail");
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
fn uses_every_allowed_logical_cpu_without_rayon_num_threads() {
    let allowed = std::thread::available_parallelism().unwrap().get();
    let stderr = stderr_of(None);
    assert_eq!(field(&stderr, "threads").parse::<usize>().unwrap(), allowed);
    assert_eq!(field(&stderr, "from"), "\"default\"");
}

#[test]
fn uses_rayon_num_threads_when_given() {
    let stderr = stderr_of(Some("3"));
    assert_eq!(field(&stderr, "threads"), "3");
    assert_eq!(field(&stderr, "from"), "\"RAYON_NUM_THREADS\"");
}

#[test]
fn rayon_num_threads_that_is_not_a_positive_count_is_an_error() {
    for value in ["abc", "0", ""] {
        let stderr = stderr_of(Some(value));
        assert!(
            stderr.contains("error: RAYON_NUM_THREADS"),
            "RAYON_NUM_THREADS={value:?} should stop lattice:\n{stderr}"
        );
        assert!(
            !stderr.contains("threads="),
            "no pool for {value:?}:\n{stderr}"
        );
    }
}
