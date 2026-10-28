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

/// 入力の無いパスで lattice を実行し（入力を読む前に止まる）、ログの `threads=` の数を返す
fn logged_threads(rayon_num_threads: Option<&str>) -> usize {
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
    let stderr = String::from_utf8(output.stderr).unwrap();
    let field = stderr
        .split_whitespace()
        .find_map(|word| word.strip_prefix("threads="))
        .unwrap_or_else(|| panic!("no `threads=` in the log:\n{stderr}"));
    field.parse().unwrap()
}

#[test]
fn uses_every_allowed_logical_cpu_without_rayon_num_threads() {
    let allowed = std::thread::available_parallelism().unwrap().get();
    assert_eq!(logged_threads(None), allowed);
}

#[test]
fn uses_rayon_num_threads_when_given() {
    assert_eq!(logged_threads(Some("3")), 3);
}
