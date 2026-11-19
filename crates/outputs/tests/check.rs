//! 形式を変える修正の検査（scripts/check-readable-outputs.py）は、本番の読み手（read-outputs）で dvc.lock の指す上流の出力を
//! 読み、読めない出力とそれを書いたステージを出して止まる。open の issue に結び付けた出力だけを通す（ADR-0258）。
#![allow(
    clippy::unwrap_used,
    clippy::disallowed_macros,
    clippy::panic,
    reason = "a test: a failed run should stop it with the message"
)]

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use rmpv::Value;
use xlc_format::medium::{Medium, MediumState};

const MEDIUM_MD5: &str = "11111111111111111111111111111111";

const ROOT_YAML: &str = "\
stages:
  build:
    cmd: scripts/build-binaries.sh
    deps:
      - crates
";

const LOCK: &str = "\
schema: '2.0'
stages:
  medium:
    cmd: ../../bin/lattice medium.toml
    outs:
    - path: medium.msgpack
      hash: md5
      md5: 11111111111111111111111111111111
  axial:
    cmd: ../../bin/lattice axial.toml
    deps:
    - path: medium.msgpack
      hash: md5
      md5: 11111111111111111111111111111111
    - path: ../../binary-forms/lattice.md5
      hash: md5
      md5: 22222222222222222222222222222222
    outs:
    - path: axial.out
      hash: md5
      md5: 33333333333333333333333333333333
";

const OPEN_ISSUE: &str = "\
---
type: Issue
status: open
---

# 媒質を作り直す
";

fn medium() -> Vec<u8> {
    Medium {
        bounds: vec![2.0e7, 0.625, 1.0e-5],
        axes: vec![],
        states: vec![MediumState {
            total: vec![0.5, 1.0],
            scatter: vec![vec![vec![0.4, 0.05], vec![0.0, 0.9]]],
            fission: vec![0.0, 0.02],
            nu_fission: vec![0.0, 0.05],
            chi: vec![1.0, 0.0],
        }],
    }
    .to_msgpack()
    .unwrap()
}

/// 形式の名前を変えた読み手から見た古い出力（名前だけが違う）
fn renamed_medium() -> Vec<u8> {
    let mut value = rmpv::decode::read_value(&mut &medium()[..]).unwrap();
    let Value::Map(entries) = &mut value else {
        panic!("not a map")
    };
    for (k, v) in entries.iter_mut() {
        if k.as_str() == Some("format") {
            *v = Value::from("xlc-medium-0");
        }
    }
    let mut out = Vec::new();
    rmpv::encode::write_value(&mut out, &value).unwrap();
    out
}

struct Repo {
    _tmp: tempfile_dir::Dir,
    root: PathBuf,
    cache: PathBuf,
}

/// 依存を足さずに、テストごとの一時ディレクトリを作って消す
mod tempfile_dir {
    use std::path::PathBuf;

    pub struct Dir(pub PathBuf);

    impl Dir {
        pub fn new(name: &str) -> Self {
            let path =
                std::env::temp_dir().join(format!("xlc-outputs-{name}-{}", std::process::id()));
            let _ = std::fs::remove_dir_all(&path);
            std::fs::create_dir_all(&path).unwrap();
            Self(path)
        }
    }

    impl Drop for Dir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
}

fn git(root: &Path, args: &[&str]) {
    // 利用者の git の hook（core.hooksPath）を一時のリポジトリに掛けない
    let status = Command::new("git")
        .arg("-C")
        .arg(root)
        .args(["-c", "core.hooksPath=/dev/null"])
        .args(args)
        .status()
        .unwrap();
    assert!(status.success(), "git {args:?}");
}

fn write(path: &Path, bytes: &[u8]) {
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    std::fs::write(path, bytes).unwrap();
}

impl Repo {
    fn new(name: &str, medium: &[u8]) -> Self {
        let tmp = tempfile_dir::Dir::new(name);
        let root = tmp.0.join("repo");
        let cache = tmp.0.join("cache");
        write(&root.join("dvc.yaml"), ROOT_YAML.as_bytes());
        write(&root.join("examples/a/dvc.lock"), LOCK.as_bytes());
        write(
            &root.join("docs/issues/0001-rebuild.md"),
            OPEN_ISSUE.as_bytes(),
        );
        git(&root, &["init", "-q", "-b", "main"]);
        git(&root, &["add", "."]);
        git(
            &root,
            &[
                "-c",
                "user.name=t",
                "-c",
                "user.email=t@t",
                "commit",
                "-qm",
                "base",
            ],
        );
        write(
            &cache.join(format!(
                "files/md5/{}/{}",
                &MEDIUM_MD5[..2],
                &MEDIUM_MD5[2..]
            )),
            medium,
        );
        Self {
            _tmp: tmp,
            root,
            cache,
        }
    }

    fn check(&self) -> Output {
        let script =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../scripts/check-readable-outputs.py");
        Command::new("python3")
            .arg(script)
            .arg("--all")
            .arg("--cache")
            .arg(&self.cache)
            .arg("--reader")
            .arg(env!("CARGO_BIN_EXE_read-outputs"))
            .current_dir(&self.root)
            .output()
            .unwrap()
    }
}

fn shown(output: &Output) -> String {
    format!(
        "{}{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    )
}

#[test]
fn passes_when_the_readers_read_every_locked_output() {
    let repo = Repo::new("readable", &medium());
    let output = repo.check();
    assert_eq!(output.status.code(), Some(0), "{}", shown(&output));
}

#[test]
fn stops_on_an_output_the_readers_cannot_read_and_names_its_stage() {
    let repo = Repo::new("renamed", &renamed_medium());
    let output = repo.check();
    assert_eq!(output.status.code(), Some(1), "{}", shown(&output));
    let text = shown(&output);
    assert!(
        text.contains(&format!(
            "UNREADABLE: examples/a/medium.msgpack (md5 {MEDIUM_MD5}, written by examples/a/dvc.yaml:medium)"
        )),
        "{text}"
    );
    assert!(text.contains("xlc-medium-0"), "{text}");
}

#[test]
fn passes_an_unreadable_output_tied_to_an_open_issue() {
    let repo = Repo::new("allowed", &renamed_medium());
    write(
        &repo.root.join("scripts/unreadable-outputs.toml"),
        format!(
            "[[output]]\nmd5 = \"{MEDIUM_MD5}\"\nissue = \"#0001\"\npath = \"examples/a/medium.msgpack\"\n"
        )
        .as_bytes(),
    );
    let output = repo.check();
    assert_eq!(output.status.code(), Some(0), "{}", shown(&output));
    assert!(
        shown(&output).contains("ALLOWED: examples/a/medium.msgpack"),
        "{}",
        shown(&output)
    );
}
