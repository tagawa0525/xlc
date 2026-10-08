//! リポジトリの文書の整合性の検査（#0013）。
//!
//! git が追跡しているファイルを一度だけ読み、xlc_docs::Docs の検査に違反がないことを確かめる。閉じた issue の
//! `resolution: コミット <hash>` は、git で HEAD の祖先のコミットを指すことを確かめる。違反はすべて並べて示す。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use xlc_docs::Docs;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

/// 中身を読むファイルの拡張子（文書と、参照を書くソースとスクリプト）。
const TEXT_EXTENSIONS: [&str; 5] = ["md", "rs", "py", "toml", "sh"];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join("..")
}

/// リポジトリの根で git を走らせる。git がなければエラーにする。
fn git(root: &Path, args: &[&str]) -> Result<Output> {
    Command::new("git")
        .arg("-C")
        .arg(root)
        .args(args)
        .output()
        .map_err(|e| format!("cannot run git {args:?} (is git installed?): {e}").into())
}

/// git が追跡しているファイル（`/` 区切りのパス）。
fn tracked_files(root: &Path) -> Result<Vec<String>> {
    let output = git(root, &["ls-files", "-z"])?;
    if !output.status.success() {
        return Err(format!(
            "git ls-files failed: {}",
            String::from_utf8_lossy(&output.stderr)
        )
        .into());
    }
    Ok(String::from_utf8(output.stdout)?
        .split('\0')
        .filter(|p| !p.is_empty())
        .map(str::to_string)
        .collect())
}

/// `hash` が HEAD の祖先のコミットを指すか。違えば理由を返す。
fn commit_problem(root: &Path, hash: &str) -> Result<Option<String>> {
    let commit = format!("{hash}^{{commit}}");
    if !git(root, &["rev-parse", "--quiet", "--verify", &commit])?
        .status
        .success()
    {
        return Ok(Some(format!("commit {hash} does not exist")));
    }
    let ancestor = git(root, &["merge-base", "--is-ancestor", hash, "HEAD"])?;
    Ok(match ancestor.status.code() {
        Some(0) => None,
        Some(1) => Some(format!(
            "commit {hash} is not an ancestor of HEAD (rebased? update the hash)"
        )),
        _ => {
            return Err(format!(
                "git merge-base failed: {}",
                String::from_utf8_lossy(&ancestor.stderr)
            )
            .into());
        }
    })
}

#[test]
fn docs_are_consistent() -> Result<()> {
    let root = repo_root();
    let mut texts = BTreeMap::new();
    let mut others = Vec::new();
    for path in tracked_files(&root)? {
        let extension = path.rsplit_once('.').map_or("", |(_, e)| e);
        // 追跡していても作業ツリーで消したファイルは、読まずに数えない
        let file = root.join(&path);
        if !file.exists() {
            continue;
        }
        if TEXT_EXTENSIONS.contains(&extension) {
            texts.insert(path, fs::read_to_string(file)?);
        } else {
            others.push(path);
        }
    }
    assert!(
        texts.contains_key("CLAUDE.md") && texts.contains_key("docs/adr/README.md"),
        "the repository root was not found at {}",
        root.display()
    );

    let docs = Docs::new(texts, others);
    let mut violations: Vec<String> = docs.check().iter().map(ToString::to_string).collect();
    for (path, hash) in docs.resolution_commits() {
        if let Some(problem) = commit_problem(&root, &hash)? {
            violations.push(format!("{path}: `resolution` names {problem}"));
        }
    }
    assert!(
        violations.is_empty(),
        "{} documentation violation(s):\n{}",
        violations.len(),
        violations.join("\n")
    );
    Ok(())
}
