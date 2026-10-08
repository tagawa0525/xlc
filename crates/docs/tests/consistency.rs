//! リポジトリの文書の整合性の検査（#0013）。
//!
//! リポジトリのファイルを一度だけ読み、xlc_docs::Docs の検査に違反がないことを確かめる。違反はすべて並べて示す。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

use std::collections::BTreeMap;
use std::fs;
use std::io;
use std::path::Path;

use xlc_docs::Docs;

/// 中身を読むファイルの拡張子（文書と、参照を書くソースとスクリプト）。
const TEXT_EXTENSIONS: [&str; 5] = ["md", "rs", "py", "toml", "sh"];
/// 読まない根のディレクトリ。references/ は git で追跡しない（docs/environment.md）ので、文書からリンクしない。
const SKIPPED_TOP_DIRS: [&str; 2] = ["target", "references"];

/// `dir` の下のファイルを集める。隠しファイルと隠しディレクトリは除く。
fn collect(
    root: &Path,
    dir: &Path,
    texts: &mut BTreeMap<String, String>,
    others: &mut Vec<String>,
) -> io::Result<()> {
    for entry in fs::read_dir(dir)? {
        let path = entry?.path();
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if name.starts_with('.') || (dir == root && SKIPPED_TOP_DIRS.contains(&name)) {
            continue;
        }
        if path.is_dir() {
            collect(root, &path, texts, others)?;
            continue;
        }
        let relative = path
            .strip_prefix(root)
            .unwrap_or(&path)
            .components()
            .map(|c| c.as_os_str().to_string_lossy())
            .collect::<Vec<_>>()
            .join("/");
        let extension = path.extension().and_then(|e| e.to_str()).unwrap_or("");
        if TEXT_EXTENSIONS.contains(&extension) {
            texts.insert(relative, fs::read_to_string(&path)?);
        } else {
            others.push(relative);
        }
    }
    Ok(())
}

#[test]
fn docs_are_consistent() -> io::Result<()> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join("..");
    let mut texts = BTreeMap::new();
    let mut others = Vec::new();
    collect(&root, &root, &mut texts, &mut others)?;
    assert!(
        texts.contains_key("CLAUDE.md") && texts.contains_key("docs/adr/README.md"),
        "the repository root was not found at {}",
        root.display()
    );

    let violations = Docs::new(texts, others).check();
    assert!(
        violations.is_empty(),
        "{} documentation violation(s):\n{}",
        violations.len(),
        violations
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join("\n")
    );
    Ok(())
}
