//! リポジトリの文書の整合性の検査（#0013）。
//!
//! docs/adr/README.md と docs/issues/README.md が決めている約束を確かめる。違反はすべて集めてから、まとめて
//! 失敗として示す。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};

use xlc_docs::{
    FrontMatter, adr_refs, issue_refs, numbered_file, parse_front_matter, relative_links,
};

const ADR_KEYS: [&str; 5] = ["status", "date", "supersedes", "superseded-by", "issues"];
const ADR_STATUSES: [&str; 5] = [
    "proposed",
    "accepted",
    "rejected",
    "superseded",
    "deprecated",
];
const ISSUE_KEYS: [&str; 7] = [
    "status",
    "labels",
    "opened",
    "closed",
    "milestone",
    "blocked-by",
    "resolution",
];
const ISSUE_STATUSES: [&str; 2] = ["open", "closed"];
/// issue の本文で、関係する ADR を挙げる行の書き出し（docs/issues/README.md）。
const RELATED_ADRS: &str = "関係する ADR:";

type Result<T> = std::result::Result<T, Box<dyn Error>>;

/// 番号つきの ADR か issue の1件。
struct Doc {
    /// リポジトリの根からのパス（表示用）。
    path: String,
    number: u32,
    front: FrontMatter,
    text: String,
}

impl Doc {
    fn status(&self) -> &str {
        self.items("status").first().copied().unwrap_or("")
    }

    fn items(&self, key: &str) -> Vec<&str> {
        self.front.get(key).map(|v| v.items()).unwrap_or_default()
    }
}

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join("..")
}

fn relative(root: &Path, path: &Path) -> String {
    path.strip_prefix(root)
        .unwrap_or(path)
        .display()
        .to_string()
}

/// `dir` の下の、拡張子が `extensions` のどれかのファイル（隠しディレクトリと target は除く）。
fn files_under(dir: &Path, extensions: &[&str], out: &mut Vec<PathBuf>) -> Result<()> {
    for entry in fs::read_dir(dir)? {
        let path = entry?.path();
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if name.starts_with('.') || name == "target" {
            continue;
        }
        if path.is_dir() {
            files_under(&path, extensions, out)?;
        } else if path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| extensions.contains(&e))
        {
            out.push(path);
        }
    }
    Ok(())
}

/// `dir` の `NNNN-name.md`（README.md を除く）を読む。
fn load_numbered(root: &Path, dir: &str, violations: &mut Vec<String>) -> Result<Vec<Doc>> {
    let mut docs = Vec::new();
    let mut paths: Vec<PathBuf> = fs::read_dir(root.join(dir))?
        .map(|e| e.map(|e| e.path()))
        .collect::<std::result::Result<_, _>>()?;
    paths.sort();
    for path in paths {
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if name == "README.md" {
            continue;
        }
        let shown = relative(root, &path);
        let Some(number) = numbered_file(name) else {
            violations.push(format!("{shown}: file name is not NNNN-name.md"));
            continue;
        };
        let text = fs::read_to_string(&path)?;
        match parse_front_matter(&text) {
            Ok(front) => docs.push(Doc {
                path: shown,
                number,
                front,
                text,
            }),
            Err(e) => violations.push(format!("{shown}: {e}")),
        }
    }
    Ok(docs)
}

/// 番号の重複（1）。
fn check_unique(docs: &[Doc], violations: &mut Vec<String>) {
    let mut by_number: BTreeMap<u32, Vec<&str>> = BTreeMap::new();
    for doc in docs {
        by_number.entry(doc.number).or_default().push(&doc.path);
    }
    for (number, paths) in by_number {
        if paths.len() > 1 {
            violations.push(format!(
                "number {number:04} is used by {}",
                paths.join(", ")
            ));
        }
    }
}

/// 前付けのキーと `status` の値（2）。
fn check_keys(docs: &[Doc], keys: &[&str], statuses: &[&str], violations: &mut Vec<String>) {
    for doc in docs {
        for key in keys {
            if doc.front.get(key).is_none() {
                violations.push(format!("{}: front matter has no `{key}`", doc.path));
            }
        }
        for key in doc.front.keys() {
            if !keys.contains(&key) {
                violations.push(format!("{}: unexpected front matter key `{key}`", doc.path));
            }
        }
        if !statuses.contains(&doc.status()) {
            violations.push(format!(
                "{}: status {:?} is not one of {statuses:?}",
                doc.path,
                doc.status()
            ));
        }
    }
}

/// `ADR-NNNN` か `#NNNN` の一つだけの文字列なら、その番号。
fn single_ref(item: &str, refs: fn(&str) -> Vec<u32>, prefix: &str) -> Option<u32> {
    match refs(item).as_slice() {
        [n] if item == format!("{prefix}{n:04}") => Some(*n),
        _ => None,
    }
}

/// 前付けの `key` の要素を番号にする。形が違うものと、存在しない番号は違反にする。
fn front_refs(
    doc: &Doc,
    key: &str,
    refs: fn(&str) -> Vec<u32>,
    prefix: &str,
    existing: &BTreeSet<u32>,
    violations: &mut Vec<String>,
) -> BTreeSet<u32> {
    let mut numbers = BTreeSet::new();
    for item in doc.items(key) {
        match single_ref(item, refs, prefix) {
            Some(n) if existing.contains(&n) => {
                numbers.insert(n);
            }
            Some(_) => violations.push(format!("{}: `{key}` names missing {item}", doc.path)),
            None => violations.push(format!(
                "{}: `{key}` item {item:?} is not {prefix}NNNN",
                doc.path
            )),
        }
    }
    numbers
}

/// 置き換えの状態とリンク（3）。
fn check_supersession(adrs: &[Doc], violations: &mut Vec<String>) {
    let existing: BTreeSet<u32> = adrs.iter().map(|d| d.number).collect();
    let mut supersedes: BTreeMap<u32, BTreeSet<u32>> = BTreeMap::new();
    let mut superseded_by: BTreeMap<u32, BTreeSet<u32>> = BTreeMap::new();
    for adr in adrs {
        let by = front_refs(
            adr,
            "superseded-by",
            adr_refs,
            "ADR-",
            &existing,
            violations,
        );
        let superseded = adr.status() == "superseded";
        if superseded && adr.items("superseded-by").is_empty() {
            violations.push(format!(
                "{}: superseded but `superseded-by` is empty",
                adr.path
            ));
        }
        if !superseded && !adr.items("superseded-by").is_empty() {
            violations.push(format!(
                "{}: status {:?} but `superseded-by` is set",
                adr.path,
                adr.status()
            ));
        }
        superseded_by.insert(adr.number, by);
        let old = front_refs(adr, "supersedes", adr_refs, "ADR-", &existing, violations);
        supersedes.insert(adr.number, old);
    }
    for (new, olds) in &supersedes {
        for old in olds {
            if !superseded_by.get(old).is_some_and(|s| s.contains(new)) {
                violations.push(format!(
                    "ADR-{new:04} supersedes ADR-{old:04}, but ADR-{old:04} is not superseded-by ADR-{new:04}"
                ));
            }
        }
    }
    for (old, news) in &superseded_by {
        for new in news {
            if !supersedes.get(new).is_some_and(|s| s.contains(old)) {
                violations.push(format!(
                    "ADR-{old:04} is superseded-by ADR-{new:04}, but ADR-{new:04} does not supersede ADR-{old:04}"
                ));
            }
        }
    }
}

/// 閉じた issue の日付と解決（4）と、`blocked-by` の参照先。
fn check_issue_states(issues: &[Doc], violations: &mut Vec<String>) {
    let existing: BTreeSet<u32> = issues.iter().map(|d| d.number).collect();
    for issue in issues {
        let closed = issue.status() == "closed";
        for key in ["closed", "resolution"] {
            let set = !issue.items(key).is_empty();
            if closed && !set {
                violations.push(format!("{}: closed but `{key}` is empty", issue.path));
            }
            if !closed && set {
                violations.push(format!(
                    "{}: status {:?} but `{key}` is set",
                    issue.path,
                    issue.status()
                ));
            }
        }
        front_refs(issue, "blocked-by", issue_refs, "#", &existing, violations);
    }
}

/// ADR の前付けの `issues` と、issue の「関係する ADR:」の行が両方向にそろっているか（5）。
fn check_adr_issue_links(adrs: &[Doc], issues: &[Doc], violations: &mut Vec<String>) {
    let existing: BTreeSet<u32> = issues.iter().map(|d| d.number).collect();
    let mut from_adrs: BTreeMap<u32, BTreeSet<u32>> = BTreeMap::new();
    for adr in adrs {
        for issue in front_refs(adr, "issues", issue_refs, "#", &existing, violations) {
            from_adrs.entry(issue).or_default().insert(adr.number);
        }
    }
    for issue in issues {
        let lines: Vec<&str> = issue
            .text
            .lines()
            .filter_map(|l| l.strip_prefix(RELATED_ADRS))
            .collect();
        let [line] = lines.as_slice() else {
            violations.push(format!(
                "{}: expected one `{RELATED_ADRS}` line, found {}",
                issue.path,
                lines.len()
            ));
            continue;
        };
        let listed: BTreeSet<u32> = adr_refs(line).into_iter().collect();
        let none = line.trim_start().starts_with("なし");
        if listed.is_empty() != none {
            violations.push(format!(
                "{}: `{RELATED_ADRS}` must list ADRs or say なし: {line:?}",
                issue.path
            ));
        }
        let expected = from_adrs.remove(&issue.number).unwrap_or_default();
        for adr in listed.difference(&expected) {
            violations.push(format!(
                "{}: lists ADR-{adr:04}, whose `issues` does not name #{:04}",
                issue.path, issue.number
            ));
        }
        for adr in expected.difference(&listed) {
            violations.push(format!(
                "{}: ADR-{adr:04} names #{:04} in `issues`, but `{RELATED_ADRS}` does not list it",
                issue.path, issue.number
            ));
        }
    }
}

/// 本文の参照 `ADR-NNNN` か `#NNNN` の参照先が存在するか（6）。
fn check_mentions(
    root: &Path,
    files: &[PathBuf],
    refs: fn(&str) -> Vec<u32>,
    prefix: &str,
    existing: &BTreeSet<u32>,
    violations: &mut Vec<String>,
) -> Result<()> {
    for path in files {
        let text = fs::read_to_string(path)?;
        for (i, line) in text.lines().enumerate() {
            for n in refs(line) {
                if !existing.contains(&n) {
                    violations.push(format!(
                        "{}:{}: {prefix}{n:04} does not exist",
                        relative(root, path),
                        i + 1
                    ));
                }
            }
        }
    }
    Ok(())
}

/// Markdown の相対リンクの行き先が存在するか（7）。
fn check_links(root: &Path, files: &[PathBuf], violations: &mut Vec<String>) -> Result<()> {
    for path in files {
        let text = fs::read_to_string(path)?;
        let dir = path.parent().unwrap_or(root);
        for (i, line) in text.lines().enumerate() {
            for link in relative_links(line) {
                if !dir.join(&link).exists() {
                    violations.push(format!(
                        "{}:{}: link target {link} does not exist",
                        relative(root, path),
                        i + 1
                    ));
                }
            }
        }
    }
    Ok(())
}

#[test]
fn docs_are_consistent() -> Result<()> {
    let root = repo_root();
    let mut violations = Vec::new();

    let adrs = load_numbered(&root, "docs/adr", &mut violations)?;
    let issues = load_numbered(&root, "docs/issues", &mut violations)?;
    check_unique(&adrs, &mut violations);
    check_unique(&issues, &mut violations);
    check_keys(&adrs, &ADR_KEYS, &ADR_STATUSES, &mut violations);
    check_keys(&issues, &ISSUE_KEYS, &ISSUE_STATUSES, &mut violations);
    check_supersession(&adrs, &mut violations);
    check_issue_states(&issues, &mut violations);
    check_adr_issue_links(&adrs, &issues, &mut violations);

    let claude_md = root.join("CLAUDE.md");
    let mut docs_md = Vec::new();
    files_under(&root.join("docs"), &["md"], &mut docs_md)?;
    let mut examples_md = Vec::new();
    files_under(&root.join("examples"), &["md"], &mut examples_md)?;

    // ADR の参照: docs、examples（文書とスクリプト）、crates のソース、CLAUDE.md。crates/docs は、
    // 存在しない番号を検査の材料として書くので除く
    let mut adr_files = vec![claude_md.clone()];
    adr_files.extend(docs_md.iter().cloned());
    files_under(
        &root.join("examples"),
        &["md", "py", "rs", "toml", "sh"],
        &mut adr_files,
    )?;
    let docs_crate = root.join("crates").join("docs");
    let mut crate_files = Vec::new();
    files_under(&root.join("crates"), &["rs", "toml"], &mut crate_files)?;
    adr_files.extend(
        crate_files
            .into_iter()
            .filter(|p| !p.starts_with(&docs_crate)),
    );
    let adr_numbers: BTreeSet<u32> = adrs.iter().map(|d| d.number).collect();
    check_mentions(
        &root,
        &adr_files,
        adr_refs,
        "ADR-",
        &adr_numbers,
        &mut violations,
    )?;

    // issue の参照: docs、examples の検証の記録、CLAUDE.md
    let mut issue_files = vec![claude_md.clone()];
    issue_files.extend(docs_md.iter().cloned());
    issue_files.extend(
        examples_md
            .iter()
            .filter(|p| p.file_name().is_some_and(|n| n == "verification.md"))
            .cloned(),
    );
    let issue_numbers: BTreeSet<u32> = issues.iter().map(|d| d.number).collect();
    check_mentions(
        &root,
        &issue_files,
        issue_refs,
        "#",
        &issue_numbers,
        &mut violations,
    )?;

    // 相対リンク: docs、examples、CLAUDE.md の Markdown
    let mut link_files = vec![claude_md];
    link_files.extend(docs_md);
    link_files.extend(examples_md);
    check_links(&root, &link_files, &mut violations)?;

    assert!(
        violations.is_empty(),
        "{} documentation violation(s):\n{}",
        violations.len(),
        violations.join("\n")
    );
    Ok(())
}
