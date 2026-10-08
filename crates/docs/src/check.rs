//! 文書の約束の検査（#0013）。リポジトリのファイルをパスと中身の組で受け取り、違反をすべて返す。
//!
//! 約束は docs/adr/README.md と docs/issues/README.md が決めている。

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::{
    FrontMatter, Link, Value, adr_refs, issue_refs, links, markdown_lines, numbered_file,
    parse_front_matter,
};

/// 約束に反しているところ。
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Violation {
    /// リポジトリの根からのパス（`/` 区切り）。
    pub path: String,
    /// 行番号（1 から）。ファイル全体についての違反なら `None`。
    pub line: Option<usize>,
    pub message: String,
}

impl fmt::Display for Violation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.line {
            Some(line) => write!(f, "{}:{line}: {}", self.path, self.message),
            None => write!(f, "{}: {}", self.path, self.message),
        }
    }
}

/// 検査するリポジトリのファイル。
#[derive(Debug, Clone, Default)]
pub struct Docs {
    texts: BTreeMap<String, String>,
    paths: BTreeSet<String>,
}

impl Docs {
    /// `texts` は中身を読むファイル（リポジトリの根からの `/` 区切りのパスと中身）、`other_paths` は中身を読まない
    /// ほかのファイル（リンク先として、あるかどうかだけを見る）。
    pub fn new(
        texts: BTreeMap<String, String>,
        other_paths: impl IntoIterator<Item = String>,
    ) -> Self {
        let mut paths: BTreeSet<String> = other_paths.into_iter().collect();
        paths.extend(texts.keys().cloned());
        Docs { texts, paths }
    }

    /// すべての約束を確かめ、違反を並べて返す。
    pub fn check(&self) -> Vec<Violation> {
        let mut out = Violations::default();
        let lines: BTreeMap<&str, Vec<(usize, &str)>> = self
            .texts
            .iter()
            .map(|(path, text)| (path.as_str(), scanned_lines(path, text, &mut out)))
            .collect();
        let adrs = self.numbered(ADR_DIR, &lines, &mut out);
        let issues = self.numbered(ISSUE_DIR, &lines, &mut out);
        let adr_numbers: BTreeSet<u32> = adrs.iter().map(|d| d.number).collect();
        let issue_numbers: BTreeSet<u32> = issues.iter().map(|d| d.number).collect();

        check_unique(&adrs, &mut out);
        check_unique(&issues, &mut out);
        check_front(&adrs, &ADR_KEYS, &ADR_STATUSES, "date", &mut out);
        check_front(&issues, &ISSUE_KEYS, &ISSUE_STATUSES, "opened", &mut out);
        check_headings(&adrs, &mut out);
        check_supersession(&adrs, &adr_numbers, &mut out);
        check_issue_states(&issues, &issue_numbers, &mut out);
        check_adr_issue_links(&adrs, &issues, &issue_numbers, &mut out);
        check_mentions(&lines, &adr_numbers, &issue_numbers, &mut out);
        self.check_links(&lines, &mut out);

        let mut found = out.0;
        found.sort();
        // 閉じていないフロントマターは、フロントマターと本文の両方の読み取りで見つかる
        found.dedup();
        found
    }

    /// `dir` の直下の番号つきの文書（README.md を除く）。
    fn numbered<'a>(
        &'a self,
        dir: &str,
        lines: &BTreeMap<&'a str, Vec<(usize, &'a str)>>,
        out: &mut Violations,
    ) -> Vec<Doc<'a>> {
        let mut docs = Vec::new();
        for (path, text) in &self.texts {
            let Some(name) = path.strip_prefix(dir) else {
                continue;
            };
            if name.contains('/') || !name.ends_with(".md") || name == "README.md" {
                continue;
            }
            let Some(number) = numbered_file(name) else {
                out.file(path, "file name is not NNNN-kebab-case.md");
                continue;
            };
            match parse_front_matter(text) {
                Ok(front) => docs.push(Doc {
                    path,
                    number,
                    front,
                    lines: lines.get(path.as_str()).cloned().unwrap_or_default(),
                }),
                Err(e) => out.file(path, e),
            }
        }
        docs
    }

    /// Markdown の相対リンク先があるか。絶対パスと、http、https、mailto 以外のスキームは違反にする。
    fn check_links(&self, lines: &BTreeMap<&str, Vec<(usize, &str)>>, out: &mut Violations) {
        for (path, lines) in lines {
            if !scanned_for_links(path) {
                continue;
            }
            let dir = path.rsplit_once('/').map_or("", |(dir, _)| dir);
            for (number, line) in lines {
                for link in links(line) {
                    match link {
                        Link::Relative(target) => {
                            if !resolve(dir, &target).is_some_and(|p| self.exists(&p)) {
                                out.line(
                                    path,
                                    *number,
                                    format!("link target {target} does not exist"),
                                );
                            }
                        }
                        Link::Absolute(target) => {
                            out.line(path, *number, format!("absolute link {target}"))
                        }
                        Link::OtherScheme(target) => out.line(
                            path,
                            *number,
                            format!("link with unsupported scheme {target}"),
                        ),
                    }
                }
            }
        }
    }

    /// ファイルか、ファイルを含むディレクトリがあるか（空のパスはリポジトリの根）。
    fn exists(&self, path: &str) -> bool {
        if path.is_empty() || self.paths.contains(path) {
            return true;
        }
        let dir = format!("{path}/");
        self.paths
            .range(dir.clone()..)
            .next()
            .is_some_and(|p| p.starts_with(&dir))
    }
}

const ADR_DIR: &str = "docs/adr/";
const ISSUE_DIR: &str = "docs/issues/";
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

/// 集めている違反。
#[derive(Default)]
struct Violations(Vec<Violation>);

impl Violations {
    fn file(&mut self, path: &str, message: impl ToString) {
        self.push(path, None, message);
    }

    fn line(&mut self, path: &str, line: usize, message: impl ToString) {
        self.push(path, Some(line), message);
    }

    fn push(&mut self, path: &str, line: Option<usize>, message: impl ToString) {
        self.0.push(Violation {
            path: path.to_string(),
            line,
            message: message.to_string(),
        });
    }
}

/// 番号つきの ADR か issue の1件。
struct Doc<'a> {
    path: &'a str,
    number: u32,
    front: FrontMatter,
    /// フロントマターとコードブロックの外の行。
    lines: Vec<(usize, &'a str)>,
}

impl Doc<'_> {
    /// 一つの文字列の `status`。ないか一覧なら空。
    fn status(&self) -> &str {
        match self.front.get("status") {
            Some(Value::Scalar(s)) => s,
            _ => "",
        }
    }

    fn items(&self, key: &str) -> Vec<&str> {
        self.front.get(key).map(Value::items).unwrap_or_default()
    }
}

fn extension(path: &str) -> &str {
    let name = path.rsplit('/').next().unwrap_or(path);
    name.rsplit_once('.').map_or("", |(_, e)| e)
}

/// 参照 `ADR-NNNN` と `#NNNN` を確かめるファイル: CLAUDE.md、docs の Markdown、examples の文書とスクリプト、
/// crates のソースと Cargo.toml。
fn scanned_for_refs(path: &str) -> bool {
    let ext = extension(path);
    path == "CLAUDE.md"
        || (path.starts_with("docs/") && ext == "md")
        || (path.starts_with("examples/") && ["md", "py", "rs", "toml", "sh"].contains(&ext))
        || (path.starts_with("crates/") && ["rs", "toml"].contains(&ext))
}

/// 相対リンクを確かめるファイル: CLAUDE.md と、docs と examples の Markdown。
fn scanned_for_links(path: &str) -> bool {
    path == "CLAUDE.md"
        || ((path.starts_with("docs/") || path.starts_with("examples/")) && extension(path) == "md")
}

/// 参照とリンクを探す行。Markdown はフロントマターとコードブロックを除き、ほかのファイルはすべての行。
fn scanned_lines<'a>(path: &str, text: &'a str, out: &mut Violations) -> Vec<(usize, &'a str)> {
    if extension(path) != "md" {
        return text.lines().enumerate().map(|(i, l)| (i + 1, l)).collect();
    }
    markdown_lines(text).unwrap_or_else(|e| {
        out.file(path, e);
        Vec::new()
    })
}

/// `dir` からの相対パス `target` を、根からのパスにする。根より上に出れば `None`。
fn resolve(dir: &str, target: &str) -> Option<String> {
    let mut parts: Vec<&str> = dir.split('/').filter(|p| !p.is_empty()).collect();
    for part in target.split('/') {
        match part {
            "" | "." => {}
            ".." => {
                parts.pop()?;
            }
            _ => parts.push(part),
        }
    }
    Some(parts.join("/"))
}

/// フロントマターの値の表示（違反の文に入れる）。
fn show(value: &Value) -> String {
    match value {
        Value::Empty => "empty".to_string(),
        Value::Scalar(s) => format!("{s:?}"),
        Value::List(items) => format!("{items:?}"),
    }
}

/// `YYYY-MM-DD`（月は 01〜12、日は 01〜31）か。
fn is_date(s: &str) -> bool {
    let bytes = s.as_bytes();
    let digits = |r: std::ops::Range<usize>| bytes[r].iter().all(u8::is_ascii_digit);
    if bytes.len() != 10
        || bytes[4] != b'-'
        || bytes[7] != b'-'
        || !digits(0..4)
        || !digits(5..7)
        || !digits(8..10)
    {
        return false;
    }
    let month: u32 = s[5..7].parse().unwrap_or(0);
    let day: u32 = s[8..10].parse().unwrap_or(0);
    (1..=12).contains(&month) && (1..=31).contains(&day)
}

/// 番号の重複（1）。
fn check_unique(docs: &[Doc], out: &mut Violations) {
    let mut by_number: BTreeMap<u32, Vec<&str>> = BTreeMap::new();
    for doc in docs {
        by_number.entry(doc.number).or_default().push(doc.path);
    }
    for (number, paths) in by_number {
        if let [first, rest @ ..] = paths.as_slice()
            && !rest.is_empty()
        {
            out.file(
                first,
                format!("number {number:04} is also used by {}", rest.join(", ")),
            );
        }
    }
}

/// フロントマターのキー、`status` の値、日付（2）。
fn check_front(
    docs: &[Doc],
    keys: &[&str],
    statuses: &[&str],
    date_key: &str,
    out: &mut Violations,
) {
    for doc in docs {
        for key in keys {
            if doc.front.get(key).is_none() {
                out.file(doc.path, format!("front matter has no `{key}`"));
            }
        }
        for key in doc.front.keys() {
            if !keys.contains(&key) {
                out.file(doc.path, format!("unexpected front matter key `{key}`"));
            }
        }
        match doc.front.get("status") {
            Some(Value::Scalar(s)) if statuses.contains(&s.as_str()) => {}
            Some(value) => out.file(
                doc.path,
                format!(
                    "`status` must be one of {}, not {}",
                    statuses.join(", "),
                    show(value)
                ),
            ),
            None => {}
        }
        check_date(doc, date_key, out);
    }
}

/// `key` があれば、一つの `YYYY-MM-DD` か。
fn check_date(doc: &Doc, key: &str, out: &mut Violations) {
    match doc.front.get(key) {
        Some(Value::Scalar(s)) if is_date(s) => {}
        Some(value) => out.file(
            doc.path,
            format!("`{key}` is not YYYY-MM-DD: {}", show(value)),
        ),
        None => {}
    }
}

/// ADR の最初の見出しが `# ADR-NNNN: 題名` で、番号がファイルの番号と同じか。
fn check_headings(adrs: &[Doc], out: &mut Violations) {
    for adr in adrs {
        let prefix = format!("# ADR-{:04}: ", adr.number);
        let heading = adr.lines.iter().find(|(_, l)| l.starts_with("# "));
        let ok = heading.is_some_and(|(_, l)| {
            l.strip_prefix(&prefix)
                .is_some_and(|t| !t.trim().is_empty())
        });
        if !ok {
            let message = format!("first heading must be `{}<title>`", prefix);
            match heading {
                Some((number, _)) => out.line(adr.path, *number, message),
                None => out.file(adr.path, message),
            }
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

/// フロントマターの `key` の要素を番号にする。形が違うものと、存在しない番号は違反にする。
fn front_refs(
    doc: &Doc,
    key: &str,
    refs: fn(&str) -> Vec<u32>,
    prefix: &str,
    existing: &BTreeSet<u32>,
    out: &mut Violations,
) -> BTreeSet<u32> {
    let mut numbers = BTreeSet::new();
    for item in doc.items(key) {
        match single_ref(item, refs, prefix) {
            Some(n) if existing.contains(&n) => {
                numbers.insert(n);
            }
            Some(_) => out.file(
                doc.path,
                format!("`{key}` names {item}, which does not exist"),
            ),
            None => out.file(
                doc.path,
                format!("`{key}` item {item:?} is not {prefix}NNNN"),
            ),
        }
    }
    numbers
}

/// 置き換えの状態と、両方向のリンク（3）。
fn check_supersession(adrs: &[Doc], existing: &BTreeSet<u32>, out: &mut Violations) {
    let mut supersedes: BTreeMap<u32, (&str, BTreeSet<u32>)> = BTreeMap::new();
    let mut superseded_by: BTreeMap<u32, (&str, BTreeSet<u32>)> = BTreeMap::new();
    for adr in adrs {
        let superseded = adr.status() == "superseded";
        let by_set = !adr.items("superseded-by").is_empty();
        if superseded && !by_set {
            out.file(
                adr.path,
                "status is superseded but `superseded-by` is empty",
            );
        }
        if !superseded && by_set {
            out.file(
                adr.path,
                "`superseded-by` is set but status is not superseded",
            );
        }
        let by = front_refs(adr, "superseded-by", adr_refs, "ADR-", existing, out);
        superseded_by.insert(adr.number, (adr.path, by));
        let old = front_refs(adr, "supersedes", adr_refs, "ADR-", existing, out);
        supersedes.insert(adr.number, (adr.path, old));
    }
    for (new, (path, olds)) in &supersedes {
        for old in olds {
            if !superseded_by.get(old).is_some_and(|(_, s)| s.contains(new)) {
                out.file(
                    path,
                    format!("supersedes ADR-{old:04}, but ADR-{old:04} is not superseded-by ADR-{new:04}"),
                );
            }
        }
    }
    for (old, (path, news)) in &superseded_by {
        for new in news {
            if !supersedes.get(new).is_some_and(|(_, s)| s.contains(old)) {
                out.file(
                    path,
                    format!("superseded-by ADR-{new:04}, but ADR-{new:04} does not supersede ADR-{old:04}"),
                );
            }
        }
    }
}

/// 閉じた issue の日付と解決（4）と、`blocked-by` の参照先。
fn check_issue_states(issues: &[Doc], existing: &BTreeSet<u32>, out: &mut Violations) {
    for issue in issues {
        let closed = issue.status() == "closed";
        for key in ["closed", "resolution"] {
            let set = !issue.items(key).is_empty();
            if closed && !set {
                out.file(issue.path, format!("status is closed but `{key}` is empty"));
            }
            if !closed && set {
                out.file(
                    issue.path,
                    format!("status is {} but `{key}` is set", issue.status()),
                );
            }
        }
        if closed && !issue.items("closed").is_empty() {
            check_date(issue, "closed", out);
        }
        front_refs(issue, "blocked-by", issue_refs, "#", existing, out);
    }
}

/// ADR のフロントマターの `issues` と、issue の「関係する ADR:」の行が両方向にそろっているか（5）。
fn check_adr_issue_links(
    adrs: &[Doc],
    issues: &[Doc],
    existing: &BTreeSet<u32>,
    out: &mut Violations,
) {
    let mut from_adrs: BTreeMap<u32, BTreeSet<u32>> = BTreeMap::new();
    for adr in adrs {
        for issue in front_refs(adr, "issues", issue_refs, "#", existing, out) {
            from_adrs.entry(issue).or_default().insert(adr.number);
        }
    }
    for issue in issues {
        let related: Vec<(usize, &str)> = issue
            .lines
            .iter()
            .filter_map(|(n, l)| l.strip_prefix(RELATED_ADRS).map(|rest| (*n, rest)))
            .collect();
        let [(number, line)] = related.as_slice() else {
            out.file(
                issue.path,
                format!(
                    "expected one `{RELATED_ADRS}` line outside code blocks, found {}",
                    related.len()
                ),
            );
            continue;
        };
        let listed: BTreeSet<u32> = adr_refs(line).into_iter().collect();
        if listed.is_empty() != line.trim_start().starts_with("なし") {
            out.line(
                issue.path,
                *number,
                format!("`{RELATED_ADRS}` must list ADRs or say なし"),
            );
        }
        let expected = from_adrs.get(&issue.number).cloned().unwrap_or_default();
        for adr in listed.difference(&expected) {
            out.line(
                issue.path,
                *number,
                format!(
                    "lists ADR-{adr:04}, whose `issues` does not name #{:04}",
                    issue.number
                ),
            );
        }
        for adr in expected.difference(&listed) {
            out.line(
                issue.path,
                *number,
                format!(
                    "ADR-{adr:04} names #{:04} in `issues`, but `{RELATED_ADRS}` does not list it",
                    issue.number
                ),
            );
        }
    }
}

/// 本文の参照 `ADR-NNNN` と `#NNNN` の参照先があるか（6）。
fn check_mentions(
    lines: &BTreeMap<&str, Vec<(usize, &str)>>,
    adrs: &BTreeSet<u32>,
    issues: &BTreeSet<u32>,
    out: &mut Violations,
) {
    for (path, lines) in lines {
        if !scanned_for_refs(path) {
            continue;
        }
        for (number, line) in lines {
            for n in adr_refs(line).into_iter().filter(|n| !adrs.contains(n)) {
                out.line(path, *number, format!("ADR-{n:04} does not exist"));
            }
            for n in issue_refs(line).into_iter().filter(|n| !issues.contains(n)) {
                out.line(path, *number, format!("#{n:04} does not exist"));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // 検査の対象のソースに参照と同じ形の文字列を書かないよう、参照は組み立てる
    fn adr(n: u32) -> String {
        format!("ADR-{n:04}")
    }

    fn issue(n: u32) -> String {
        format!("#{n:04}")
    }

    /// ADR の本文。見出しは 9 行目。
    fn adr_text(
        n: u32,
        status: &str,
        supersedes: &str,
        superseded_by: &str,
        issues: &str,
    ) -> String {
        let lines = [
            "---".to_string(),
            format!("status: {status}"),
            "date: 2026-10-08".to_string(),
            format!("supersedes: {supersedes}"),
            format!("superseded-by: {superseded_by}"),
            format!("issues: {issues}"),
            "---".to_string(),
            String::new(),
            format!("# {}: decide something", adr(n)),
            String::new(),
            "## 背景".to_string(),
        ];
        lines.map(|l| l.trim_end().to_string()).join("\n") + "\n"
    }

    /// issue の本文。「関係する ADR:」は 15 行目。
    fn issue_text(
        status: &str,
        closed: &str,
        blocked_by: &str,
        resolution: &str,
        related: &str,
    ) -> String {
        let lines = [
            "---".to_string(),
            format!("status: {status}"),
            "labels: [docs]".to_string(),
            "opened: 2026-10-08".to_string(),
            format!("closed: {closed}"),
            "milestone:".to_string(),
            format!("blocked-by: {blocked_by}"),
            format!("resolution: {resolution}"),
            "---".to_string(),
            String::new(),
            "# Something to do".to_string(),
            String::new(),
            "## 問題".to_string(),
            String::new(),
            format!("関係する ADR: {related}"),
        ];
        lines.map(|l| l.trim_end().to_string()).join("\n") + "\n"
    }

    const ADR1: &str = "docs/adr/0001-keep-a.md";
    const ADR2: &str = "docs/adr/0002-keep-b.md";
    const ISSUE1: &str = "docs/issues/0001-open-one.md";
    const ISSUE2: &str = "docs/issues/0002-closed-two.md";

    /// 約束をすべて満たす小さなリポジトリ。ADR-1 は ADR-2 に置き換えられ、ADR-2 は issue 1 に結ぶ。
    fn base() -> BTreeMap<String, String> {
        let mut files = BTreeMap::new();
        let mut add = |path: &str, text: String| files.insert(path.to_string(), text);
        add(ADR1, adr_text(1, "superseded", "[]", &adr(2), "[]"));
        add(
            ADR2,
            adr_text(
                2,
                "accepted",
                &format!("[{}]", adr(1)),
                "",
                &format!("[\"{}\"]", issue(1)),
            ),
        );
        add(
            ISSUE1,
            issue_text(
                "open",
                "",
                "[]",
                "",
                &format!("[{}](../adr/0002-keep-b.md)", adr(2)),
            ),
        );
        add(
            ISSUE2,
            issue_text(
                "closed",
                "2026-10-08",
                &format!("[\"{}\"]", issue(1)),
                "wontfix: not needed",
                "なし",
            ),
        );
        // コードブロックとコードスパンの中の参照とリンクは読まない
        add(
            "docs/adr/README.md",
            format!(
                "# ADR\n\n```markdown\nsupersedes: [{}]  # {}\n[example](missing.md)\n```\n\n\
                 Write `[x](nowhere.md)` and see [plans](../plans/a.md#part).\n",
                adr(9999),
                issue(9999)
            ),
        );
        add(
            "docs/plans/a.md",
            "# Plan\n\n[back][adrs]\n\n[adrs]: ../adr/\n".to_string(),
        );
        add(
            "CLAUDE.md",
            format!(
                "See [{}](docs/adr/0002-keep-b.md), {} and [patch](docs/patches/a.patch).\n",
                adr(2),
                issue(2)
            ),
        );
        add(
            "crates/x/src/lib.rs",
            format!("//! {} and {}\n", adr(1), issue(1)),
        );
        add("crates/x/Cargo.toml", format!("# {}\n", adr(2)));
        add("examples/v/run.py", format!("# {} {}\n", adr(2), issue(2)));
        add("examples/v/run.sh", format!("# {}\n", issue(1)));
        files
    }

    fn check(files: BTreeMap<String, String>) -> Vec<String> {
        let docs = Docs::new(files, ["docs/patches/a.patch".to_string()]);
        let mut found: Vec<String> = docs.check().iter().map(ToString::to_string).collect();
        found.sort();
        found
    }

    fn assert_violations(files: BTreeMap<String, String>, expected: &[String]) {
        let mut expected = expected.to_vec();
        expected.sort();
        assert_eq!(check(files), expected);
    }

    fn edit(files: &mut BTreeMap<String, String>, path: &str, from: &str, to: &str) {
        let text = files.get_mut(path).unwrap();
        assert!(text.contains(from), "{path} has no {from:?}");
        *text = text.replacen(from, to, 1);
    }

    fn add(files: &mut BTreeMap<String, String>, path: &str, text: String) {
        files.insert(path.to_string(), text);
    }

    #[test]
    fn consistent_docs_pass() {
        assert_violations(base(), &[]);
    }

    #[test]
    fn finds_duplicate_adr_numbers() {
        let mut files = base();
        let copy = files[ADR2].clone();
        add(&mut files, "docs/adr/0002-keep-c.md", copy);
        assert_violations(
            files,
            &[format!(
                "{ADR2}: number 0002 is also used by docs/adr/0002-keep-c.md"
            )],
        );
    }

    #[test]
    fn finds_duplicate_issue_numbers() {
        let mut files = base();
        let copy = files[ISSUE2].clone();
        add(&mut files, "docs/issues/0002-closed-again.md", copy);
        assert_violations(
            files,
            &[format!(
                "docs/issues/0002-closed-again.md: number 0002 is also used by {ISSUE2}"
            )],
        );
    }

    #[test]
    fn finds_bad_file_names() {
        let mut files = base();
        let copy = files[ADR1].clone();
        add(&mut files, "docs/adr/0003-Keep_C.md", copy);
        add(
            &mut files,
            "docs/issues/12-short.md",
            "# short\n".to_string(),
        );
        assert_violations(
            files,
            &[
                "docs/adr/0003-Keep_C.md: file name is not NNNN-kebab-case.md".to_string(),
                "docs/issues/12-short.md: file name is not NNNN-kebab-case.md".to_string(),
            ],
        );
    }

    #[test]
    fn reports_unreadable_front_matter() {
        let mut files = base();
        add(
            &mut files,
            "docs/issues/0003-no-front-matter.md",
            "# title\n".to_string(),
        );
        assert_violations(
            files,
            &["docs/issues/0003-no-front-matter.md: the file does not start with a `---` front matter"
                .to_string()],
        );
    }

    #[test]
    fn finds_missing_and_unexpected_keys() {
        let mut files = base();
        edit(&mut files, ADR2, "date: 2026-10-08\n", "title: x\n");
        edit(&mut files, ISSUE1, "milestone:\n", "");
        assert_violations(
            files,
            &[
                format!("{ADR2}: front matter has no `date`"),
                format!("{ADR2}: unexpected front matter key `title`"),
                format!("{ISSUE1}: front matter has no `milestone`"),
            ],
        );
    }

    #[test]
    fn finds_bad_statuses() {
        let mut files = base();
        edit(&mut files, ADR2, "status: accepted", "status: bogus");
        edit(&mut files, ISSUE1, "status: open", "status: [open]");
        assert_violations(
            files,
            &[
                format!(
                    "{ADR2}: `status` must be one of proposed, accepted, rejected, superseded, \
                     deprecated, not \"bogus\""
                ),
                format!("{ISSUE1}: `status` must be one of open, closed, not [\"open\"]"),
            ],
        );
    }

    #[test]
    fn finds_bad_dates() {
        let mut files = base();
        edit(&mut files, ADR2, "date: 2026-10-08", "date: 2026-13-08");
        edit(&mut files, ISSUE1, "opened: 2026-10-08", "opened: 8 Oct");
        edit(
            &mut files,
            ISSUE2,
            "closed: 2026-10-08",
            "closed: [2026-10-08]",
        );
        edit(&mut files, ADR1, "date: 2026-10-08", "date: 2026-02-30x");
        assert_violations(
            files,
            &[
                format!("{ADR1}: `date` is not YYYY-MM-DD: \"2026-02-30x\""),
                format!("{ADR2}: `date` is not YYYY-MM-DD: \"2026-13-08\""),
                format!("{ISSUE1}: `opened` is not YYYY-MM-DD: \"8 Oct\""),
                format!("{ISSUE2}: `closed` is not YYYY-MM-DD: [\"2026-10-08\"]"),
            ],
        );
    }

    #[test]
    fn checks_adr_headings() {
        let mut files = base();
        edit(
            &mut files,
            ADR1,
            &format!("# {}: decide", adr(1)),
            "# Decide",
        );
        edit(
            &mut files,
            ADR2,
            &format!("# {}: decide something", adr(2)),
            "text",
        );
        // コードブロックの中の見出しは数えない
        let fenced = format!("```\n# {}: in a fence\n```\n\n# {}: real", adr(3), adr(1));
        add(
            &mut files,
            "docs/adr/0003-keep-c.md",
            adr_text(3, "rejected", "[]", "", "[]")
                .replace(&format!("# {}: decide something", adr(3)), &fenced),
        );
        assert_violations(
            files,
            &[
                format!("{ADR1}:9: first heading must be `# {}: <title>`", adr(1)),
                format!("{ADR2}: first heading must be `# {}: <title>`", adr(2)),
                format!(
                    "docs/adr/0003-keep-c.md:13: first heading must be `# {}: <title>`",
                    adr(3)
                ),
            ],
        );
    }

    #[test]
    fn requires_superseded_by_on_superseded_adrs() {
        let mut files = base();
        edit(
            &mut files,
            ADR1,
            &format!("superseded-by: {}", adr(2)),
            "superseded-by:",
        );
        assert_violations(
            files,
            &[
                format!("{ADR1}: status is superseded but `superseded-by` is empty"),
                format!(
                    "{ADR2}: supersedes {}, but {} is not superseded-by {}",
                    adr(1),
                    adr(1),
                    adr(2)
                ),
            ],
        );
    }

    #[test]
    fn rejects_superseded_by_on_other_adrs() {
        let mut files = base();
        edit(
            &mut files,
            ADR2,
            "superseded-by:",
            &format!("superseded-by: {}", adr(1)),
        );
        assert_violations(
            files,
            &[
                format!("{ADR2}: `superseded-by` is set but status is not superseded"),
                format!(
                    "{ADR2}: superseded-by {}, but {} does not supersede {}",
                    adr(1),
                    adr(1),
                    adr(2)
                ),
            ],
        );
    }

    #[test]
    fn checks_supersession_targets() {
        let mut files = base();
        edit(
            &mut files,
            ADR2,
            &format!("supersedes: [{}]", adr(1)),
            &format!("supersedes: [{}, {}, adr-3]", adr(1), adr(99)),
        );
        assert_violations(
            files,
            &[
                format!(
                    "{ADR2}: `supersedes` names {}, which does not exist",
                    adr(99)
                ),
                format!("{ADR2}: `supersedes` item \"adr-3\" is not ADR-NNNN"),
            ],
        );
    }

    #[test]
    fn checks_closed_and_open_issues() {
        let mut files = base();
        edit(&mut files, ISSUE2, "closed: 2026-10-08", "closed:");
        edit(
            &mut files,
            ISSUE2,
            "resolution: wontfix: not needed",
            "resolution:",
        );
        edit(&mut files, ISSUE1, "resolution:", "resolution: wontfix: x");
        assert_violations(
            files,
            &[
                format!("{ISSUE2}: status is closed but `closed` is empty"),
                format!("{ISSUE2}: status is closed but `resolution` is empty"),
                format!("{ISSUE1}: status is open but `resolution` is set"),
            ],
        );
    }

    #[test]
    fn accepts_resolution_forms() {
        let hash40 = "0123456789abcdef0123456789abcdef01234567";
        for resolution in [
            adr(2),
            "コミット 36fbabc".to_string(),
            format!("コミット {hash40}"),
            "wontfix: not needed".to_string(),
        ] {
            let mut files = base();
            edit(
                &mut files,
                ISSUE2,
                "resolution: wontfix: not needed",
                &format!("resolution: {resolution}"),
            );
            assert_violations(files, &[]);
        }
    }

    #[test]
    fn rejects_other_resolutions() {
        let hash41 = "0123456789abcdef0123456789abcdef012345678";
        for resolution in [
            "crates/docs（cargo test）".to_string(),
            "コミット 4f9e59".to_string(),
            format!("コミット {hash41}"),
            "コミット 4F9E59A".to_string(),
            "コミット xyz1234".to_string(),
            "コミット".to_string(),
            "wontfix:".to_string(),
            format!("{} と {}", adr(1), adr(2)),
        ] {
            let mut files = base();
            edit(
                &mut files,
                ISSUE2,
                "resolution: wontfix: not needed",
                &format!("resolution: {resolution}"),
            );
            assert_violations(
                files,
                &[format!(
                    "{ISSUE2}: `resolution` must be ADR-NNNN, `コミット <hash>` or `wontfix: <reason>`, \
                     not {resolution:?}"
                )],
            );
        }
        let mut files = base();
        edit(
            &mut files,
            ISSUE2,
            "resolution: wontfix: not needed",
            "resolution: [wontfix: x]",
        );
        assert_violations(
            files,
            &[format!(
                "{ISSUE2}: `resolution` must be ADR-NNNN, `コミット <hash>` or `wontfix: <reason>`, \
                 not [\"wontfix: x\"]"
            )],
        );
    }

    #[test]
    fn checks_resolving_adrs_exist() {
        let mut files = base();
        edit(
            &mut files,
            ISSUE2,
            "resolution: wontfix: not needed",
            &format!("resolution: {}", adr(99)),
        );
        assert_violations(
            files,
            &[format!(
                "{ISSUE2}: `resolution` names {}, which does not exist",
                adr(99)
            )],
        );
    }

    #[test]
    fn checks_blocking_issues() {
        let mut files = base();
        edit(
            &mut files,
            ISSUE2,
            &format!("blocked-by: [\"{}\"]", issue(1)),
            &format!("blocked-by: [\"{}\", \"{}\"]", issue(1), issue(99)),
        );
        assert_violations(
            files,
            &[format!(
                "{ISSUE2}: `blocked-by` names {}, which does not exist",
                issue(99)
            )],
        );
    }

    #[test]
    fn finds_adrs_listed_only_on_the_issue() {
        let mut files = base();
        edit(
            &mut files,
            ISSUE1,
            "関係する ADR: ",
            &format!("関係する ADR: {}、", adr(1)),
        );
        assert_violations(
            files,
            &[format!(
                "{ISSUE1}:15: lists {}, whose `issues` does not name {}",
                adr(1),
                issue(1)
            )],
        );
    }

    #[test]
    fn finds_adrs_listed_only_in_front_matter() {
        let mut files = base();
        edit(
            &mut files,
            ADR1,
            "issues: []",
            &format!("issues: [\"{}\"]", issue(2)),
        );
        assert_violations(
            files,
            &[format!(
                "{ISSUE2}:15: {} names {} in `issues`, but `関係する ADR:` does not list it",
                adr(1),
                issue(2)
            )],
        );
    }

    #[test]
    fn requires_one_related_adr_line_outside_fences() {
        let mut files = base();
        edit(
            &mut files,
            ISSUE2,
            "関係する ADR: なし",
            "```\n関係する ADR: なし\n```",
        );
        let twice = files[ISSUE1].clone() + "\n関係する ADR: なし\n";
        add(&mut files, ISSUE1, twice);
        assert_violations(
            files,
            &[
                format!("{ISSUE1}: expected one `関係する ADR:` line outside code blocks, found 2"),
                format!("{ISSUE2}: expected one `関係する ADR:` line outside code blocks, found 0"),
            ],
        );
    }

    #[test]
    fn requires_adrs_or_none_on_the_related_line() {
        let mut files = base();
        edit(
            &mut files,
            ISSUE2,
            "関係する ADR: なし",
            "関係する ADR: 未定",
        );
        assert_violations(
            files,
            &[format!(
                "{ISSUE2}:15: `関係する ADR:` must list ADRs or say なし"
            )],
        );
    }

    #[test]
    fn finds_missing_references() {
        let mut files = base();
        add(
            &mut files,
            "examples/v/verification.md",
            format!("{} and {}\n\n```\n{}\n```\n", adr(99), issue(99), adr(98)),
        );
        edit(
            &mut files,
            "crates/x/src/lib.rs",
            "\n",
            &format!("\n// {} {}\n", issue(77), adr(77)),
        );
        edit(
            &mut files,
            "examples/v/run.py",
            "\n",
            &format!("\n# {}\n", adr(66)),
        );
        edit(&mut files, "CLAUDE.md", "\n", &format!("\n{}\n", issue(55)));
        assert_violations(
            files,
            &[
                format!("examples/v/verification.md:1: {} does not exist", adr(99)),
                format!("examples/v/verification.md:1: {} does not exist", issue(99)),
                format!("crates/x/src/lib.rs:2: {} does not exist", adr(77)),
                format!("crates/x/src/lib.rs:2: {} does not exist", issue(77)),
                format!("examples/v/run.py:2: {} does not exist", adr(66)),
                format!("CLAUDE.md:2: {} does not exist", issue(55)),
            ],
        );
    }

    #[test]
    fn reports_unclosed_fences() {
        let mut files = base();
        add(
            &mut files,
            "docs/plans/b.md",
            "# B\n\n```\ncode\n".to_string(),
        );
        assert_violations(
            files,
            &["docs/plans/b.md: code fence opened on line 3 is not closed".to_string()],
        );
    }

    #[test]
    fn checks_links() {
        let mut files = base();
        add(
            &mut files,
            "docs/plans/c.md",
            "[gone](gone.md) [abs](/docs/plans/a.md) [ftp](ftp://host/x) [mail](mailto:a@b) \
             [web](https://x.org/a.md) [here](#top) [up](../../CLAUDE.md#x) [dir](../adr)\n\
             [ref]: ../adr/missing.md\n\
             [ok]: a.md\n\
             `[code](nope.md)` and ``[a](b.md)``\n\
             [escape](../../../outside.md)\n"
                .to_string(),
        );
        add(
            &mut files,
            "examples/v/notes.md",
            "[x](run.py) [y](none.py)\n".to_string(),
        );
        assert_violations(
            files,
            &[
                "docs/plans/c.md:1: link target gone.md does not exist".to_string(),
                "docs/plans/c.md:1: absolute link /docs/plans/a.md".to_string(),
                "docs/plans/c.md:1: link with unsupported scheme ftp://host/x".to_string(),
                "docs/plans/c.md:2: link target ../adr/missing.md does not exist".to_string(),
                "docs/plans/c.md:5: link target ../../../outside.md does not exist".to_string(),
                "examples/v/notes.md:1: link target none.py does not exist".to_string(),
            ],
        );
    }
}
