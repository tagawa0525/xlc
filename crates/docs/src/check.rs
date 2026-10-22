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
        let adrs = self.numbered(ADR_DIR, &mut out);
        let issues = self.numbered(ISSUE_DIR, &mut out);
        let adr_numbers: BTreeSet<u32> = adrs.iter().map(|d| d.number).collect();
        let issue_numbers: BTreeSet<u32> = issues.iter().map(|d| d.number).collect();

        check_unique(&adrs, &mut out);
        check_unique(&issues, &mut out);
        check_front(
            &adrs,
            &ADR_KEYS,
            &ADR_OPTIONAL_KEYS,
            &ADR_STATUSES,
            "date",
            &mut out,
        );
        check_front(
            &issues,
            &ISSUE_KEYS,
            &[],
            &ISSUE_STATUSES,
            "opened",
            &mut out,
        );
        check_headings(&adrs, &mut out);
        check_supersession(&adrs, &adr_numbers, &mut out);
        check_amendment(&adrs, &adr_numbers, &mut out);
        for adr in &adrs {
            front_refs(adr, "requires", adr_refs, "ADR-", &adr_numbers, &mut out);
        }
        let adr_issues = adr_issues(&adrs, &issue_numbers, &mut out);
        check_issue_states(&issues, &issue_numbers, &adr_issues, &mut out);
        check_adr_issue_links(&adr_issues, &issues, &mut out);
        check_mentions(&lines, &adr_numbers, &issue_numbers, &mut out);
        self.check_links(&lines, &mut out);

        let mut found = out.0;
        found.sort();
        found
    }

    /// 閉じた issue の `resolution` が `コミット <hash>` の形のとき、その issue のパスと hash。git で hash が
    /// 指すコミットを確かめるのは呼ぶ側（tests/consistency.rs）。
    pub fn resolution_commits(&self) -> Vec<(String, String)> {
        self.numbered(ISSUE_DIR, &mut Violations::default())
            .iter()
            .filter(|issue| issue.status() == "closed")
            .filter_map(|issue| match issue.front.get("resolution") {
                Some(Value::Scalar(s)) => {
                    commit_hash(s).map(|h| (issue.path.to_string(), h.to_string()))
                }
                _ => None,
            })
            .collect()
    }

    /// `dir` の直下の番号つきの文書（README.md を除く）。
    fn numbered(&self, dir: &str, out: &mut Violations) -> Vec<Doc<'_>> {
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
                    // 閉じていないコードブロックは、行を読むほかの検査（scanned_lines）が違反にする
                    lines: markdown_lines(text).ok(),
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
const ADR_KEYS: [&str; 6] = [
    "status",
    "date",
    "supersedes",
    "superseded-by",
    "issues",
    "requires",
];
/// 書いたときだけ検査するキー。一部の置き換えは両方向に結ぶ（`check_amendment`）。
const ADR_OPTIONAL_KEYS: [&str; 2] = ["amends", "amended-by"];
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
    /// コードブロックの外の行（フロントマターを含む）。コードブロックが閉じていなければ `None` で、行を見る検査を飛ばす。
    lines: Option<Vec<(usize, &'a str)>>,
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

/// 参照とリンクを探す行。Markdown はコードブロックを除き、ほかのファイルはすべての行。コードブロックが閉じて
/// いなければ違反にし、そのファイルの行は見ない。
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

/// 暦にある日付の `YYYY-MM-DD` か（月の日数と閏年を見る）。
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
    let (Ok(year), Ok(month), Ok(day)) = (
        s[0..4].parse::<u32>(),
        s[5..7].parse::<u32>(),
        s[8..10].parse::<u32>(),
    ) else {
        return false;
    };
    let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
    let days = match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 if leap => 29,
        2 => 28,
        _ => return false,
    };
    (1..=days).contains(&day)
}

/// `コミット <hash>` の hash（7〜40 桁の英小文字の16進）。形が違えば `None`。
fn commit_hash(resolution: &str) -> Option<&str> {
    resolution.strip_prefix("コミット ").filter(|h| {
        (7..=40).contains(&h.len())
            && h.bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    })
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
    optional_keys: &[&str],
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
            if !keys.contains(&key) && !optional_keys.contains(&key) {
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
        let Some(lines) = &adr.lines else {
            continue;
        };
        let heading = lines.iter().find(|(_, l)| l.starts_with("# "));
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

/// 一部の置き換えの両方向のリンク。`amends` に挙げた ADR は `amended-by` にその ADR を挙げ、逆も同じ。
fn check_amendment(adrs: &[Doc], existing: &BTreeSet<u32>, out: &mut Violations) {
    let mut amends: BTreeMap<u32, (&str, BTreeSet<u32>)> = BTreeMap::new();
    let mut amended_by: BTreeMap<u32, (&str, BTreeSet<u32>)> = BTreeMap::new();
    for adr in adrs {
        let by = front_refs(adr, "amended-by", adr_refs, "ADR-", existing, out);
        amended_by.insert(adr.number, (adr.path, by));
        let old = front_refs(adr, "amends", adr_refs, "ADR-", existing, out);
        amends.insert(adr.number, (adr.path, old));
    }
    for (new, (path, olds)) in &amends {
        for old in olds {
            if !amended_by.get(old).is_some_and(|(_, s)| s.contains(new)) {
                out.file(
                    path,
                    format!("amends ADR-{old:04}, but ADR-{old:04} is not amended-by ADR-{new:04}"),
                );
            }
        }
    }
    for (old, (path, news)) in &amended_by {
        for new in news {
            if !amends.get(new).is_some_and(|(_, s)| s.contains(old)) {
                out.file(
                    path,
                    format!(
                        "amended-by ADR-{new:04}, but ADR-{new:04} does not amend ADR-{old:04}"
                    ),
                );
            }
        }
    }
}

/// 閉じた issue の日付と解決（4）と、`blocked-by` の参照先。
fn check_issue_states(
    issues: &[Doc],
    existing: &BTreeSet<u32>,
    adr_issues: &BTreeMap<u32, BTreeSet<u32>>,
    out: &mut Violations,
) {
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
        if closed {
            check_resolution(issue, adr_issues, out);
        }
        front_refs(issue, "blocked-by", issue_refs, "#", existing, out);
    }
}

/// `resolution` が `ADR-NNNN`（ある ADR）、`コミット <hash>`（7〜40 桁の英小文字の16進）、`wontfix: 理由` のどれか
/// （docs/issues/README.md）。ADR なら、その ADR の `issues` がこの issue を挙げていること（docs/adr/README.md）。
/// 空なら何もしない（空のことは状態の検査が見る）。
fn check_resolution(issue: &Doc, adr_issues: &BTreeMap<u32, BTreeSet<u32>>, out: &mut Violations) {
    let value = match issue.front.get("resolution") {
        None | Some(Value::Empty) => return,
        Some(Value::Scalar(s)) => Some(s.as_str()),
        Some(Value::List(_)) => None,
    };
    if let Some(n) = value.and_then(|s| single_ref(s, adr_refs, "ADR-")) {
        match adr_issues.get(&n) {
            None => out.file(
                issue.path,
                format!("`resolution` names ADR-{n:04}, which does not exist"),
            ),
            Some(named) if !named.contains(&issue.number) => out.file(
                issue.path,
                format!(
                    "`resolution` names ADR-{n:04}, whose `issues` does not name #{:04}",
                    issue.number
                ),
            ),
            Some(_) => {}
        }
        return;
    }
    let valid = value.is_some_and(|s| {
        commit_hash(s).is_some()
            || s.strip_prefix("wontfix:")
                .is_some_and(|reason| !reason.trim().is_empty())
    });
    if !valid {
        let shown = issue.front.get("resolution").map(show).unwrap_or_default();
        out.file(
            issue.path,
            format!(
                "`resolution` must be ADR-NNNN, `コミット <hash>` or `wontfix: <reason>`, not {shown}"
            ),
        );
    }
}

/// ADR ごとの、フロントマターの `issues` が挙げる issue。形が違うものと、ない issue は違反にする。
fn adr_issues(
    adrs: &[Doc],
    existing: &BTreeSet<u32>,
    out: &mut Violations,
) -> BTreeMap<u32, BTreeSet<u32>> {
    adrs.iter()
        .map(|adr| {
            let issues = front_refs(adr, "issues", issue_refs, "#", existing, out);
            (adr.number, issues)
        })
        .collect()
}

/// ADR のフロントマターの `issues` と、issue の「関係する ADR:」の行が両方向にそろっているか（5）。
fn check_adr_issue_links(
    adr_issues: &BTreeMap<u32, BTreeSet<u32>>,
    issues: &[Doc],
    out: &mut Violations,
) {
    let mut from_adrs: BTreeMap<u32, BTreeSet<u32>> = BTreeMap::new();
    for (adr, named) in adr_issues {
        for issue in named {
            from_adrs.entry(*issue).or_default().insert(*adr);
        }
    }
    for issue in issues {
        let Some(lines) = &issue.lines else {
            continue;
        };
        let related: Vec<(usize, &str)> = lines
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
mod tests;
