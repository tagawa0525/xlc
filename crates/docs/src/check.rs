//! 文書の約束の検査（#0013）。リポジトリのファイルをパスと中身の組で受け取り、違反をすべて返す。
//!
//! 約束は docs/adr/README.md と docs/issues/README.md が決めている。

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

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
    /// ほかのファイル（リンクの行き先として、あるかどうかだけを見る）。
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
        Vec::new()
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
