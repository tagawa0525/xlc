//! 文書（docs/adr、docs/issues ほか）の整合性の検査（#0013）。
//!
//! 約束（番号、フロントマター、置き換えのリンク、issue と ADR の対応、参照先、相対リンク）の検査は [`check`] にあり、
//! tests/consistency.rs がリポジトリのファイルを読んで実行する（ワークスペースの `cargo test`）。ここには、
//! Markdown のフロントマター、番号つきのファイル名、`ADR-NNNN` と `#NNNN` の参照、リンクを読む小さな関数を置く。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

pub mod check;

pub use check::{Docs, Violation};

use std::fmt;

/// フロントマターを読めなかった理由。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// 先頭の行が `---` でない。
    MissingFrontMatter,
    /// フロントマターを閉じる `---` がない。
    UnclosedFrontMatter,
    /// `key: value` の形でない行（フロントマターの中の行番号は 1 から）。
    MalformedLine { line: usize, text: String },
    /// 同じキーが2回ある。
    DuplicateKey(String),
    /// 閉じていない引用符や角括弧。
    MalformedValue { key: String, value: String },
    /// 閉じていないコードブロック（開いた行の番号は 1 から）。
    UnclosedFence { line: usize },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::MissingFrontMatter => {
                write!(f, "the file does not start with a `---` front matter")
            }
            Error::UnclosedFrontMatter => write!(f, "the front matter has no closing `---`"),
            Error::MalformedLine { line, text } => {
                write!(f, "front matter line {line} is not `key: value`: {text:?}")
            }
            Error::DuplicateKey(key) => write!(f, "front matter key `{key}` appears twice"),
            Error::MalformedValue { key, value } => {
                write!(f, "front matter value of `{key}` is malformed: {value:?}")
            }
            Error::UnclosedFence { line } => {
                write!(f, "code fence opened on line {line} is not closed")
            }
        }
    }
}

impl std::error::Error for Error {}

/// フロントマターの値。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Value {
    /// `key:` のあとが空。
    Empty,
    /// 裸の文字列か、引用符で囲んだ文字列（引用符は外す）。
    Scalar(String),
    /// `[a, "b"]` の形の一覧（各要素の引用符は外す）。
    List(Vec<String>),
}

impl Value {
    /// 値に含まれる文字列。空なら空、一つの文字列なら1つ、一覧ならその要素。
    pub fn items(&self) -> Vec<&str> {
        match self {
            Value::Empty => Vec::new(),
            Value::Scalar(s) => vec![s.as_str()],
            Value::List(items) => items.iter().map(String::as_str).collect(),
        }
    }
}

/// フロントマターのキーと値（書かれた順）。
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FrontMatter {
    entries: Vec<(String, Value)>,
}

impl FrontMatter {
    /// キーの値。キーがなければ `None`。
    pub fn get(&self, key: &str) -> Option<&Value> {
        self.entries.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }

    /// 書かれた順のキー。
    pub fn keys(&self) -> impl Iterator<Item = &str> {
        self.entries.iter().map(|(k, _)| k.as_str())
    }
}

/// Markdown の先頭の `---` で囲んだフロントマターを読む。
///
/// 各行は `key: value` で、値は空、裸の文字列、引用符で囲んだ文字列、`[a, "b"]` の一覧のどれか。引用符の外の
/// `#`（行頭か空白のあと）から行末まではコメントとして捨てる。空行、重複したキー、閉じていない引用符や角括弧は
/// エラーにする。
pub fn parse_front_matter(text: &str) -> Result<FrontMatter, Error> {
    let mut lines = text.lines();
    if lines.next() != Some("---") {
        return Err(Error::MissingFrontMatter);
    }
    let mut entries: Vec<(String, Value)> = Vec::new();
    for (i, line) in lines.enumerate() {
        if line == "---" {
            return Ok(FrontMatter { entries });
        }
        let malformed = || Error::MalformedLine {
            line: i + 1,
            text: line.to_string(),
        };
        let (key, value) = strip_comment(line).split_once(':').ok_or_else(malformed)?;
        let key = key.trim();
        if key.is_empty() || key.contains(char::is_whitespace) {
            return Err(malformed());
        }
        if entries.iter().any(|(k, _)| k == key) {
            return Err(Error::DuplicateKey(key.to_string()));
        }
        let value = parse_value(key, value.trim())?;
        entries.push((key.to_string(), value));
    }
    Err(Error::UnclosedFrontMatter)
}

/// 引用符の外で、行頭か空白のあとにある `#` から行末までを除く。
fn strip_comment(line: &str) -> &str {
    let mut in_quote = false;
    let mut previous: Option<char> = None;
    for (i, c) in line.char_indices() {
        match c {
            '"' => in_quote = !in_quote,
            '#' if !in_quote && previous.is_none_or(char::is_whitespace) => {
                return line[..i].trim_end();
            }
            _ => {}
        }
        previous = Some(c);
    }
    line.trim_end()
}

fn parse_value(key: &str, value: &str) -> Result<Value, Error> {
    let malformed = || Error::MalformedValue {
        key: key.to_string(),
        value: value.to_string(),
    };
    if value.is_empty() {
        return Ok(Value::Empty);
    }
    if let Some(rest) = value.strip_prefix('[') {
        let inner = rest.strip_suffix(']').ok_or_else(malformed)?.trim();
        if inner.is_empty() {
            return Ok(Value::List(Vec::new()));
        }
        return inner
            .split(',')
            .map(|item| match unquote(item.trim()) {
                Some(s) if !s.is_empty() => Ok(s.to_string()),
                _ => Err(malformed()),
            })
            .collect::<Result<Vec<_>, _>>()
            .map(Value::List);
    }
    unquote(value)
        .map(|s| Value::Scalar(s.to_string()))
        .ok_or_else(malformed)
}

/// `"text"` なら中身を、引用符を含まない裸の文字列ならそのままを返す。引用符が対にならなければ `None`。
fn unquote(s: &str) -> Option<&str> {
    let inner = match s.strip_prefix('"') {
        Some(rest) => rest.strip_suffix('"')?,
        None => s,
    };
    (!inner.contains('"')).then_some(inner)
}

/// `NNNN-name.md` の形のファイル名から番号を取り出す。形が違えば `None`。
pub fn numbered_file(name: &str) -> Option<u32> {
    let (number, rest) = name.strip_suffix(".md")?.split_once('-')?;
    if number.len() != 4 || !number.bytes().all(|b| b.is_ascii_digit()) || rest.is_empty() {
        return None;
    }
    number.parse().ok()
}

/// 英数字か `_`（ASCII のみ。日本語の文に続けて書いた参照は、参照として数える）。
fn is_word(c: char) -> bool {
    c.is_ascii_alphanumeric() || c == '_'
}

/// `text[start..]` の先頭がちょうど4桁の数字で、そのあとに英数字が続かなければ、その数を返す。
fn four_digits_at(text: &str, start: usize) -> Option<u32> {
    let digits = text.get(start..start + 4)?;
    if !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    if text[start + 4..].chars().next().is_some_and(is_word) {
        return None;
    }
    digits.parse().ok()
}

/// `ADR-NNNN`（4桁）の参照の番号を、現れた順に返す。
pub fn adr_refs(text: &str) -> Vec<u32> {
    const PREFIX: &str = "ADR-";
    text.match_indices(PREFIX)
        .filter(|(i, _)| !text[..*i].chars().next_back().is_some_and(is_word))
        .filter_map(|(i, _)| four_digits_at(text, i + PREFIX.len()))
        .collect()
}

/// issue の参照 `#NNNN`（ちょうど4桁で、前に `#` も英数字もない）の番号を、現れた順に返す。
pub fn issue_refs(text: &str) -> Vec<u32> {
    text.match_indices('#')
        .filter(|(i, _)| {
            !text[..*i]
                .chars()
                .next_back()
                .is_some_and(|c| c == '#' || is_word(c))
        })
        .filter_map(|(i, _)| four_digits_at(text, i + 1))
        .collect()
}

/// Markdown のリンク先。外部の URL（`http:`、`https:`、`mailto:`）と、同じ文書の中の `#anchor` だけの
/// ものは含めない。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Link {
    /// 相対パス（`#anchor` は外す）。
    Relative(String),
    /// `/` で始まる絶対パス。
    Absolute(String),
    /// `http:`、`https:`、`mailto:` 以外のスキーム（`ftp:`、`file:` など）。
    OtherScheme(String),
}

/// 1行の中のリンク先を、現れた順に返す。インラインのリンク `[text](target)` と、参照の定義
/// `[label]: target` を読む。コードスパンの中は読まない。
pub fn links(_line: &str) -> Vec<Link> {
    Vec::new()
}

/// 1行からコードスパン（対になったバッククォートで囲んだ部分）を除く。
pub fn strip_code_spans(line: &str) -> String {
    line.to_string()
}

/// Markdown の、フロントマターとコードブロック（```` ``` ```` か `~~~` で囲んだ部分）の外の行を、行番号（1 から）とともに
/// 返す。閉じていないフロントマターとコードブロックはエラーにする。
pub fn markdown_lines(_text: &str) -> Result<Vec<(usize, &str)>, Error> {
    Ok(Vec::new())
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

    fn scalar(s: &str) -> Value {
        Value::Scalar(s.to_string())
    }

    fn list(items: &[&str]) -> Value {
        Value::List(items.iter().map(|s| s.to_string()).collect())
    }

    #[test]
    fn reads_front_matter_values() {
        let text = format!(
            "---\n\
             status: superseded\n\
             date: 2026-10-08\n\
             supersedes: [{a2}, {a3}]\n\
             superseded-by: {a7}\n\
             issues: [\"{i1}\"]\n\
             closed:\n\
             title: \"quoted text\"\n\
             ---\n\
             \n\
             # {a1}: title\n",
            a1 = adr(1),
            a2 = adr(2),
            a3 = adr(3),
            a7 = adr(7),
            i1 = issue(1),
        );
        let fm = parse_front_matter(&text).unwrap();
        assert_eq!(fm.get("status"), Some(&scalar("superseded")));
        assert_eq!(fm.get("date"), Some(&scalar("2026-10-08")));
        assert_eq!(fm.get("supersedes"), Some(&list(&[&adr(2), &adr(3)])));
        assert_eq!(fm.get("superseded-by"), Some(&scalar(&adr(7))));
        assert_eq!(fm.get("issues"), Some(&list(&[&issue(1)])));
        assert_eq!(fm.get("closed"), Some(&Value::Empty));
        assert_eq!(fm.get("title"), Some(&scalar("quoted text")));
        assert_eq!(fm.get("missing"), None);
        assert_eq!(
            fm.keys().collect::<Vec<_>>(),
            [
                "status",
                "date",
                "supersedes",
                "superseded-by",
                "issues",
                "closed",
                "title"
            ]
        );
    }

    #[test]
    fn strips_comments_but_not_quoted_hashes() {
        let text = format!(
            "---\n\
             status: accepted      # proposed | accepted\n\
             supersedes: []        # 例: [{a2}]\n\
             superseded-by:        # 例: {a7}\n\
             issues: [\"{i1}\", \"{i14}\"] # 例\n\
             ---\n",
            a2 = adr(2),
            a7 = adr(7),
            i1 = issue(1),
            i14 = issue(14),
        );
        let fm = parse_front_matter(&text).unwrap();
        assert_eq!(fm.get("status"), Some(&scalar("accepted")));
        assert_eq!(fm.get("supersedes"), Some(&list(&[])));
        assert_eq!(fm.get("superseded-by"), Some(&Value::Empty));
        assert_eq!(fm.get("issues"), Some(&list(&[&issue(1), &issue(14)])));
    }

    #[test]
    fn value_items_flatten_scalars_and_lists() {
        assert!(Value::Empty.items().is_empty());
        assert_eq!(scalar("a").items(), ["a"]);
        assert_eq!(list(&["a", "b"]).items(), ["a", "b"]);
    }

    #[test]
    fn rejects_malformed_front_matter() {
        assert_eq!(
            parse_front_matter("# title\n"),
            Err(Error::MissingFrontMatter)
        );
        assert_eq!(
            parse_front_matter("---\nstatus: open\n"),
            Err(Error::UnclosedFrontMatter)
        );
        assert!(matches!(
            parse_front_matter("---\nno colon here\n---\n"),
            Err(Error::MalformedLine { line: 1, .. })
        ));
        assert!(matches!(
            parse_front_matter("---\nstatus: open\n\n---\n"),
            Err(Error::MalformedLine { line: 2, .. })
        ));
        assert_eq!(
            parse_front_matter("---\nstatus: open\nstatus: closed\n---\n"),
            Err(Error::DuplicateKey("status".to_string()))
        );
        assert!(matches!(
            parse_front_matter("---\nissues: [\"x\"\n---\n"),
            Err(Error::MalformedValue { .. })
        ));
        assert!(matches!(
            parse_front_matter("---\ntitle: \"open quote\n---\n"),
            Err(Error::MalformedValue { .. })
        ));
    }

    #[test]
    fn reads_numbered_file_names() {
        assert_eq!(numbered_file("0013-docs-consistency-check.md"), Some(13));
        assert_eq!(numbered_file("0001-a.md"), Some(1));
        assert_eq!(numbered_file("0090-solve-full-core.md"), Some(90));
        assert_eq!(numbered_file("0005-k-offset-2d.md"), Some(5));
        assert_eq!(numbered_file("README.md"), None);
        assert_eq!(numbered_file("013-short.md"), None);
        assert_eq!(numbered_file("00013-long.md"), None);
        assert_eq!(numbered_file("0013-name.txt"), None);
        assert_eq!(numbered_file("0013.md"), None);
        assert_eq!(numbered_file("0013-.md"), None);
    }

    #[test]
    fn requires_kebab_case_names() {
        assert_eq!(numbered_file("0013-Docs-check.md"), None);
        assert_eq!(numbered_file("0013-docs_check.md"), None);
        assert_eq!(numbered_file("0013-docs--check.md"), None);
        assert_eq!(numbered_file("0013--docs.md"), None);
        assert_eq!(numbered_file("0013-docs-.md"), None);
        assert_eq!(numbered_file("0013-docs check.md"), None);
        assert_eq!(numbered_file("0013-文書.md"), None);
    }

    #[test]
    fn finds_adr_references() {
        let text = format!(
            "{} と [{}](../adr/0068-x.md)、({})",
            adr(2),
            adr(68),
            adr(93)
        );
        assert_eq!(adr_refs(&text), [2, 68, 93]);
        // 4桁でないもの、語の途中のもの、書式の説明の `ADR-NNNN` は数えない
        assert!(adr_refs("ADR-NNNN ADR-001 ADR-00012 XADR-0001 ADR-0001x").is_empty());
    }

    #[test]
    fn finds_issue_references() {
        let text = format!(
            "{}と ({})、[\"{}\", \"{}\"]",
            issue(13),
            issue(1),
            issue(2),
            issue(14)
        );
        assert_eq!(issue_refs(&text), [13, 1, 2, 14]);
        // 見出しの記号、`##`、5桁、16進の色、語の途中のものは数えない
        assert!(issue_refs("## 0001 ##0001 #00012 #0000ff a#0001 #001").is_empty());
    }

    fn relative(path: &str) -> Link {
        Link::Relative(path.to_string())
    }

    #[test]
    fn finds_inline_links() {
        let text = "[題](../adr/0068-solve.md)、[計画](../plans/m3.md#節) と \
                    ![図](fig.png)、[外](https://example.com/a.md)、[外](http://x)、\
                    [メール](mailto:a@b)、[同じ文書](#anchor)、[題](a.md \"title\")";
        assert_eq!(
            links(text),
            [
                relative("../adr/0068-solve.md"),
                relative("../plans/m3.md"),
                relative("fig.png"),
                relative("a.md"),
            ]
        );
        assert!(links("[text] (not a link)").is_empty());
    }

    #[test]
    fn reports_absolute_links_and_other_schemes() {
        assert_eq!(
            links("[a](/docs/a.md) [b](ftp://host/b) [c](file:///c.md) [d](C:/d.md)"),
            [
                Link::Absolute("/docs/a.md".to_string()),
                Link::OtherScheme("ftp://host/b".to_string()),
                Link::OtherScheme("file:///c.md".to_string()),
                Link::OtherScheme("C:/d.md".to_string()),
            ]
        );
    }

    #[test]
    fn finds_reference_definitions() {
        assert_eq!(
            links("[label]: ../plans/m3.md#節"),
            [relative("../plans/m3.md")]
        );
        assert_eq!(links("   [x y]: a.md \"title\""), [relative("a.md")]);
        assert_eq!(links("[web]: https://example.com"), Vec::<Link>::new());
        assert_eq!(links("[abs]: /a.md"), [Link::Absolute("/a.md".to_string())]);
        // 脚注と、4つ以上の空白で字下げしたもの（コード）は定義でない
        assert!(links("[^1]: note.md").is_empty());
        assert!(links("    [x]: a.md").is_empty());
    }

    #[test]
    fn ignores_links_in_code_spans() {
        assert!(links("`[x](nowhere.md)` and ``[a](`b`.md)``").is_empty());
        assert_eq!(links("`code` then [a](a.md)"), [relative("a.md")]);
        // 対にならないバッククォートはそのまま
        assert_eq!(links("` [a](a.md)"), [relative("a.md")]);
    }

    #[test]
    fn strips_code_spans() {
        assert_eq!(strip_code_spans("a `b` c"), "a  c");
        assert_eq!(strip_code_spans("a ``b ` c`` d"), "a  d");
        assert_eq!(strip_code_spans("a `b"), "a `b");
        assert_eq!(strip_code_spans("no code"), "no code");
    }

    #[test]
    fn skips_front_matter_and_fences() {
        let text = "---\n\
                    status: open\n\
                    ---\n\
                    one\n\
                    ```markdown\n\
                    hidden\n\
                    ~~~\n\
                    still hidden\n\
                    ```\n\
                    two\n\
                    ~~~~\n\
                    hidden\n\
                    ```\n\
                    ~~~~~\n\
                    three\n";
        assert_eq!(
            markdown_lines(text).unwrap(),
            [(4, "one"), (10, "two"), (15, "three")]
        );
        assert_eq!(
            markdown_lines("a\n---\nb\n").unwrap(),
            [(1, "a"), (2, "---"), (3, "b")]
        );
    }

    #[test]
    fn rejects_unclosed_blocks() {
        assert_eq!(
            markdown_lines("a\n\n```\ncode\n"),
            Err(Error::UnclosedFence { line: 3 })
        );
        assert_eq!(
            markdown_lines("````\n```\n"),
            Err(Error::UnclosedFence { line: 1 })
        );
        assert_eq!(
            markdown_lines("---\nstatus: open\n"),
            Err(Error::UnclosedFrontMatter)
        );
    }
}
