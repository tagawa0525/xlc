//! 文書（docs/adr、docs/issues ほか）の整合性を確かめるための読み取り（#0013）。
//!
//! 検査そのものは tests/consistency.rs にあり、ワークスペースの `cargo test` で走る。ここには、Markdown のフロントマター、
//! 番号つきのファイル名、`ADR-NNNN` と `#NNNN` の参照、相対リンクを読む小さな関数を置く。

// テストでは assert! 系を使う（本番のコードでは clippy.toml で禁止している）
#![cfg_attr(test, allow(clippy::disallowed_macros))]

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
pub fn parse_front_matter(_text: &str) -> Result<FrontMatter, Error> {
    Ok(FrontMatter::default())
}

/// `NNNN-name.md` の形のファイル名から番号を取り出す。形が違えば `None`。
pub fn numbered_file(_name: &str) -> Option<u32> {
    None
}

/// `ADR-NNNN`（4桁）の参照の番号を、現れた順に返す。
pub fn adr_refs(_text: &str) -> Vec<u32> {
    Vec::new()
}

/// issue の参照 `#NNNN`（ちょうど4桁で、前に `#` も英数字もない）の番号を、現れた順に返す。
pub fn issue_refs(_text: &str) -> Vec<u32> {
    Vec::new()
}

/// Markdown の相対リンク `[text](path)` の行き先を、現れた順に返す。`http:`、`https:`、`mailto:` と、
/// 同じ文書の中の `#anchor` だけのものは除き、パスの `#anchor` は外す。
pub fn relative_links(_text: &str) -> Vec<String> {
    Vec::new()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scalar(s: &str) -> Value {
        Value::Scalar(s.to_string())
    }

    fn list(items: &[&str]) -> Value {
        Value::List(items.iter().map(|s| s.to_string()).collect())
    }

    #[test]
    fn reads_front_matter_values() {
        let text = "---\n\
                    status: superseded\n\
                    date: 2026-10-08\n\
                    supersedes: [ADR-0002, ADR-0003]\n\
                    superseded-by: ADR-0007\n\
                    issues: [\"#0001\"]\n\
                    closed:\n\
                    title: \"quoted text\"\n\
                    ---\n\
                    \n\
                    # ADR-0001: title\n";
        let fm = parse_front_matter(text).unwrap();
        assert_eq!(fm.get("status"), Some(&scalar("superseded")));
        assert_eq!(fm.get("date"), Some(&scalar("2026-10-08")));
        assert_eq!(fm.get("supersedes"), Some(&list(&["ADR-0002", "ADR-0003"])));
        assert_eq!(fm.get("superseded-by"), Some(&scalar("ADR-0007")));
        assert_eq!(fm.get("issues"), Some(&list(&["#0001"])));
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
        let text = "---\n\
                    status: accepted      # proposed | accepted\n\
                    supersedes: []        # 例: [ADR-0002]\n\
                    superseded-by:        # 例: ADR-0007\n\
                    issues: [\"#0001\", \"#0014\"] # 例\n\
                    ---\n";
        let fm = parse_front_matter(text).unwrap();
        assert_eq!(fm.get("status"), Some(&scalar("accepted")));
        assert_eq!(fm.get("supersedes"), Some(&list(&[])));
        assert_eq!(fm.get("superseded-by"), Some(&Value::Empty));
        assert_eq!(fm.get("issues"), Some(&list(&["#0001", "#0014"])));
    }

    #[test]
    fn value_items_flatten_scalars_and_lists() {
        assert!(Value::Empty.items().is_empty());
        assert_eq!(scalar("ADR-0007").items(), ["ADR-0007"]);
        assert_eq!(
            list(&["ADR-0068", "ADR-0074"]).items(),
            ["ADR-0068", "ADR-0074"]
        );
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
            parse_front_matter("---\nissues: [\"#0001\"\n---\n"),
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
        assert_eq!(numbered_file("README.md"), None);
        assert_eq!(numbered_file("013-short.md"), None);
        assert_eq!(numbered_file("00013-long.md"), None);
        assert_eq!(numbered_file("0013-name.txt"), None);
        assert_eq!(numbered_file("0013.md"), None);
        assert_eq!(numbered_file("0013-.md"), None);
    }

    #[test]
    fn finds_adr_references() {
        assert_eq!(
            adr_refs("ADR-0002 と [ADR-0068](../adr/0068-x.md)、(ADR-0093)"),
            [2, 68, 93]
        );
        // 4桁でないもの、語の途中のもの、書式の説明の `ADR-NNNN` は数えない
        assert!(adr_refs("ADR-NNNN ADR-001 ADR-00012 XADR-0001 ADR-0001x").is_empty());
    }

    #[test]
    fn finds_issue_references() {
        assert_eq!(
            issue_refs("#0013 と (#0001)、[\"#0002\", \"#0014\"]"),
            [13, 1, 2, 14]
        );
        // 見出しの記号、`##`、5桁、16進の色、語の途中のものは数えない
        assert!(issue_refs("## 0001 ##0001 #00012 #0000ff a#0001 #001").is_empty());
    }

    #[test]
    fn finds_relative_links() {
        let text = "[ADR-0068](../adr/0068-solve.md)、[計画](../plans/m3.md#節) と \
                    ![図](fig.png)、[外](https://example.com/a.md)、[外](http://x)、\
                    [メール](mailto:a@b)、[同じ文書](#anchor)、[題](a.md \"title\")";
        assert_eq!(
            relative_links(text),
            ["../adr/0068-solve.md", "../plans/m3.md", "fig.png", "a.md"]
        );
        assert!(relative_links("[text] (not a link) and `[x]`").is_empty());
    }
}
