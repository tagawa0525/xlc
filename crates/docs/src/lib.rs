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
/// 各行は `key: value` で、値は空、裸の文字列、引用符で囲んだ文字列、`[a, "b"]` の一覧のどれか。`#` は値の
/// 一部として読む（コメントは書かない）。空行、重複したキー、閉じていない引用符や角括弧はエラーにする。
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
        let (key, value) = line.trim_end().split_once(':').ok_or_else(malformed)?;
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
    if number.len() != 4 || !number.bytes().all(|b| b.is_ascii_digit()) || !is_kebab_case(rest) {
        return None;
    }
    number.parse().ok()
}

/// 英小文字と数字の語を `-` でつないだもの（`[a-z0-9]+(-[a-z0-9]+)*`）か。
fn is_kebab_case(name: &str) -> bool {
    name.split('-').all(|word| {
        !word.is_empty()
            && word
                .bytes()
                .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit())
    })
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
pub fn links(line: &str) -> Vec<Link> {
    let line = strip_code_spans(line);
    let mut found: Vec<Link> = reference_definition(&line).into_iter().collect();
    found.extend(line.match_indices("](").filter_map(|(i, open)| {
        let rest = &line[i + open.len()..];
        // `[text](target "title")` の題は除く
        rest[..rest.find(')')?]
            .split_whitespace()
            .next()
            .and_then(classify_link)
    }));
    found
}

/// 参照の定義 `[label]: target` の行き先。脚注 `[^label]:` と、4つ以上の空白で字下げした行（コード）は除く。
fn reference_definition(line: &str) -> Option<Link> {
    let indent = line.len() - line.trim_start_matches(' ').len();
    if indent > 3 {
        return None;
    }
    let rest = line[indent..].strip_prefix('[')?;
    let (label, target) = rest.split_once("]:")?;
    if label.is_empty() || label.starts_with('^') || label.contains(['[', ']']) {
        return None;
    }
    target.split_whitespace().next().and_then(classify_link)
}

/// リンク先を分ける。外部の URL と、`#anchor` だけのものは `None`。
fn classify_link(target: &str) -> Option<Link> {
    const EXTERNAL: [&str; 4] = ["http:", "https:", "mailto:", "#"];
    if EXTERNAL.iter().any(|p| target.starts_with(p)) {
        return None;
    }
    if target.starts_with('/') {
        return Some(Link::Absolute(target.to_string()));
    }
    if has_scheme(target) {
        return Some(Link::OtherScheme(target.to_string()));
    }
    let path = target.split('#').next().unwrap_or("");
    (!path.is_empty()).then(|| Link::Relative(path.to_string()))
}

/// `scheme:` で始まるか（英字のあとに英数字か `+ - .` が続き、`:` で終わる）。
fn has_scheme(target: &str) -> bool {
    target.split_once(':').is_some_and(|(scheme, _)| {
        scheme.starts_with(|c: char| c.is_ascii_alphabetic())
            && scheme
                .chars()
                .all(|c| c.is_ascii_alphanumeric() || matches!(c, '+' | '-' | '.'))
    })
}

/// 1行からコードスパン（対になったバッククォートで囲んだ部分）を除く。
pub fn strip_code_spans(line: &str) -> String {
    let mut out = String::new();
    let mut rest = line;
    while let Some(start) = rest.find('`') {
        out.push_str(&rest[..start]);
        let ticks = backtick_run(&rest[start..]);
        let body = &rest[start + ticks..];
        match closing_run(body, ticks) {
            Some(end) => rest = &body[end + ticks..],
            None => {
                // 対にならないバッククォートは文字のまま
                out.push_str(&rest[start..start + ticks]);
                rest = body;
            }
        }
    }
    out.push_str(rest);
    out
}

/// 先頭に続くバッククォートの数。
fn backtick_run(text: &str) -> usize {
    text.len() - text.trim_start_matches('`').len()
}

/// `text` の中で、ちょうど `ticks` 個続くバッククォートの位置。
fn closing_run(text: &str, ticks: usize) -> Option<usize> {
    let mut from = 0;
    while let Some(offset) = text[from..].find('`') {
        let at = from + offset;
        let run = backtick_run(&text[at..]);
        if run == ticks {
            return Some(at);
        }
        from = at + run;
    }
    None
}

/// 4つ以上の空白で字下げしていない行の、字下げを除いた残り。
fn unindented(line: &str) -> Option<&str> {
    let rest = line.trim_start_matches(' ');
    (line.len() - rest.len() <= 3).then_some(rest)
}

/// 先頭に続く `c` の数。
fn run_of(text: &str, c: char) -> usize {
    text.len() - text.trim_start_matches(c).len()
}

/// コードブロックを開く囲い（CommonMark）の記号と数。3つ以上の ```` ` ```` か `~` で、字下げは3つまで。
/// ```` ` ```` の囲いの情報の文字列にはバッククォートを書けない（```` ```cargo test``` ```` は囲いでない）。
fn opening_fence(line: &str) -> Option<(char, usize)> {
    let rest = unindented(line)?;
    ['`', '~'].into_iter().find_map(|c| {
        let n = run_of(rest, c);
        (n >= 3 && !(c == '`' && rest[n..].contains('`'))).then_some((c, n))
    })
}

/// `c` が `n` 個以上並び、ほかに何もない、字下げが3つまでの行か（閉じる囲い）。
fn closes_fence(line: &str, c: char, n: usize) -> bool {
    unindented(line).is_some_and(|rest| {
        let run = run_of(rest, c);
        run >= n && rest[run..].trim().is_empty()
    })
}

/// Markdown の、コードブロック（CommonMark の ```` ``` ```` か `~~~` の囲い）の外の行を、行番号（1 から）とともに
/// 返す。フロントマターは読み分けない（ほかの行と同じに返す）。閉じていないコードブロックはエラーにする。
pub fn markdown_lines(text: &str) -> Result<Vec<(usize, &str)>, Error> {
    let mut out = Vec::new();
    // 開いているコードブロックの記号、数、開いた行
    let mut fence: Option<(char, usize, usize)> = None;
    for (i, line) in text.lines().enumerate() {
        let number = i + 1;
        match fence {
            None => match opening_fence(line) {
                Some((c, n)) => fence = Some((c, n, number)),
                None => out.push((number, line)),
            },
            Some((c, n, _)) => {
                if closes_fence(line, c, n) {
                    fence = None;
                }
            }
        }
    }
    match fence {
        Some((_, _, line)) => Err(Error::UnclosedFence { line }),
        None => Ok(out),
    }
}

#[cfg(test)]
mod tests;
