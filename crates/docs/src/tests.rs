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
fn keeps_hashes_in_values() {
    let text = format!(
        "---\n\
             status: accepted # note\n\
             resolution: wontfix: see {i1}\n\
             issues: [\"{i1}\", \"{i14}\"]\n\
             ---\n",
        i1 = issue(1),
        i14 = issue(14),
    );
    let fm = parse_front_matter(&text).unwrap();
    assert_eq!(fm.get("status"), Some(&scalar("accepted # note")));
    assert_eq!(
        fm.get("resolution"),
        Some(&scalar(&format!("wontfix: see {}", issue(1))))
    );
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
fn skips_fences_but_not_front_matter() {
    // フロントマターは読み分けない（番号つきの文書のフロントマターも、参照を探す行に入れる）
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
        [
            (1, "---"),
            (2, "status: open"),
            (3, "---"),
            (4, "one"),
            (10, "two"),
            (15, "three")
        ]
    );
    assert_eq!(
        markdown_lines("---\nstatus: open\n").unwrap(),
        [(1, "---"), (2, "status: open")]
    );
}

#[test]
fn follows_commonmark_fences() {
    // 囲いの記号のあとにバッククォートを含む行は、コードブロックを開かない
    assert_eq!(
        markdown_lines("```cargo test``` で走る\nx\n").unwrap(),
        [(1, "```cargo test``` で走る"), (2, "x")]
    );
    // 4つ以上の空白で字下げした囲いは、囲いでない
    assert_eq!(
        markdown_lines("    ```\nx\n").unwrap(),
        [(1, "    ```"), (2, "x")]
    );
    assert!(markdown_lines("   ```\nx\n   ```\n").unwrap().is_empty());
    // 閉じる囲いは、同じ記号が開いたときと同じ数以上並び、ほかに何もない行
    assert_eq!(
        markdown_lines("```\n``` not closing\n```\ny\n").unwrap(),
        [(4, "y")]
    );
    assert_eq!(
        markdown_lines("```rust\nx\n    ```\n```\ny\n").unwrap(),
        [(5, "y")]
    );
    // `~~~` の情報の文字列にはバッククォートを書ける
    assert_eq!(
        markdown_lines("~~~ info with ` tick\nx\n~~~\ny\n").unwrap(),
        [(4, "y")]
    );
}

#[test]
fn rejects_unclosed_fences() {
    assert_eq!(
        markdown_lines("a\n\n```\ncode\n"),
        Err(Error::UnclosedFence { line: 3 })
    );
    assert_eq!(
        markdown_lines("````\n```\n"),
        Err(Error::UnclosedFence { line: 1 })
    );
}
