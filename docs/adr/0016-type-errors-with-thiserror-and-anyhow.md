---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0016: エラーはライブラリ側を thiserror、コマンド側を anyhow とする

## 背景

M2 までのエラーは `String` だった。エラーの種類で分岐できず、原因の連鎖も残らない。

## 検討した案

- `String` のエラー（M2 までの形）
- ライブラリ側を thiserror の型付きのエラー、コマンド側を anyhow（採用）

## 決定と理由

エラーは、ライブラリ側を thiserror の型付きのエラー、コマンド側を anyhow で文脈を付けて表示する。事実上の標準の
クレートを使う（[ADR-0006](0006-use-de-facto-crates.md)）。

## 確認

型の使い分けそのものを確かめるテストはない。エラーの原因の連鎖の表示は `cell_errors_keep_the_cause_in_the_chain_only`、
`material_errors_keep_the_cause_in_the_chain_only`、`toml_errors_are_not_repeated_in_the_chain` で確かめる。内部の検査の
補助関数には、メッセージを `String` で返して型付きのエラーの variant に包むものがある（例: crates/format/src/library.rs の
`check`）。

## 結果

`String` のエラーは M2 の中で移行した。panic しないことは [ADR-0017](0017-forbid-panics-in-production.md)。
