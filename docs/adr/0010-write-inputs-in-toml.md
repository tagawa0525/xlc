---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
---

# ADR-0010: 入力ファイル（人が書くもの）は TOML とする

## 背景

段階間のファイルと入力ファイルは、すべてのクレートと検証用のツールに広がるので、後から変えにくい。M2 の途中で
ファイルの約束事を決めた。

- 手で書きやすく、コメントが書けること
- Rust の serde 対応が安定していること

## 検討した案

- YAML: Rust での serde 対応の保守が途絶えている
- JSON: コメントが書けない
- XML: 手書きに向かない
- TOML（採用）

## 決定と理由

入力ファイル（人が書くもの）は TOML とする。コメントが書け、手書きの誤りに強く、Rust の serde 対応が安定しているため。

## 確認

入力の読み込みのテスト（`parses_lattice`、`parses_pin_cell`、`rejects_bad_input` など）。

## 結果

- 既存の自前のテキスト形式の入力は M2 の中で移行した
- 読み込みの厳密さは [ADR-0013](0013-read-files-strictly.md)
