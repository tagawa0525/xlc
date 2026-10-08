---
status: open
labels: [format]
opened: 2026-10-08
closed:
milestone: M3
blocked-by: []
resolution:
---

# 多群ライブラリのファイルの型に `Option` がある（[ADR-0013](../adr/0013-read-files-strictly.md) と食い違う）

## 問題

[ADR-0013](../adr/0013-read-files-strictly.md) は、データファイルの型に `Option` を使わないと決めている（serde は欠けた `Option` の項目を黙って `None` に
するため）。しかし多群ライブラリの `NuclideData.resolved` は `Option<Resolved>` である（crates/format/src/library.rs）。
`deny_unknown_fields` は付いているが、欠けた項目には効かない。`resolved` の項目がない核種の map を
`rmpv::ext::from_value` で `NuclideData` に読むと、エラーにならず `resolved = None` になることを確かめた（2026-10-08、
作業用の小さなプログラムで確認。リポジトリのテストはない）。分離共鳴の範囲を持つ核種のファイルからこの項目が欠けても、
黙って「範囲なし」として読まれる。同じファイルの `pointwise_range` は、空の配列で「なし」を表して `Option` を避けている。

関係する ADR: [ADR-0013](../adr/0013-read-files-strictly.md)

## 完了の条件

`resolved` を `Option` でない形（`pointwise_range` と同じく空の配列、または種類を表す enum）にして、項目が欠けた
ファイルをエラーにするテストを加える。ほかのデータファイルの型（少数群定数など）にも `Option` がないことを確かめる。

## 経緯

- 2026-10-08: ADR のレビューで見つかった
