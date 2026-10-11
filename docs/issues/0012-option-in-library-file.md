---
status: closed
labels: [format]
opened: 2026-10-08
closed: 2026-10-11
milestone: M3
blocked-by: []
resolution: コミット 41c8491
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
- 2026-10-11: `resolved` に `deserialize_with` を付け、項目の欠けをエラーにした（コミット 41c8491）。型は `Option` のまま、
  なしは nil で書く。形式（map のキーと nil）は変わらず、既存のライブラリ（examples/vera/vera.msgpack）は読める。
  `fewgroup.rs` の `Pins` の `Option` も同じ方式（table.rs の `nil_or_pins`）。table.rs の `AxisFile` の
  `density` / `b10` は軸の名前で必須の項目が決まり、欠けは別のエラーになる（`reading_rejects_wrong_axis_fields`）。table の `pins` の欠けは `states_need_the_pins_key` が拒む。ADR-0013 の「帰結」の食い違いの記述は、
  ADR が確定済みなのでそのまま残す
