---
status: accepted
date: 2026-10-11
requires: [ADR-0013]
supersedes: []
superseded-by:
issues: ["#0012"]
---

# ADR-0129: 「なし」を表す項目は Option にし、欠けを拒んで nil を明示させる

## 背景

[ADR-0013](0013-read-files-strictly.md) は、欠けた項目を黙って受けないために、データファイルの型に `Option` を使わない
と決めた（serde は欠けた `Option` の項目を黙って `None` にする）。実装は食い違っていた。多群ライブラリの
`NuclideData.resolved`（分離共鳴の範囲。crates/format/src/library.rs）と、少数群定数の表の `pins`
（crates/format/src/table.rs）は `Option` で、#0012 で `deserialize_with`（`nil_or_resolved`、`nil_or_pins`）を付けて、
項目が欠けたらエラー、nil を明示したときだけ `None` とした。

判断の軸:

- 欠けた項目を黙って受けない（ADR-0013 の目的）
- ファイルの形式（MessagePack の map のキーと nil）を変えない。既存のライブラリ（examples/vera/vera.msgpack）を
  読めるままにする

## 検討した案

- 空の配列で「なし」を表す: `pointwise_range` はこの形である。しかし `resolved` と `pins` は配列でない構造体で、
  配列にすると map のキーと値の形が変わる。「要素が 0 個か 1 個」という制約が型に出ず、2 個以上のファイルを別に
  拒む必要がある
- 種類を表す enum にする: 「なし」と「あり」を型で区別できるが、MessagePack での表現（タグつきの map など）が
  変わり、既存のファイルを作り直す必要がある。項目が 2 通りしかないものに対して、読み込みの規則が増える
- `Option` + `deserialize_with` で nil を明示させる（採用）: `deserialize_with` を付けると、serde は欠けた項目を
  `None` にせず、欠けをエラーにする。書く側は `None` を nil で書く。形式は変わらない

空の配列と enum を採らなかった理由の経緯は #0012 にある。

## 決定と理由

「なし」を表す項目は `Option` にしてよい。ただし必ず `deserialize_with` で `Option::deserialize` を呼ぶ関数を付けて欠けを
拒み、ファイルには nil を明示して書く。serde の `default` は使わない（欠けを `None` にしてしまうため）。

これは ADR-0013 の「データファイルの型に `Option` を使わない」の部分を置き換える。ADR-0013 のほかの部分（未知の項目、
欠けた項目、範囲外の値の拒否など）は変わらない。ADR 全体の置き換えではないので、ADR-0013 のフロントマターは変えない。
欠けた項目を黙って受けないという目的は、この方式でも満たされる。

## 帰結

- ADR-0013 の「確認」と「帰結」にある、`Option` の項目の欠けを拒まず、`resolved` が ADR と食い違うという記述は、
  #0012 の前の状態であり、今は当てはまらない
- 新しい `Option` の項目を足すときは、欠けを拒むテストを同じ変更に入れる

## 確認

- `states_need_the_pins_key`（crates/format/src/table/tests.rs）: `pins` のキーがない表を拒む
- `nuclide_without_the_resolved_field_is_rejected`（crates/format/src/library.rs）: `resolved` のキーがない核種を拒む
