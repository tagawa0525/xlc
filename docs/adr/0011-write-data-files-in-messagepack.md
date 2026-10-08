---
status: accepted
date: 2026-10-06
supersedes: [ADR-0009]
superseded-by:
issues: []
---

# ADR-0011: データファイル（プログラムが書いて読むもの）は MessagePack とする

## 背景

最初は自前のテキスト形式とした（[ADR-0009](0009-write-text-intermediate-files.md)）。多群ライブラリや少数群定数は人が読む必要がなく、f64 を誤差なく
受け渡したい。

- 数値を誤差なく受け渡せること
- 検証用のスクリプト（Python）から読めること
- Rust と Python での実績

## 検討した案

- 自前のテキスト形式（[ADR-0009](0009-write-text-intermediate-files.md)）: 10 進表記との変換が入る
- JSON: 読み込み側の実装と設定によっては、最下位のビットがずれる
- CBOR: MessagePack と同じ性質を持つが、Rust と Python での実績が少ない
- HDF5: 大きな表の一部だけを読める。今は要らない
- MessagePack（採用）

## 決定と理由

データファイル（多群ライブラリ、少数群定数など）は MessagePack とする。人が読む必要はなく、検証用のスクリプトから
読めればよい（Python の `msgpack`）。f64 を 8 バイトの IEEE 754 のまま持つので、10 進表記との変換による誤差が
原理的に起きない。CBOR も同じ性質を持つが、Rust と Python での実績が多い MessagePack を選んだ。

serde で型を定義しておけば、形式の変更は読み書きの部分だけで済む。

## 確認

`round_trips_bit_exactly`（書いて読んだ値がビット単位で同じ）、`written_as_a_map_with_the_format_name`。

## 結果

- 大きな表を一部だけ読む必要が出たら（燃焼の分岐テーブルを扱う M4 が目安）、HDF5 への移行を検討する。分岐テーブルは
  作り直しが重いので、この判断は M4 の分岐テーブルの設計と同時に行う
- 読み込みの厳密さは [ADR-0013](0013-read-files-strictly.md)、後方互換は [ADR-0012](0012-drop-backward-compatibility.md)
