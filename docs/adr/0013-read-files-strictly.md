---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: ["#0012"]
requires: []
---

# ADR-0013: ファイルの読み込みは厳密にし、推測で補わない

## 背景

入力ファイルは手で書き、データファイルは段階をまたいで渡る。誤ったファイルを黙って読むと、結果の誤りの原因を
追えなくなる。

## 検討した案

- 欠けた項目を既定値で補い、未知の項目を無視する
- すぐにエラーにする（採用）

## 決定と理由

読み込みでは、未知の項目、欠けた項目、範囲外の値をすぐにエラーにする（serde の `deny_unknown_fields`）。推測で
補ったり、既定値で黙って続けたりしない。データファイルは形式の名前（種類）を持ち、別の種類のファイルを渡されたら
エラーにする。

- MessagePack はフィールド名付きの map 形式で書き、読むときも map 以外を拒否する。配列形式ではフィールド名がなく、
  未知の項目を検出できず、順番の変わったファイルを黙って誤読するため
- データファイルの型に `Option` を使わない（serde は欠けた `Option` の項目を黙って `None` にする）。入力ファイルでは、
  省略に意味がある項目（出力しない、帯をセルの外周までとする、など）だけ省略を許す

## 確認

`unknown_and_missing_fields_are_rejected`、`array_form_is_rejected`、`nested_structs_in_array_form_are_rejected`、
`other_formats_are_rejected`、`trailing_bytes_are_rejected`、`non_finite_values_are_rejected`。欠けた `Option` の項目を
拒むことは確かめておらず、実際に拒まない（#0012）。

## 帰結

- 形式を変えたら古いファイルは読めなくなる（[ADR-0012](0012-drop-backward-compatibility.md)）
- 多群ライブラリの `NuclideData.resolved` は `Option` で、項目が欠けても `None` として読まれる。この ADR と
  食い違っている（#0012）
