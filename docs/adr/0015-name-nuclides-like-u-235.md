---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0015: 核種の名前は `U-235`、`Am-242m` の形とする

## 背景

核種の名前は、入力ファイル、データファイル、検証用のツールのすべてに現れる。

## 検討した案

- IUPAC, IUPAP の上付きの記号 ²³⁵U や名前 uranium-235: 規格が定めるのはこの2つだけで、上付きを使わない書き方の
  規格はない
- GNDS や OpenMC の `U235`, `Am242_m1`
- 元素記号と質量数をハイフンでつなぐ `U-235`（採用）

## 決定と理由

核種の名前は、元素記号と質量数をハイフンでつなぐ `U-235`, `H-1` の形とし、元素記号は先頭だけ大文字にする。
準安定状態は `Am-242m` と書き、2番目以降は `Hf-178m2` のように番号を付ける（`m` は `m1` と同じ。番号は GNDS の
metastableIndex に合わせる）。この形は名前を記号で略しただけに見えるので、原子力分野の外でも読める。

## 確認

`nuclide_names_follow_the_convention`、`library_nuclide_names_are_checked`。

## 結果

GNDS や OpenMC の `U235`, `Am242_m1` とは、受け渡すときに文字列の変換で対応させる。
