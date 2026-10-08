---
status: accepted
date: 2026-10-08
supersedes: []
superseded-by:
issues: ["#0008"]
---

# ADR-0093: 格子計算の外周の辺に、アルベドの境界条件を選べるようにする

## 背景

格子計算の外周の境界条件は、鏡面反射と真空だった（[ADR-0060](0060-track-whole-assembly-cyclically.md)）。計画には書かれていないが、アルベドの境界を加えた
ブランチ `feat/albedo-boundary` のマージ（21e52b5、2026-10-08）の本文に理由が書かれている。

## 検討した案

- 鏡面反射と真空だけ（[ADR-0060](0060-track-whole-assembly-cyclically.md)）: 一部だけ漏れる問題を表せない
- アルベド α の境界を加え、真空（α = 0）と鏡面反射（α = 1）をその特別な場合とする（採用）

## 決定と理由

一部だけ漏れる問題には、全反射と真空だけでは足りず、アルベドがその一般の場合（真空 0、反射 1）になるため（マージ
21e52b5 の本文）。

- アルベドの辺は鏡面反射と同じくトラックをつなぎ、出ていく角度中性子束を α（戻る分）と 1 − α（漏れ）に分ける。
  中性子収支は作り方から保たれる（45e0894）
- 0 ≤ α ≤ 1 を確かめ、α = 1 と 0 は鏡面反射と真空として持つ
- 入力では、辺ごとに `"reflective"`、`"vacuum"`、`{ albedo = α }` のどれかを書く（52f0c27）
- ADF は鏡面反射の辺だけで作る（それ以外では均質な中性子束が平らにならないため）

計画に、アルベドを使う問題（どの検証の対象に使うか）の記録はない。

## 確認

`sides_may_have_an_albedo`、`albedo_sides_lie_between_vacuum_and_reflection`、`albedo_must_be_a_fraction`、
`adf_needs_reflective_sides`。

## 結果

- Dancoff 係数の計算での真空やアルベドの境界の扱いは決めていない（#0008）
- トラッキングと鏡面反射・真空の境界は [ADR-0060](0060-track-whole-assembly-cyclically.md)
