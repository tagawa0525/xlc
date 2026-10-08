---
status: superseded
date: 2026-10-09
supersedes: []
superseded-by: ADR-0094
issues: []
---

# ADR-0093: 格子計算の外周の境界条件は、辺ごとに鏡面反射、真空、アルベドから選ぶ

## 背景

M2（2026-10-06）では、外周の境界条件を鏡面反射と真空（C5G7 のため）とした（docs/plans/m2-assembly.md、d59903e）。
その後、アルベドの境界を加えた。計画には書かれていないが、ブランチ `feat/albedo-boundary` のマージ（e6a3356、
2026-10-08）の本文に理由が書かれている。

## 検討した案

- 鏡面反射と真空だけ（M2 の形）: 一部だけ漏れる問題を表せない
- アルベド α の境界を加え、真空（α = 0）と鏡面反射（α = 1）をその特別な場合とする（採用）

## 決定と理由

外周の境界条件は、辺ごとに鏡面反射、真空、アルベドから選ぶ。鏡面反射は cyclic tracking（[ADR-0060](0060-track-whole-assembly-cyclically.md)）でトラックを
つなぎ、真空は入射を 0 とする。一部だけ漏れる問題には、全反射と真空だけでは足りず、アルベドがその一般の場合（真空 0、
反射 1）になるため（マージ e6a3356 の本文）。

- アルベドの辺は鏡面反射と同じくトラックをつなぎ、出ていく角度中性子束を α（戻る分）と 1 − α（漏れ）に分ける。
  中性子の釣り合いは作り方から保たれる（a5c2b6d）
- 0 ≤ α ≤ 1 を確かめ、α = 1 と 0 は鏡面反射と真空として持つ
- 入力では、辺ごとに `"reflective"`、`"vacuum"`、`{ albedo = α }` のどれかを書く（e64c7c3）
- ADF は鏡面反射の辺だけで作る（それ以外では均質な中性子束が平らにならないため）

計画に、アルベドを使う問題（どの検証の対象に使うか）の記録はない。

## 確認

`lattice_tracks_reflect_at_the_outer_boundary`、`vacuum_sides_have_no_links`、`sides_may_have_an_albedo`、
`albedo_sides_lie_between_vacuum_and_reflection`、`albedo_must_be_a_fraction`、`adf_needs_reflective_sides`。C5G7 との比較
（examples/c5g7/verification.md）。

## 結果

- Dancoff 係数の計算での真空やアルベドの境界の扱いは決めていない（#0008）
- トラッキングは [ADR-0060](0060-track-whole-assembly-cyclically.md)
