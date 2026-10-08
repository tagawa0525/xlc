---
status: accepted
date: 2026-10-07
supersedes: []
superseded-by:
issues: []
---

# ADR-0058: 非分離共鳴の LRF=1 は LRF=2 に読み替え、l=3 は剛体球の式で計算して warn する

## 背景

Fe-58 の非分離共鳴は LRF=1 で、ENDF が非分離共鳴に定めない l=3 の系列を持つ。NJOY（UNRESR, RECONR, PURR）は l≥2 を
l=2 の透過因子と位相シフトで計算してしまう（UNRESR では l=3 を l=2 にすると7桁一致した）。

## 検討した案

- NJOY と同じく l=2 の式で計算する: 物理として誤り
- 分離共鳴と同じ l=3 の剛体球の式を使い、仕様外であることを warn で知らせる（採用）

## 決定と理由

LRF=1 は、LRF=2 と同じ形（パラメータをエネルギーに線形補間）に読み替えて、同じ自己遮蔽の計算（[ADR-0048](0048-shield-unresolved-like-unresr.md)）を使う。
l=3 は分離共鳴と同じ l=3 の剛体球の式を使い、仕様外であることを warn で知らせる。NJOY の値とは 350 keV で全断面積が
1.5%、捕獲が 42% 違うが（docs/plans/m2-assembly.md、a3d4433、2026-10-07）、Zircaloy-4 の中の Fe-58 はごくわずかなので
影響はない。

## 確認

`reads_unresolved_lrf1_with_fission_widths`、`reads_unresolved_lrf1_without_fission_widths`。NJOY UNRESR で l=3 を l=2 に
すると7桁一致することで、差の原因を確かめた。

## 結果

NJOY と比べるときは、Fe-58 の非分離共鳴で上の差が残る。
