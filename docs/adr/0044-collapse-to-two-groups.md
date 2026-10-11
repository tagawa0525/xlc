---
status: superseded
date: 2026-10-09
supersedes: []
superseded-by: ADR-0098
issues: []
requires: []
---

# ADR-0044: 少数群は境界 0.625 eV の2群とし、拡散係数は D = 1/(3Σt) とする

## 背景

格子計算から炉心計算に渡す少数群定数の群と、拡散係数の作り方を決める必要があった。

## 検討した案

計画には、ほかの案を比べた記録はない。

## 決定と理由

2群（境界 0.625 eV）に縮約し、拡散係数は D = 1/(3Σt) とする。

## 確認

`two_group_constants_reproduce_k_inf`、`one_group_collapse_preserves_k` は、縮約が k∞ を保つことを確かめる。
D = 1/(3Σt) の選び方そのものを確かめる比較はない。

## 帰結

なし。
