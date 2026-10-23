---
status: proposed
date: 2026-10-23
requires: [ADR-0097, ADR-0187, ADR-0190]
supersedes: []
superseded-by:
amends: [ADR-0187]
issues: []
---

# ADR-0192: core は燃焼度の軸を、3 点の放物線の傾きの 3 次 Hermite で内挿し、燃焼度 0 からの区間だけ線形にする

## 背景

[ADR-0187](0187-tabulate-burnup-branches-from-history.md) は、core が燃焼度の軸を隣の 2 点で線形に内挿する（案 I1）とした。[ADR-0190](0190-advance-node-burnups-by-midpoint-steps.md) の確認（2A の集合体 1 体）では、
履歴の 40 点から選んだ 17 点のテーブルで、点の間の k∞ の誤差が −32.8〜+179.7 pcm になった。ホウ素に直すと 15〜20 ppm で、Problem 9 の
サイクルの追跡に使えない。

（下書き。測定のあとで書く）

## 検討した案

## 決定と理由

これは ADR-0187 の「燃焼度の内挿（案 I1）」の部分を置き換える。ほかの部分は変わらない。

## 確認

## 帰結
