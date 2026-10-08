---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
---

# ADR-0049: 幅のゆらぎの期待値は、閉じた式と適応積分で求め、Hwang の求積表を使わない

## 背景

非分離共鳴の自己遮蔽（[ADR-0048](0048-shield-unresolved-like-unresr.md)）は、幅のゆらぎ（χ² 分布）の期待値を使う。NJOY は Hwang の求積表（10 点の Gauss
求積とその表）でこれを求める。

- 精度（評価済みデータだけで決まり、求積の誤差を持ち込まないこと）
- 計算量

## 検討した案

- NJOY と同じ Hwang の求積表: 典型的な幅の比で、自由度 1 で最大 0.9%、自由度 2 で最大 0.5% の誤差がある
- 10 点の Gauss 求積: 最大 8% の誤差になる
- 閉じた式と、許容誤差を指定した適応積分（採用）

## 決定と理由

NJOY の求積表は使わず、許容誤差を指定した適応積分で求める。中性子幅の期待値は閉じた式で求める（E[1/(x+c)] は自由度
1, 3 では erfc、2, 4 では E1 で書ける。c が大きいところは 1/c の級数）。核分裂幅と競合幅は適応 Gauss-Kronrod 積分で
求める。誤差の値は docs/plans/m1-nuclear-data.md（1fead04、2026-10-06）による。積分は非分離共鳴の格子点でだけ行うので、
計算量は問題にならない。

## 確認

`width_ratios_match_integration`、`reproduces_moments_of_chi_square`、`nested_expectations_converge_in_either_order`、
`three_fluctuating_widths_converge`。NJOY UNRESR との比較では最大 4e-4 ずれ、その原因は UNRESR の求積表による
（examples/m1/README.md、d6b405e、2026-10-06）。

## 結果

NJOY と比べるときは、求積表の誤差の分だけ違う。
