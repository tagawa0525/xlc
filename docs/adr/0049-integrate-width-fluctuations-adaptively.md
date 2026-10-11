---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0049: 幅のゆらぎの期待値は、閉じた式と適応積分で求め、NJOY の Hwang の求積を使わない

## 背景

非分離共鳴の自己遮蔽（[ADR-0048](0048-shield-unresolved-like-unresr.md)）は、幅のゆらぎ（χ² 分布）の期待値を使う。NJOY の UNRESR（と RECONR の `gnrl`）は、
これを R. N. Hwang の求積で求める。自由度 1〜4 のそれぞれに、χ² 分布に合わせた 10 点の節点と重みの表を持つもので
（`references/NJOY2016/src/unresr.f90` の `qp` と `qw`。冒頭の説明に「The accurate Hwang quadrature set is used for the
integrals」とある）、ふつうの Gauss 求積の点ではない。

- 精度（評価済みデータだけで決まり、求積の誤差を持ち込まないこと）
- 計算量

## 検討した案

誤差の値は、計画に記録したもの（docs/plans/m1-nuclear-data.md、d8e8098、2026-10-06。コミットの本文に同じ値がある）。
どの幅の比で測ったかの詳細は記録にない。

- NJOY と同じ Hwang の求積の表: 典型的な幅の比で、自由度 1 で最大 0.9%、自由度 2 で最大 0.5% の誤差がある
- ふつうの 10 点の Gauss 求積（Hwang の表とは別の、χ² 分布に合わせていない求積）: 最大 8% の誤差になる
- 閉じた式と、許容誤差を指定した適応積分（採用）

## 決定と理由

NJOY の求積の表は使わず、許容誤差を指定した適応積分で求める（2026-10-06）。中性子幅の期待値は閉じた式で求める
（E[1/(x+c)] は自由度 1, 3 では erfc、2, 4 では E1 で書ける。c が大きいところは 1/c の級数。260245f、2026-10-06）。
核分裂幅と競合幅は適応 Gauss-Kronrod 積分で求める。積分は非分離共鳴の格子点でだけ行うので、計算量は問題にならない
（U-236 で競合幅もゆらぐところでは重くなる。examples/m1/verification.md の「計算時間について」）。

## 確認

`width_ratios_match_integration`、`reproduces_moments_of_chi_square`、`nested_expectations_converge_in_either_order`、
`three_fluctuating_widths_converge`。NJOY UNRESR との比較では最大 4e-4 ずれ、その原因は UNRESR の求積表による
（examples/m1/verification.md、46873af、2026-10-06）。上の 0.9%、0.5%、8% を出したスクリプトやテストはリポジトリにない。
NJOY の説明書は手元にないので、Hwang の求積の中身は NJOY2016 のソースで確かめた。

## 結果

NJOY と比べるときは、求積表の誤差の分だけ違う。
