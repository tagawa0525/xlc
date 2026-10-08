---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
---

# ADR-0039: 群定数化の重み関数は NJOY の iwt=4 と同じものにする

## 背景

群構造は WIMS 69 群とした（[ADR-0022](0022-use-wims-69-groups.md)）。群平均の重み関数を決める必要があった。

- NJOY の結果と同じ条件で比べられること

## 検討した案

計画には、ほかの重み関数を比べた記録はない。

## 決定と理由

重み関数は NJOY の iwt=4 と同じもの（Maxwell 分布 + 1/E + 核分裂スペクトル）にする。同じ条件の NJOY（RECONR, GROUPR）
の結果と比べて検証できるようにするためである。

## 確認

`weight_only_in_maxwellian_region`、`weight_only_in_one_over_e_region`、`group_spanning_weight_breakpoint`。
NJOY GROUPR との比較（`crates/xs/examples/compare_njoy_gendf.rs`、examples/m1/README.md）。

## 結果

分離共鳴の範囲の σ0 ごとの重みは、のちに減速計算で作る（[ADR-0072](0072-weight-shielded-data-by-flux-calculator.md)）。その範囲の外と熱中性子の領域はこの重み関数のまま。
