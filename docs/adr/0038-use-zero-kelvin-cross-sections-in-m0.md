---
status: superseded
date: 2026-10-06
supersedes: []
superseded-by: [ADR-0046, ADR-0048, ADR-0052]
issues: []
requires: []
---

# ADR-0038: M0 の断面積は 0 K・無限希釈とし、熱中性子の散乱だけ 293.6 K の自由ガスで扱う

## 背景

M0 は3段階を最小構成でつなぐ縦切りで（[ADR-0030](0030-build-a-vertical-slice-first.md)）、物理の精度は M1 以降で上げる。温度と共鳴の扱いをどこまで
簡略にできるかを決める必要があった。

## 検討した案

- すべて 0 K とする: 標的の熱運動を無視すると上方散乱がなくなり、熱スペクトルが成り立たない
- 0 K とし、熱中性子の散乱だけ自由ガスで扱う（採用）

## 決定と理由

- 分離共鳴領域は Reich-Moore（LRF=3）の共鳴パラメータから 0 K の断面積を再構成する。非分離共鳴領域は、無限希釈の
  平均断面積が MF3 に入っている（LSSF=1）のでそれを使う
- ドップラー拡がりは M1 で扱い、M0 の断面積は 0 K とする。熱中性子の散乱だけは、標的の熱運動を無視すると上方散乱が
  なくなり熱スペクトルが成り立たないため、293.6 K の自由ガス模型で扱う。熱中性子散乱則 S(α,β) は M1 で扱う

## 確認

M0 の例で、同じ条件の NJOY（RECONR, GROUPR）の結果と比べた（examples/m0/verification.md）。

## 結果

M1 で、ドップラー拡がり（[ADR-0046](0046-broaden-with-sigma1.md)）、非分離共鳴の自己遮蔽（[ADR-0048](0048-shield-unresolved-like-unresr.md)）、H in H2O の S(α,β)（[ADR-0052](0052-scatter-h-in-h2o-by-sab.md)）に
置き換えた。4 eV 未満の弾性散乱を自由ガスで扱うことは、S(α,β) を持たない核種で残る（[ADR-0081](0081-use-mf4-angular-distributions.md)）。
