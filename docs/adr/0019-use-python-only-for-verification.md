---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
---

# ADR-0019: Python は検証用のツールにだけ使い、環境は flake.nix で固定する

## 背景

OpenMC の入力の作成や比較のスクリプトには Python が便利である。本体の計算とは分けたい。

## 検討した案

計画には、ほかの案を比べた記録はない。

## 決定と理由

Python は検証用のツール（OpenMC の入力の作成、比較のスクリプト）にだけ使い、本体の計算には使わない。
開発と検証の環境は `flake.nix` の devShell で固定する。

## 確認

約束事で、確かめるテストはない。

## 結果

OpenMC は検証にだけ使う外部ツールとする（[ADR-0031](0031-verify-against-openmc-with-same-data.md)）。
