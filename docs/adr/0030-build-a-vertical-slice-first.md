---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
---

# ADR-0030: 最初に3つの段階を最小構成で端から端までつなぐ

## 背景

3つの段階（[ADR-0001](0001-split-into-three-stage-crates.md)）をどの順に作るかを決める必要があった。

## 検討した案

計画には、ほかの順序を比べた記録はない。

## 決定と理由

最初に、3つの段階を最小構成で端から端までつなぐ（縦切り、M0）。段階間のデータ形式の設計ミスに早く気づくためである。
その後で各段階を作り込む。M0 で決めた簡略化は、後のマイルストーンで置き換える前提のものとする。

## 確認

M0 の完了条件（ENDF から作った多群ライブラリで、ピンセルの k∞ と均質な炉心の k_eff を、3段階のコマンドを順に実行して
求められること）。examples/m0/README.md。

## 結果

M0 の簡略化は [ADR-0038](0038-use-zero-kelvin-cross-sections-in-m0.md)〜[ADR-0045](0045-expand-nem-to-fourth-order.md)。
