---
type: Issue
status: open
labels: [xs]
opened: 2026-11-20
closed:
milestone:
blocked-by: []
resolution:
---

# MF6 の LAW=1 と MF5 の表の出射分布の規格化を確かめていない

## 問題

- MF6 の LAW=7 は、入射エネルギーごとの ∫∫ f dμ dE' が 1 から 1e-5 を超えてずれたら止める（ADR-0271、#0323）。LAW=1 と MF5 の表の出射分布
  （`OutgoingDistributions`）には同じ検査が無く、和が 1 に足りない表は中性子を黙って失う（どの評価で起きるかは調べていない）

関係する ADR: [ADR-0271](../adr/0271-group-mf6-law6-and-law7-by-endf-102.md)

## 完了の条件

- [ ] 4 つの版の表の規格化のずれを測り、検査を足すかを決める

## 経緯

- 2026-11-20: #0311 のレビューで、LAW=7 だけに規格化の検査があると指摘された（#0323 の範囲の外として起こした）
