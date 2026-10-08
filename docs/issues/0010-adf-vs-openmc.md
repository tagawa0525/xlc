---
status: open
labels: [lattice, verification]
opened: 2026-10-08
closed:
milestone: M3
blocked-by: []
resolution:
---

# ADF を OpenMC と比べていない

## 問題

集合体の面ごとの不連続因子（ADF）は、反射境界の場合だけ作っているが、OpenMC などの参照と比べていない。

関係する ADR: [ADR-0026](../adr/0026-output-two-group-tables-from-lattice.md)、[ADR-0084](../adr/0084-make-reflector-constants-with-fuel.md)

## 完了の条件

参照の作り方（OpenMC で面の中性子束と均質化した解）を決め、2B などで比べる。非反射境界の ADF（GET）を作るときに
あわせて確かめる。

## 経緯

- 2026-10-08: メモから移した
