---
status: open
labels: [lattice, shielding, verification]
opened: 2026-10-08
closed:
milestone: M2
blocked-by: []
resolution:
---

# 集合体全体を UFG の MOC で解く参照の道具がない

## 問題

[ADR-0068](../adr/0068-solve-pinwise-slowing-down.md) は、ピンごとの等価なセルの誤差を、集合体全体を UFG（ポイントワイズ）の
MOC で直接解く参照の道具で、案内管や集合体の端のそばのピンで測るとしていた（M2 の進め方の 6）。この道具は作って
いない。等価なセルの誤差は、OpenMC との比較（examples/vera/verification.md）で見ているだけである。

関係する ADR: なし

## 完了の条件

参照の道具を作り、VERA 2B などで等価なセルの誤差をピンごとに測るか、OpenMC との比較で足りるとして ADR に
記録し、M2 の計画の項目を閉じる。

## 経緯

- 2026-10-08: 計画に進み具合を書くときに、M2 の項目の残りとして起こした
