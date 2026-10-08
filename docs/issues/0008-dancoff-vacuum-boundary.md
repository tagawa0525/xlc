---
status: open
labels: [lattice, shielding]
opened: 2026-10-08
closed:
milestone: M3
blocked-by: []
resolution:
---

# Dancoff 係数の計算での真空境界（C6）

## 問題

Dancoff 係数の MOC は、反射境界の格子を前提にしている。反射体の定数を作る格子（燃料の集合体と反射体を並べたもの）
では、真空やアルベドの境界（[ADR-0093](../adr/0093-choose-lattice-boundary-conditions.md)）の外から入る中性子の扱いを決めていない。

関係する ADR: [ADR-0084](../adr/0084-make-reflector-constants-with-fuel.md)

## 完了の条件

M3 の計画（docs/plans/m3-core-hzp.md）の進め方の 5 で反射体の定数（[ADR-0084](../adr/0084-make-reflector-constants-with-fuel.md)）を作るときに、境界の条件での Dancoff 係数の扱いを決め、テストで確かめる。

## 経緯

- 2026-10-08: 自己遮蔽の前提の洗い出しで挙がった。反射体の段階で扱う（ユーザーと合意）
