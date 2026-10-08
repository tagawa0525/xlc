---
status: accepted
date: 2026-10-08
supersedes: []
superseded-by:
issues: []
---

# ADR-0008: 自前の小さな密行列の解法は faer に置き換えない

## 背景

線形代数は faer とした（[ADR-0007](0007-use-faer-for-linear-algebra.md)）。NEM、UFG、SPH、共鳴の計算、CMFD の群×群のブロックには、自前の密行列の
解法がある。これも faer にするかを、呼ばれる大きさで測って決めた。

- 速さ

## 検討した案

- 自前の密行列の解法を faer に置き換える
- 自前のままにする（採用）

## 決定と理由

自前の密行列の解法は faer にしない。測った値（docs/plans/architecture.md、67fd8fb、2026-10-08）:

- NEM の 3×3、UFG と SPH のリング数の大きさ、共鳴の複素 3×3 まで: 毎回分解して解く大きさ 20 までの行列では、
  faer は呼び出しの手間で自前より 1.4〜17 倍遅い
- CMFD の群×群のブロック（1回分解し、反復ごとに解く）: 69 群の前進・後退代入では faer が 1.2 倍速いが、
  7 群では自前の方が 2 倍速い。差が小さく、群数で逆転するので、faer にする理由にならない

## 確認

上の測定。自前の解法の正しさは `lu_solves_with_pivoting`、`solves_a_system_that_needs_pivoting`、
`solve_rejects_singular_and_overflowing_systems` などで確かめている。

## 結果

疎行列の LU 分解とクリロフ法は faer を使う（[ADR-0007](0007-use-faer-for-linear-algebra.md)）。
