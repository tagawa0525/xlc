---
status: open
labels: [lattice, core]
opened: 2026-11-07
closed:
milestone:
blocked-by: []
resolution:
---

# ホウ素の平均の原子質量の式が 2 か所にある

## 問題

- ホウ素の平均の原子質量 b M_B10 + (1 − b) M_B11（ADR-0106）を、crates/format/src/table.rs の `b10_density_ratio`（core の B-10 の割合の
  換算。ADR-0237）と crates/lattice/src/water.rs の `borated_water`（ホウ酸水の組成）が、それぞれ書く。原子質量の定数は 1 か所
  （`xlc_format::table::MASS_B10`、`MASS_B11`）で、式だけが 2 か所にある（crates/format/src/table/tests.rs の
  `b10_density_ratio_follows_the_atoms_per_boron_mass` が期待値として書く式は数えない）

関係する ADR: なし

## 完了の条件

- 式を 1 か所にまとめるか、まとめない理由を書いて閉じる

## 経緯

- 2026-11-07: ADR-0236 と ADR-0237 の実装のレビューで出た点のうち、外した ADR-0236 の実装に関わらないもの
