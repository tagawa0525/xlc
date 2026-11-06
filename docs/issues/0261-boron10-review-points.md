---
status: open
labels: [lattice, core]
opened: 2026-11-07
closed:
milestone:
blocked-by: []
resolution:
---

# ADR-0237 の実装のレビューで直さなかった点

## 問題

ADR-0237（B-10 の割合の ppm の換算）の実装のレビューで出た、直していない点。

- ホウ素の平均の原子質量 b M_B10 + (1 − b) M_B11（ADR-0106）を、crates/format/src/table.rs の `b10_density_ratio`（core の B-10 の割合の
  換算。ADR-0237）と crates/lattice/src/water.rs の `borated_water`（ホウ酸水の組成）が、それぞれ書く。原子質量の定数は 1 か所
  （`xlc_format::table::MASS_B10`、`MASS_B11`）で、式だけが 2 か所にある（crates/format/src/table/tests.rs の
  `b10_density_ratio_follows_the_atoms_per_boron_mass` が期待値として書く式は数えない）
- examples/depletion/compare_core_b10.py の `rho_pcm` は、compare_core_cycle.py、compare_core_micro.py、compare_core_partial_rods.py、
  compare_core_restart.py にもある同じ関数の 5 つ目で、core の標準出力の `k-eff` の行の読み方（`K_LINE`）も compare_core_boron.py と重なる
- compare_core_b10.py の `summarize` は、直接の格子計算の点に当たる状態がケースに無いと、`SystemExit` でなく `KeyError` で止まる（nominal の
  状態が無いときも同じ）。`table_k` は、テーブルの軸に無い ppm か燃焼度で `ValueError` で止まる
- 原子質量の定数の置き場が分かれる: B-10 と B-11 は `xlc_format::table`、H-1 と O-16 は crates/lattice/src/water.rs
- crates/core/src/state/tests.rs の `table_b10_is_not_converted` は、割合が同じときの比が厳密に 1 なので、割合が同じでも換算する実装でも落ちない

関係する ADR: なし

## 完了の条件

- 点ごとに、直すか、直さない理由を書いて閉じる

## 経緯

- 2026-11-07: ADR-0236 と ADR-0237 の実装のレビューで出た点のうち、外した ADR-0236 の実装に関わらないもの。rebase の後の親の周の点を足した
