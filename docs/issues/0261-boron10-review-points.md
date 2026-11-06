---
status: open
labels: [lattice, core]
opened: 2026-11-07
closed:
milestone:
blocked-by: []
resolution:
---

# ADR-0236 の実装のレビューで直さなかった点

## 問題

ADR-0236（冷却材の B-10 の割合の補正）の実装のレビューで出た、直していない点。

- core の補正は拡散係数を変えない（ADR-0236 の「core の補正」）。B-10 の個数密度の差は熱群の全断面積を変えるので、漏れのある炉心では
  拡散係数と吸収が食い違う。集合体 1 体（反射境界）の確認（examples/depletion/verification.md の「冷却材の B-10 の割合がテーブルと違う状態
  （2A の集合体 1 体）」）では拡散係数が効かず、その大きさは測っていない
- `Coolant::micro`（crates/lattice/src/coolant.rs）は、ノードの出力ごとに SPH 因子、群の対応、核種の群定数を作り直す。ノードの数だけ
  繰り返す
- チェックポイントの出力の形が合わないときのエラー（crates/lattice/src/main.rs の `fits`）は、`boron10` の有無を書かない。`boron10` だけが
  違うと、両側に同じ値が出る
- examples/depletion/compare_core_b10.py は、B-10 の割合を変えた効きが 0 の点（0 ppm、テーブルと同じ割合）で 0 で割る。状態の名前の B-10 の
  割合は 3 桁（`b198`）までで、19.78 at% を書けない
- `boron10_density_ratio`（crates/format/src/table.rs）の、割合が同じときの早い返り値は、比が 1 になるので要らない。ホウ素の平均の原子質量の
  式は、`boron10_density_ratio`、crates/lattice/src/water.rs の `borated_water`、examples/depletion/compare_core_b10.py の `per_mass` の 3 か所に
  あり、compare_core_b10.py は B-10 と B-11 の原子質量も手でコピーしている
- `Coolant::micro` は、ホウ酸水の領域がなければ 0 を返す。ノードでなく格子全体の出力でも、ホウ素の軸の物質をどの領域にも置かない入力は、
  エラーでなく 0 の補正になる。0 の早い返り値の前に、使わない SPH 因子を作る
- 媒質を置いた入力は、どの領域にも使わない媒質でも `boron10` を書かない（crates/lattice/src/main.rs の `coolant_waters`）
- ホウ酸水のないノードを持つノードごとの出力を、格子計算の実行で確かめるテストがない（`Coolant::micro` を直接呼ぶテストだけ）
- #0260 のファイル名は、題名に足した媒質を含まない

関係する ADR: なし

## 完了の条件

- 上の点ごとに、直すか、直さない理由を書いて閉じる

## 経緯

- 2026-11-07: ADR-0236 の実装のレビューの 1 周目と 2 周目で出た
