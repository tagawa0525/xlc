---
status: open
labels: [lattice, examples]
opened: 2026-11-04
closed:
milestone: M3
blocked-by: []
resolution:
---

# ノードごとの制御棒の軸の付け分けの残り（p5 の棒のない状態の照合、軸の有無の判定の重複）

## 問題

この issue は、#0213 の後半（ADR-0227。ブランチ feat/0213-node-rod-axis）のレビュー（/code-review high、2026-11-04）で出た、直さずに残した指摘

- examples/vera/p5 の制御棒の分岐を解く 8 本の格子は、棒のない設計（と I 型の反射体の格子の反射体のノードの定数を使う径方向の
  反射体の列）も読む。入力は前の棒のない格子から、単一集合体の 6 本では使わない物質 ss304、aic、b4c と減速材の軸の `materials` のそれらの物質（[zircaloy, water]
  から [zircaloy, water, ss304, aic, b4c]）、I 型の反射体の格子の 2 本では aic と b4c（ss304 はバッフルで前からあり、[zircaloy, water, ss304]
  から [zircaloy, water, ss304, aic, b4c]）、置き換えるセル（`cr-aic`、`cr-b4c`）と `[branches.rod]` が増えた。棒のない値の状態が前の棒のない格子と同じになることは、2 群の合成のライブラリの小さい格子
  （`nodes_without_replaced_cells_hold_the_rodless_states`。比べる別の入力は `[branches.rod]` だけを持たず、使わない物質と置き換えるセルは
  両方にある）と、examples/core
  の集合体（ADR-0224 の確認）でしか確かめていない。p5 を再実行するとき、置き換えた格子の値 out の状態と、同じ bin とライブラリで
  解いた前の形の棒のない格子（分岐点 7cd58518 の hzp/<名前>.toml）の状態を比べる（反射体の格子を含む）。2026-11-04 の試行
  （Slurm 685）は 21-fuel を解いただけで、比べていない
- `without_rod_axis`（crates/lattice/src/branches.rs）は `micro` と `histories` をそのまま渡し、値の番号 0 が `out` であることを
  確かめない。今は、ミクロ燃焼のデータはこの後に付け、ノードの出力は履歴からの分岐と組めず（入力の検査）、値の番号 0 は
  `Branches::rod_values` が `out` にするので通らない
- ノードの範囲が格子を覆うかは、状態を解いた後の `homogenize_nodes` で確かめ、入力の読み込みの制御棒の検査が先に走る。覆わない
  並びでは「燃料のノードが置き換えるセルを含まない」のエラーになり、本当の誤りを示さない。単体テストの `no cell` の場合も、
  覆わない並びで確かめている
- 制御棒の軸の有無を決める箇所が増えた（`LatticeInput::rod_axes`、入力の検査、main.rs の軸を除く `rods > 1`。#0214 の一覧）。
  すべてのノードに軸を入れてから除く形でなく、出力ごとの有無から `tables` でテーブルの軸を作れば、除く処理と長さの照合が要らない。
  `nodes_with_replaced_cells` を入力の検査と `rod_axes` で 2 回求め、`Node` の列を複製する形も 3 か所にある
- examples/vera/p5/p5_inputs.py の `fuel_file` の `rodded` は使われず、`hzp_lattices` は同じパスを返す 2 つの呼び出しを比べる

関係する ADR: なし

## 完了の条件

p5 の再実行で値 out の状態を前の形の棒のない格子と比べ、記録に書く。ほかの項目は、直すかを決める。

## 経緯

- 2026-11-04: #0213 の後半のレビューで出た指摘のうち、落ちるテストか誤った数値を示せないものをまとめた
- 2026-11-05: 完了の条件の照合を行った（ブランチ chore/0213-rerun-p5-hzp）。前の形の棒のない格子は、置き換える前の分岐点 7cd58518 の
  `hzp/<名前>.toml` と出力のパスだけが違う入力（p5_inputs.py が `rodded` を偽にして書く `hzp/rodless-<名前>.toml`。8 本で、単一集合体の
  燃料と支持格子の層の 6 本と I 型の反射体の格子の 2 本）を、同じ bin とライブラリで解く常設のステージにし（`hzp-rodless-*`）、
  `compare_rodless.py`（ステージ `hzp-compare-rodless`）で点ごとに比べた。13 の出力の 104 状態がすべてビットまで同じ（Slurm 695。
  examples/vera/verification.md の「Problem 5 の 3D の炉心の計算」の「制御棒の分岐の格子の棒のない状態」）。ほかの項目は、直すかを
  決めていない。照合のレビューの残りは #0229
