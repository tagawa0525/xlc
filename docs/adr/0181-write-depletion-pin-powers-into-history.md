---
status: accepted
date: 2026-10-22
requires: [ADR-0096, ADR-0173, ADR-0176, ADR-0178]
supersedes: []
superseded-by:
issues: ["#0135"]
---

# ADR-0181: 格子計算の燃焼の燃焼点ごとのピン出力を、燃焼の履歴の各点に書く

## 背景

格子計算の燃焼（[ADR-0173](0173-deplete-lattice-pins-by-ce-cm.md)）の履歴は、燃焼点ごとの時間、燃焼度、k∞、領域ごとの数密度、少数群定数だけを持ち、
ピン出力を持たない。入力の `pin-powers`（セルごとの核分裂率のファイル `xlc-pin-powers`）は `[depletion]` と使えない。集合体の燃焼の検証
（[ADR-0178](0178-verify-assembly-depletion.md)）は、OpenMC の燃焼点ごとのピン出力（17×17、上の行から、燃料のピンの平均で 1）を xlc と比べる予定で、
比較のツール（examples/depletion/compare_depletion_assembly.py）は OpenMC の統計誤差しか書けていない（#0135）。

決めることは、燃焼点ごとのピン出力をどこに書くか、どの輸送の値にするか、比較で何を要約するかである。判断の軸:

- 比較のツールと DVC のステージが読む手間（examples/depletion/p2/dvc.yaml の `lattice` のステージは重く、出力の並びを変えるとステージの定義が変わる）
- 燃焼度の軸を持つ分岐のテーブル（docs/plans/m4-depletion.md の「進め方」6）での使い回し。[ADR-0096](0096-reconstruct-pin-power-analytically.md) の形状関数が燃焼度ごとに要る
- 履歴のファイルの大きさ
- OpenMC の参照解（examples/depletion/openmc_depletion_assembly.py の `bos_points`）と同じ定義

## 検討した案

置き場所:

- 案 H1: 履歴の各点に `pin_powers` を足す（採用）。比較のツールは今の履歴 1 つを読むだけで済み、p2/dvc.yaml のステージの出力も deps も変わらない。
  大きさは 289 ピンで 1 点あたり約 2.6 KB（MessagePack の float64 は 9 バイト）。数密度は 3,820 核種 × 領域で、2A の 1 ステップ（2 点、39 領域）の
  履歴が 2.7 MB（examples/depletion/lattice-2A.msgpack、2026-10-21）なので、ピン出力は 1 点の大きさの約 0.2%
- 案 F1: 燃焼点ごとの `xlc-pin-powers` のファイル（`pin-powers` の名前に点の番号を付ける）: 既存の形式と読み手（examples/vera/compare_pins.py など）を
  そのまま使えるが、40〜70 点のファイルをステージの出力と比較の deps に並べることになり、p2/dvc.yaml の重い `lattice` のステージの定義が変わる。
  燃焼度との対応をファイルの名前で保つことになる
- 案 F2: 全点のピン出力を 1 つの別のファイルにする: F1 の本数の問題はないが、ステージの出力が 1 つ増え、履歴と同じ燃焼度の列を 2 か所に持つ

書くかどうか:

- 常に書く（採用）: 履歴の大きさに比べて小さく、燃焼する格子にはかならず核分裂がある（出力の規格化が核分裂の出力を求める）
- 入力の切り替え（`[depletion]` の `pin-powers = true`）: 書かない利点がなく、入力の選択肢を 1 つ増やす

輸送の値: 燃焼点の値は、その点の組成で解いた輸送（CE/CM の予測子の始まりの輸送と、最後の点の輸送）とする（採用）。k∞ と少数群定数と同じ
輸送で、OpenMC の `bos_points`（燃焼点の始まりの輸送の statepoint `openmc_simulation_n<点>.h5`）と同じ点の定義になる。中点（修正子の
始まり）の輸送は点の組成のものでないので使わない。

比較の要約（compare_depletion_assembly.py）:

- 案 C1: 対称操作の組（xlc が 1/8 で解くときに 1 つにまとめるピンの組。ADR-0178 の OpenMC の物質の組と同じ）ごとに、xlc と OpenMC の
  ピン出力を組で平均し、組の相対差 (xlc / OpenMC − 1) の、核分裂のある組の二乗平均平方根と、絶対値が最大の組を書く（採用）。
  OpenMC の統計誤差の rms（`pin_std_rms`）は組の重みなしの二乗平均平方根なので、差の rms と同じ重みで比べられる。最大の組には、
  その組の OpenMC の統計誤差を添える
- 案 C2: 燃料のピンすべての相対差の rms（examples/vera/compare_openmc_keno.py の xlc − OpenMC と同じ）: 組の中で xlc の値は同じなので、
  組のピンの数（4 か 8）で重みを付けることになり、統計誤差の rms と重みが揃わない

## 決定と理由

- 置き場所は案 H1。履歴の各点に `pin_powers` を書く。値は既存の `pin-powers` と同じ定義（`power::cell_fission_rates`: セルごとの
  Σ_領域 Σ_g Σf φ V を、核分裂のあるセルの平均が 1 になるよう規格化）で、並びも同じ（`xlc-pin-powers` の `rows`: ピンセルのある行と列
  だけ、行は上から、各行は左から）。燃焼しない棒（案内管、計装管、AIC）と、核分裂のないセルは 0。MOC は対称操作で 1/8 を解くが、
  中性子束は格子のすべての領域に展開してあるので、展開の仕方も既存の `pin-powers` と同じ
- 常に書く。入力の `pin-powers` と `medium` は、今も `[depletion]` と使えない（ピン出力は履歴に入る）
- 燃焼点の値は、その点の組成で解いた輸送（予測子の始まりと最後の点）
- 比較の要約は案 C1。燃焼点ごとに、組の相対差の rms、最大（位置と OpenMC の統計誤差）、組の数を `pins` に書く。記録の表は
  examples/verification_tables.py の `depletion-<問題>-pins`（燃焼度ごとの rms、最大 ± 統計誤差と位置、OpenMC の統計誤差の rms）で描く
  （[ADR-0176](0176-generate-verification-tables.md)）。燃焼度の点の数が問題で違う（40 点と 70 点）ので、問題ごとに 1 つの表にする
- 燃焼度の軸を持つテーブル（「進め方」6）で再構成に使うのは、ピン出力でなく群ごとの形状関数（ADR-0096）で、`[depletion]` の
  `few-group` の `pins = true` が燃焼点ごとに少数群定数の `pins` として履歴に書く。燃焼点ごとの値を履歴の同じ点に置く今の形
  （k∞、少数群定数、形状関数）に、ピン出力も揃える

## 確認

- 単体のテスト（crates/lattice/src/depletion/tests.rs）: `pin_powers_are_reported_at_each_point`（2×2 の格子の 1 つのセルを水にした
  格子で、各点のピン出力が 2×2、核分裂のあるセルの平均が 1、水のセルが 0、対称操作で移り合うピンが同じ値、水に面したピンが高い）、
  `the_pin_powers_of_a_point_come_from_its_own_composition`（初めの点のピン出力が、最初のステップの長さによらない。中点の輸送の値なら
  ステップの長さで変わる）
- 比較のツール（examples/depletion/test_compare_depletion_assembly.py）: `test_pin_powers_are_compared_as_means_over_each_symmetry_orbit`
  （組の平均の相対差の rms と最大、位置、最大の組の統計誤差）、`test_pins_with_fission_in_other_cells_stop_the_comparison`、
  `test_a_history_without_pin_powers_stops_the_comparison`、`test_pin_maps_of_other_shapes_stop_the_comparison`、
  `test_the_report_has_the_pin_power_differences`。表（examples/test_verification_tables.py の `DrawAssemblyPins`）: 燃焼度ごとの行と、問題ごとの表の名前
- 動作: 2A の 1 ステップ（examples/depletion/lattice-2A.toml と同じ入力。コミット 1746edc5、2026-10-22）で、履歴の 2 点に 17×17 のピン出力が
  入り（核分裂のある 264 ピンの平均が 1、案内管と計装管が 0、対角線と上下の鏡映で 3e-15 以内の対称）、k∞ は main の履歴と同じ。比較のツールは、
  この履歴と OpenMC の結果の形をした合成のファイルを読んだ（OpenMC の集合体の燃焼の結果はまだない。詳細は #0135 の経緯）

## 帰結

- 燃焼の履歴の形式が変わる。今までの履歴（`pin_powers` のないもの）は compare_depletion_assembly.py で読めない（止まる）。ピンセルの比較
  （compare_depletion.py）はピン出力を読まないので、今までの 1C の履歴もそのまま読める
- 組の rms は、examples/vera の集合体の xlc − OpenMC（燃料のピンすべての rms）と重みが違う。ADR-0032 の目標（rms 0.2%）と比べるときは、
  組で数えたこと（17×17 の 1/8 で核分裂のある組は、2A で 39）を記録に書く
