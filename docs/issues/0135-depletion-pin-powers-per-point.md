---
status: open
labels: [lattice, depletion]
opened: 2026-10-15
closed:
milestone: M4
blocked-by: []
resolution:
---

# 格子計算の燃焼が燃焼点ごとのピン出力を出さず、集合体の燃焼のピン出力の比較ができない

## 問題

docs/plans/m4-depletion.md の「参照解と比べる量」は、格子計算の燃焼を k∞、ピン出力、主要な核種の数密度で OpenMC の燃焼と比べる
（ピン出力は ADR-0032 の rms 0.2%）。OpenMC 側は、集合体の燃焼（examples/depletion/openmc_depletion_assembly.py。ADR-0178）が燃焼点ごとの
ピン出力（燃料のピンの平均で 1、17×17）を書く。xlc の `[depletion]` は `pin-powers` と `medium` を受け付けず
（crates/lattice/src/input.rs の「depletion writes neither `pin-powers` nor `medium`」）、履歴のファイルは燃焼点ごとの k∞、領域ごとの
数密度、少数群定数だけで、ピン出力を持たない（ADR-0173。燃焼点ごとのピン出力は決めていないとした）。そのため
examples/depletion/compare_depletion_assembly.py の要約は、ピン出力の比較を持たない。

燃焼点ごとのピン出力は、ADR-0096 の形状関数の再構成（分岐のテーブルの燃焼度の軸）にも要る見込みで、燃焼の履歴にピン出力を書くかどうか
（`pin-powers` を燃焼点ごとのファイルにするか、履歴に入れるか）は、計画の「進め方」6 の ADR と一緒に決める候補。

関係する ADR: ADR-0178

## 完了の条件

xlc の燃焼の履歴か別のファイルに、燃焼点ごとのピン出力（OpenMC と同じ並び: 上の行から、燃料のピンの平均で 1）を書き、
compare_depletion_assembly.py がピン出力の rms と最大の差を OpenMC の統計誤差とともに要約する。

## 経緯

- 2026-10-15: 集合体の燃焼の入力（ADR-0178）を作る中で、比較の道具にピン出力がないことが分かった
