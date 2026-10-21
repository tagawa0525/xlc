---
status: closed
labels: [lattice, depletion]
opened: 2026-10-15
closed: 2026-10-22
milestone: M4
blocked-by: []
resolution: ADR-0181
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

関係する ADR: ADR-0178、ADR-0181

## 完了の条件

xlc の燃焼の履歴か別のファイルに、燃焼点ごとのピン出力（OpenMC と同じ並び: 上の行から、燃料のピンの平均で 1）を書き、
compare_depletion_assembly.py がピン出力の rms と最大の差を OpenMC の統計誤差とともに要約する。

## 経緯

- 2026-10-15: 集合体の燃焼の入力（ADR-0178）を作る中で、比較のツールにピン出力がないことが分かった
- 2026-10-22: ピン出力を燃焼の履歴の各点に書くことにした（ADR-0181）。値は `pin-powers` と同じ定義で、その点の組成で解いた輸送（予測子の
  始まりと最後の点）のもの。compare_depletion_assembly.py は燃焼点ごとに、xlc と OpenMC のピン出力を対称操作の組で平均した相対差の rms、
  最大（位置と、その組の OpenMC の統計誤差）を `pins` に書き、`--report` に描く。記録の表は verification_tables.py の
  `depletion-<問題>-pins`。2A の 1 ステップ（examples/depletion/lattice-2A.toml と同じ入力、コミット 1746edc5、分岐点 fcb2289c、Slurm で 2 分 48 秒）
  で、履歴の 2 点に 17×17 のピン出力が入った（核分裂のある 264 ピンの平均が 1、案内管と計装管が 0、最小 0.940、最大 1.050（燃焼度 0）、
  対角線と上下の鏡映で 3e-15 以内の対称）。k∞（1.182836、1.145558）は main の履歴と同じ。履歴の大きさは 2 点で 5.3 KB 増えた（2.74 MB）。
  OpenMC の集合体の燃焼の結果（references/openmc/aws/m4-assembly/）はまだないので、比較のツールは、この履歴と、OpenMC の結果の形をした
  合成のファイル（数密度と k は xlc の値、ピン出力は examples/vera/p2/openmc-2A.msgpack、統計誤差はピンごとに 0.5% と置いたもの）で
  `--json` と `--report`、表を通した（数値は記録に使わない）
