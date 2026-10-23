---
status: open
labels: [verification, openmc]
opened: 2026-10-12
closed:
milestone: M3
blocked-by: []
resolution:
---

# 5A-2D の OpenMC の参照解が対角の鏡映で対称でない（集合体で最大 0.6%）

## 問題

5A-2D の装荷図と炉心の外の構造（中性子パッドは 45° を中心に置く）は対角の鏡映で対称だが、OpenMC の参照解（da3fd6ec、
2026-10-12。examples/vera/verification.md の「Problem 5-2D の参照解」）の集合体出力は、対角で移り合う位置の間で最大 0.60%
違う。統計誤差（集合体の出力の標準偏差）の 12 倍で、統計では説明できない。8 行（対称面 y = 0 の上の集合体）が H 列（対称面
x = 0 の上）より系統的に高く、外へ行くほど差が大きい（8 行 / H 列 − 1 が、G から A の向きに 0.11、0.08、0.09、0.21、0.54、0.60、
0.49%）。ピン出力では対角の鏡映の差の RMS 0.45%、最大 2.1%。extract（`openmc_core_2d.py`）の結果のディレクトリ
（examples/vera/p5/openmc-5A-2D-results。ステージ openmc-5A-2D-extract の出力）の `assemblies.msgpack` と `pins.msgpack` を転置したものと比べた（2026-10-12）。

xlc の 5A の解は 90° の回転で対称なので、この差は xlc との比較にそのまま入る（A08 と H15 の比較の差が 0.5% 違う。
examples/vera/verification.md の「Problem 5-2D の炉心の計算」、#0097）。

原因の候補（どれも確かめていない）:

- 源の収束: 非アクティブバッチ 300 で、大きな 2D の炉心の高い優勢比のため、対角の向きの傾きが残っている
- モデルの非対称: x = 0 と y = 0 の対称面で切れた集合体の扱い（extract で全体に直す係数）、メッシュの tally の境界、
  バッフルやパッドの形が対角で対称でない

関係する ADR: なし

## 完了の条件

原因を確かめ、モデルの誤りなら直して再実行する。源の収束なら、非アクティブバッチか粒子数を見直すか、差の大きさを記録の統計誤差に
加えて扱う。

## 経緯

- 2026-10-12: #0086 で xlc の 5A を OpenMC と比べたときに見つけた（feat/0086-periphery-constants）
