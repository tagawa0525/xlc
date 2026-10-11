---
status: closed
labels: [lattice, core]
opened: 2026-10-11
closed: 2026-10-11
milestone: M3
blocked-by: []
resolution: ADR-0131
---

# 真空の辺のある格子の角で、角の均質な値の大きさが負になる

## 問題

ADR-0125 は、格子の中で core と同じ角の点のバランス（CPB）を解いて、角の点ごとの均質な値の大きさ s_p を求め、
s_p が正でなければエラーにする。外周の真空の辺のある格子で、s_p の熱群が負になった。

- 燃料のピン 2 × 2 の右と上に幅 6 cm の水を置き、右と上を真空にした格子（examples/reflector/vacuum-corner の
  `w6-n*.toml`。69 群のライブラリ、自己遮蔽なし）で、右上の角の点（2 辺が真空）の s_p の熱群は、水のマスを 1 cm、0.5 cm、
  0.25 cm、0.125 cm にすると −2.284、−1.439、−0.519、+0.215 だった。1 辺が真空の右下の角の点でも 0.356〜0.789 で、
  上の辺を鏡面反射にした `w6-n6-top.toml` では、真空の辺と鏡面反射の辺の角の点で 0.039 だった。角型の反射体の格子（examples/reflector/reflector-corner.toml）では、2 辺が真空の
  角の点で（−0.217, −8.18）だった（742ce0ee、2026-10-11）
- はじめに見つけたのは、lattice のテスト `l_shaped_pins_reproduce_the_lattice` の格子の上の辺も真空にしたもの（燃料と水の
  2 群の物質）で、右上の角の点の s_p が（0.245, −0.802）だった（feat/0020-corner-cdf、2026-10-11。テストの中の一時的な出力）

原因は、CPB の真空の辺の流れを r Φ_p（辺の比 r × 角の点の不均質な中性子束）としたことにある。Φ_p は角に接する MOC の
領域の平均で、真空の近くでは点の値より大きく、均質な展開に負の角の値を求めた。展開の限界（拡散の世界では 2 辺が真空の
角でも正の値を与える）でも、#0027 と同じ仕組みでもない。切り分けの数値は ADR-0131 の背景にある。

関係する ADR: [ADR-0125](../adr/0125-derive-corner-factors-from-lattice-corner-balance.md)、[ADR-0131](../adr/0131-derive-corner-factors-with-vacuum-condition-on-homogeneous-corners.md)

## 完了の条件

真空の辺のある格子の角（5-2D の角型の反射体。#0029）で s_p を正にする扱いか、別の角の値の求め方を決めて ADR にする。
5-2D の角型の反射体の格子で、CDF を作れることを示す。

## 経緯

- 2026-10-11: ADR-0125 の実装（feat/0020-corner-cdf）のテストの格子で見つけた。I 型の反射体の格子（真空の辺が 1 つ、
  反射体の幅 21.5 cm）で当たるかは、examples/reflector の結果で確かめる
- 2026-10-11: 原因を切り分け、真空の辺の流れを均質な角の値に r × 面の ADF で課す ADR-0131 にした。小さな格子の系列で
  s_p は 0.92〜1.12、角型の反射体の格子（examples/reflector/reflector-corner.toml）で CDF を作れた
  （examples/reflector/verification.md、9ebcb470、2026-10-11）。I 型の格子はこの扱いでも作れる
