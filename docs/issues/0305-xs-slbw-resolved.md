---
type: Issue
status: open
labels: [xs, depletion]
opened: 2026-11-20
closed:
milestone:
blocked-by: []
resolution:
---

# xs が分離共鳴の SLBW（LRF=1）を読めず、ENDF/B-VII.0 と VII.1 の燃焼のライブラリを作れない

## 問題

- xs は分離共鳴（LRU=1）の LRF=1（SLBW）を読まず、「resolved LRF=1 is not supported」で止まる（crates/xs/src/endf/mf2.rs）。
  版ごとのピンセルの燃焼（#0303）の 65 核種のライブラリ（examples/depletion/endfb70/xs.toml、endfb71/xs.toml）で、`bin/xs range`
  が次の核種で止まる。VIII.0 と VIII.1 の同じ 65 核種は止まらない
  - VII.0: Pu-238、Pu-242、Am-241、Am-243、Cm-242、Xe-135、Pm-148m（7 核種）
  - VII.1: Xe-135、Pm-148m（2 核種）
  - 再現: `cd examples/depletion/endfb71 && ../../../bin/xs range xs.toml /tmp/range.toml`（最初に当たった核種で止まる。核種ごとには、
    その核種だけを書いた入力で `xs range` を実行する）
- 燃焼チェーンの全核種（ADR-0264 の「確認」の 1 つ目）では、中性子のファイルの全部を 1 つずつ `xs range` に通すと、SLBW で止まる
  ファイルは VII.0 が 393 のうち 24（Nb-93、Rh-105、Xe-135、Pm-148m、Th-230、U-237、U-239、U-240、U-241、Pu-238、Pu-242、Pu-243、
  Pu-244、Am-241、Am-243、Cm-242、Cm-243、Cm-247、Cm-248、Cf-250、Cf-251、Cf-252、Cf-253、Es-253）、VII.1 が 423 のうち 8（Nb-93、
  Rh-105、Xe-135、Pm-148m、U-239、U-240、U-241、Pu-243）。ほかに、VII.0 の天然の V と Zn は複数の同位体（NIS = 2、5）の MF2 で、
  VIII.1 の Sr-88 は分離共鳴のバックグラウンドの R 行列（KBK = 1）で止まる（2026-11-20、`bin/xs` は main 8e88c827 のもの。「no group bounds」で
  止まるファイルは、1 核種だけの入力で分離共鳴の上端が熱の打ち切りより低いためで、読めないのではない。VIII.0 の Pm-148m も同じ）。
  `xs range` は MF2 の分離共鳴の上端までしか読まないので、MF2 のほかの部分と MF3 から MF6 は確かめていない
- NJOY の RECONR は SLBW を再構成する。xs で SLBW を読むか、別の扱いにするかは決めていない

関係する ADR: ADR-0264

## 完了の条件

- [ ] VII.0 と VII.1 の版ごとのライブラリ（examples/depletion/endfb70、endfb71）の 65 核種を xs で処理できる
- [ ] 燃焼チェーンの全核種で、上の形式のどれをどう扱うかを決める（#0303 の完了の条件の 2 つ目と合わせる）

## 経緯

- 2026-11-20: 版ごとのピンセルの燃焼の準備（#0303）で、VII.0 と VII.1 のライブラリの `xs range` が SLBW で止まり、起こした
