---
status: open
labels: [xs]
opened: 2026-10-14
closed:
milestone: M4
blocked-by: []
resolution:
---

# 非分離共鳴の自己遮蔽が、低いエネルギーの小さい σ0 で NaN になる（Pm-147、Pm-148、Pm-149、Am-242m）

## 問題

燃焼チェーンの核種のライブラリ（ADR-0169、examples/depletion）を作ると、Pm-148 の段階が `Pm-148 at 565 K: non-finite value NaN` で
止まる（xlc 2441757c、2026-10-14、Slurm の排他のジョブ）。Pm-148 は分離共鳴を持たず、非分離共鳴（LRF = 2、LSSF = 0）が 1.1 eV〜100 keV
にあり、1.1 eV での平均の準位間隔 D は 3.36 eV。

565 K で、`unresolved_shielding` の σ0 = 10 b と 1 b の自己遮蔽した全断面積が、10 eV より下の格子点で NaN になる（1.2、2、4、5、8 eV で
NaN、10 eV で 13.5 b（σ0 = 10 b）。手元の確認、2441757c）。無限希釈と σ0 ≥ 30 b は有限。そのため WIMS の 25〜34 群の σ0 = 10 b、1 b の
群定数が NaN になる。`unresolved_shielding` はエラーを返さないので、棚卸しのツール（crates/chain/examples/survey_neutrons.rs）は読めると
数えた。

原因の候補（確かめていない）: 系列どうしの重なりの補正（NJOY unresl と同じ式。crates/xs/src/resonance/unresolved.rs の `shielded`）の
1 − yj が、準位間隔に比べて幅の大きい低いエネルギーで 0 以下になる。NJOY UNRESR が同じ点でどうなるかは確かめていない。

関係する ADR: なし

## 完了の条件

NaN の出どころを確かめ、NJOY UNRESR の同じ点の値と比べて扱いを決める（式の適用の範囲の外をエラーにするか、正しい式にするか）。
`unresolved_shielding` が NaN を返さない（値にするかエラーにする）ことを単体テストで確かめる。そのあとで、Pm-148 をライブラリの
最初の組（examples/depletion/xs.toml）に戻す。

## 経緯

- 2026-10-14: 燃焼チェーンの核種のライブラリの最初の組を作ったときに見つかった（feat/depletion-library）。Pm-148 を最初の組から外した
  （ほかの Pm の核種 Pm-147、Pm-148m、Pm-149 は 565 K で NaN がない）
- 2026-10-14: ADR-0172 の捕獲を書き出すと、Pm-147、Pm-149（σ0 = 1 b の捕獲）と Am-242m（σ0 = 10 b、1 b の核分裂と捕獲。WIMS の 21〜23 群）
  も NaN になった（565 K。手元の確認で、ほかの計算の排他のジョブと重なって実行してしまった。値の確認だけで時間は記録しない）。どれも
  非分離共鳴が低いエネルギー（Pm-147 は 102 eV、Pm-149 は 2.6 eV、Am-242m は 43 eV から）にある。NJOY UNRESR（unresr.f90 の 1181 行）も
  同じ重なりの式 σ = σm yy/(1 − yj) を、1 − yj を確かめずに使う。全断面積が有限で部分の反応だけが NaN になる Am-242m は、遮蔽した
  値が 0 以下の格子点で、非分離共鳴のパラメータの補間則（ln をとる則）を `UnresolvedShielding::adjust` が当てるのが候補（確かめていない）。
  4 核種を最初の組から外した。Pm-149 は Sm-149 の経路の先行核なので、燃焼の検証（「進め方」5）の前に直す
- 2026-10-14: Eu-151 も 900 K で NaN になった（燃焼のライブラリの段階、Slurm の排他のジョブ 373。565 K は有限）。最初の組から外した
