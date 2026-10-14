---
status: open
labels: [xs]
opened: 2026-10-14
closed:
milestone: M4
blocked-by: []
resolution:
---

# Pm-148 の非分離共鳴の自己遮蔽が、10 eV より下の小さい σ0 で NaN になる

## 問題

燃焼鎖の核種のライブラリ（ADR-0169、examples/depletion）を作ると、Pm-148 の段階が `Pm-148 at 565 K: non-finite value NaN` で
止まる（xlc c43f4bbc、2026-10-14、Slurm の排他のジョブ）。Pm-148 は分離共鳴を持たず、非分離共鳴（LRF = 2、LSSF = 0）が 1.1 eV〜100 keV
にあり、1.1 eV での平均の準位間隔 D は 3.36 eV。

565 K で、`unresolved_shielding` の σ0 = 10 b と 1 b の自己遮蔽した全断面積が、10 eV より下の格子点で NaN になる（1.2、2、4、5、8 eV で
NaN、10 eV で 13.5 b（σ0 = 10 b）。手元の確認、c43f4bbc）。無限希釈と σ0 ≥ 30 b は有限。そのため WIMS の 25〜34 群の σ0 = 10 b、1 b の
群定数が NaN になる。`unresolved_shielding` はエラーを返さないので、棚卸しの道具（crates/chain/examples/survey_neutrons.rs）は読めると
数えた。

原因の候補（確かめていない）: 系列どうしの重なりの補正（NJOY unresl と同じ式。crates/xs/src/resonance/unresolved.rs の `shielded`）の
1 − yj が、準位間隔に比べて幅の大きい低いエネルギーで 0 以下になる。NJOY UNRESR が同じ点でどうなるかは確かめていない。

関係する ADR: なし

## 完了の条件

NaN の出どころを確かめ、NJOY UNRESR の同じ点の値と比べて扱いを決める（式の適用の範囲の外をエラーにするか、正しい式にするか）。
`unresolved_shielding` が NaN を返さない（値にするかエラーにする）ことを単体テストで確かめる。そのあとで、Pm-148 をライブラリの
最初の組（examples/depletion/xs.toml）に戻す。

## 経緯

- 2026-10-14: 燃焼鎖の核種のライブラリの最初の組を作ったときに見つかった（feat/depletion-library）。Pm-148 を最初の組から外した
  （ほかの Pm の核種 Pm-147、Pm-148m、Pm-149 は 565 K で NaN がない）
