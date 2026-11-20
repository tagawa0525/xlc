---
type: ADR
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: ["#0308"]
requires: []
---

# ADR-0055: R-matrix limited 形式（LRF=7）は Fe-54, 57 が使う範囲だけを NJOY RECONR と同じ式で扱う

## 背景

Zircaloy-4 の天然の同位体を1つずつ処理して、足りない機能を洗い出した（[ADR-0054](0054-model-vera-materials.md)）。Fe-54 と Fe-57 の分離共鳴は
R-matrix limited 形式（LRF=7）である。除くと VERA と同じ組成で比べられなくなる。

## 検討した案

- Fe-54, 57 を除く: VERA と同じ組成で比べられない
- LRF=7 を一般に実装する: 対象核種で使わない部分まで作ることになる
- Fe-54, 57 が使う範囲だけを実装し、それ以外はエラーにする（採用）

## 決定と理由

NJOY RECONR と同じ式で実装する。Fe-54, 57 が使う範囲だけを扱う: Reich-Moore（KRM=3）、幅で与えるパラメータ（IFG=0）、
バックグラウンドの R 行列と位相シフトの指定なし（KBK=0, KPS=0）、電荷のない中性子の粒子対（剛体球の透過因子、シフト因子なし）、
核分裂の粒子対（MT=18、透過因子は 1）と、消去する捕獲のチャネル。それ以外はエラーにする。

- 捕獲のチャネルは粒子対の中身（光子、MT=102）で見分ける。NJOY は1番目の粒子対を捕獲とみなすが、ENDF は順序を
  定めていない。捕獲の幅は共鳴エネルギーの虚部（E_λ − iΓγ/2）に入れる
- 断面積は、開いた中性子と核分裂のチャネルの X = (I − iK)⁻¹K（K = P^½ R P^½）から求める。閾値より下の非弾性のチャネルは
  閉じているので除く。弾性散乱・捕獲・準位ごとの非弾性散乱（MT51〜）・核分裂を、MF3 のバックグラウンドに足す（NJOY RECONR と同じ）
- 核分裂のチャネルの透過因子は 1 とする（ENDF-102 の D.1.7）。幅は Γ = 2γ² で、LRF=3 の Reich-Moore の Γfa、Γfb と同じチャネルになる。
  粒子対の PNT は -1 と既定の 0 を受け、1 はエラーにする。2.2.1.6 は PNT=0 の既定を「fission (MT=19)」で -1、ほかの MT で +1 と書くが、
  剛体球の透過因子は核分裂に定義がなく、D.1.7 と NJOY2016 samm（`ppdefs` が MT=18 を `lpent=0` にする）に合わせる。MT=18 の粒子対が
  2 つ以上あればエラーにする（準位の MT の重なりと同じ）
- 幅から換算幅への換算は ENDF のとおり、共鳴エネルギーでの透過因子（負の共鳴は |E_λ|）を使う。共鳴エネルギーで
  閉じているチャネルに幅があれば、定義できないのでエラーにする

## 確認

NJOY RECONR の 0 K のポイントワイズの断面積と比べた。Fe-54, 57 の MT1, 2, 102, 51 は許容誤差（0.1%）の中で一致した（Fe-54 の MT51 の
閾値のすぐ上の1点だけ、NJOY が閾値を7桁目で1つ上げるために 0.5% 違う。docs/plans/m2-assembly.md、8196fbd、
2026-10-07）。ENDF/B-VIII.1 の Pu-239 の核分裂の粒子対は examples/depletion/verification.md の「版ごとのライブラリの分離共鳴の
R 行列の制限形の核分裂」。テストは `r_matrix_limited_adds_level_inelastic_scattering`、`rejects_r_matrix_limited_other_than_reich_moore`、
`rejects_charged_pairs`、`rejects_shift_factor`、`rejects_width_of_a_channel_closed_at_the_resonance`、
`fission_channel_matches_single_level_breit_wigner`、`two_fission_channels_match_reich_moore`、`rejects_penetrability_for_fission`、
`rejects_two_fission_pairs`、`r_matrix_limited_fission_goes_to_the_fission_mt`。

## 帰結

以上（ADR-0055〜[ADR-0058](0058-read-lrf1-as-lrf2-with-l3.md)）で、Zircaloy-4 の天然の同位体（Zr, Sn, Fe, Cr の 23 核種）をすべて処理できる。Zr と Sn の
ほかの同位体と Fe-56 はそのまま処理できた。

## 変更

- 2026-11-20: 「決定と理由」の扱う範囲に核分裂の粒子対を足した。前は「電荷のない中性子の粒子対（剛体球の透過因子、シフト因子なし）と、
  消去する捕獲のチャネル」と「開いた中性子のチャネルの X = (I − iK)⁻¹K」と「弾性散乱・捕獲・準位ごとの非弾性散乱（MT51〜）を、MF3 の
  バックグラウンドに足す」。ENDF/B-VIII.1 の Pu-239 が LRF=7 で核分裂の粒子対を持ち、VIII.1 の燃焼のライブラリが作れないため（#0308）
