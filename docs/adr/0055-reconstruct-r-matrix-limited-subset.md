---
type: ADR
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: ["#0308", "#0310"]
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
位相シフトの指定なし（KPS=0）、バックグラウンドの R 行列はなしか SAMMY の対数の形（KBK の記録の LBK=0 と 2）、電荷のない中性子の
粒子対（剛体球の透過因子、シフト因子は SHF=1 の粒子対だけ）、核分裂の粒子対（MT=18、透過因子は 1）と、消去する捕獲のチャネル。
それ以外はエラーにする。

- 捕獲のチャネルは粒子対の中身（光子、MT=102）で見分ける。NJOY は1番目の粒子対を捕獲とみなすが、ENDF は順序を
  定めていない。捕獲の幅は共鳴エネルギーの虚部（E_λ − iΓγ/2）に入れる
- 断面積は、開いた中性子と核分裂のチャネルの X = (I − iK)⁻¹K（K = P^½ R P^½。シフト因子は下の項）から求める。閾値より下の非弾性のチャネルは
  閉じているので除く。弾性散乱・捕獲・準位ごとの非弾性散乱（MT51〜）・核分裂を、MF3 のバックグラウンドに足す（NJOY RECONR と同じ）
- 核分裂のチャネルの透過因子は 1 とする（ENDF-102 の D.1.7）。幅は Γ = 2γ² で、LRF=3 の Reich-Moore の Γfa、Γfb と同じチャネルになる。
  粒子対の PNT は -1 と既定の 0 を受け、1 はエラーにする。2.2.1.6 は PNT=0 の既定を「fission (MT=19)」で -1、ほかの MT で +1 と書くが、
  剛体球の透過因子は核分裂に定義がなく、D.1.7 と NJOY2016 samm（`ppdefs` が MT=18 を `lpent=0` にする）に合わせる。MT=18 の粒子対が
  2 つ以上あればエラーにする（準位の MT の重なりと同じ）
- 幅から換算幅への換算は ENDF のとおり、共鳴エネルギーでの透過因子（負の共鳴は |E_λ|）を使う。共鳴エネルギーで
  閉じているチャネルに幅があれば、定義できないのでエラーにする
- シフト因子は ENDF-102 の式 (D.42)〜(D.47) の L = (S − B) + iP で入れる（B は粒子対のチャネルの BND）。K = P^½ R P^½ と
  Λ = L/P（対角）で X = (I − KΛ)⁻¹K を解き、シフト因子を計算するチャネルのない spin group は Λ = i のまま (I − iK)X = K を解く。
  SHF=0 と −1 は計算しない（−1 は ENDF-102 の 2005 年版の「計算しない」、2023 年版は 0。ENDF/B-VII.1 の Cl-35 が −1）。
  SHF=1 の捕獲と核分裂の粒子対は、透過因子と同じくシフト因子にも定義がないのでエラーにする。幅から換算幅への換算は
  シフト因子によらない（式 D.45）。共鳴のまわりの格子点は、準位のずれ −Σγ²(S − B)（S は共鳴エネルギーで求める）だけ動かした
  エネルギーに置く
- バックグラウンドの R 行列の LBK=2 は、式 (D.77) の実数を、開いたチャネルの R の対角に足す（エネルギーは実験室系。NJOY2016 samm の
  `setr` と同じ）。範囲 (ED, EU) の外は対数が定義できないのでエラーにする。消去する捕獲のチャネルのバックグラウンドはエラーにする
  （R 行列に捕獲のチャネルが無い）。LBK=1（表）と 3（Fröhner の形）、KPS≠0 は、4 つの版の評価が使わないのでエラーのままにする

## 確認

NJOY RECONR の 0 K のポイントワイズの断面積と比べた。Fe-54, 57 の MT1, 2, 102, 51 は許容誤差（0.1%）の中で一致した（Fe-54 の MT51 の
閾値のすぐ上の1点だけ、NJOY が閾値を7桁目で1つ上げるために 0.5% 違う。docs/plans/m2-assembly.md、8196fbd、
2026-10-07）。ENDF/B-VIII.1 の Pu-239 の核分裂の粒子対は examples/depletion/verification.md の「版ごとのライブラリの分離共鳴の
R 行列の制限形の核分裂」。テストは `r_matrix_limited_adds_level_inelastic_scattering`、`rejects_r_matrix_limited_other_than_reich_moore`、
`rejects_charged_pairs`、`rejects_shift_factor_for_capture`、`rejects_unknown_shift_factor_flag`、`rejects_width_of_a_channel_closed_at_the_resonance`、
`fission_channel_matches_single_level_breit_wigner`、`two_fission_channels_match_reich_moore`、`rejects_penetrability_for_fission`、
`rejects_two_fission_pairs`、`r_matrix_limited_fission_goes_to_the_fission_mt`。
ENDF/B-VIII.1 の Cu-63（シフト因子）と Sr-88（LBK=2 のバックグラウンド）は examples/depletion/verification.md の「版ごとのライブラリの分離共鳴の
R 行列の制限形のシフト因子とバックグラウンド」（#0310）。テストは `shift_factor_enters_the_level_matrix`、
`shift_factor_flag_minus_one_is_not_calculated`、`pole_with_shift_factor_is_the_limit`、`background_r_matrix_adds_to_the_channel`、
`background_r_matrix_with_shift_factor`、`peaks_of_shifted_levels_sit_at_the_observed_energy`、`reads_sammy_logarithmic_background_of_a_channel`
ほか。

## 帰結

以上（ADR-0055〜[ADR-0058](0058-read-lrf1-as-lrf2-with-l3.md)）で、Zircaloy-4 の天然の同位体（Zr, Sn, Fe, Cr の 23 核種）をすべて処理できる。Zr と Sn の
ほかの同位体と Fe-56 はそのまま処理できた。

## 変更

- 2026-11-20: 「決定と理由」の扱う範囲に核分裂の粒子対を足した。前は「電荷のない中性子の粒子対（剛体球の透過因子、シフト因子なし）と、
  消去する捕獲のチャネル」と「開いた中性子のチャネルの X = (I − iK)⁻¹K」と「弾性散乱・捕獲・準位ごとの非弾性散乱（MT51〜）を、MF3 の
  バックグラウンドに足す」。ENDF/B-VIII.1 の Pu-239 が LRF=7 で核分裂の粒子対を持ち、VIII.1 の燃焼のライブラリが作れないため（#0308）
- 2026-11-20: 「決定と理由」の扱う範囲にシフト因子（SHF=1、−1）と LBK=0、2 のバックグラウンドの R 行列を足し、その式の項を
  足した。前は「バックグラウンドの R 行列と位相シフトの指定なし（KBK=0, KPS=0）、電荷のない中性子の粒子対（剛体球の透過因子、
  シフト因子なし）」。4 つの版の LRF=7 の 28 評価のうち、ENDF/B-VIII.1 の Cu-63（弾性散乱の SHF=1）と Sr-88（7 つの spin group の
  KBK=1、LBK=2）が止まり、VII.1 の Cl-35 が SHF=−1 で止まっていたため（#0310）
  比べた案: すべての spin group を X = P^½ (I − RL)⁻¹ R P^½ の 1 つの形で解く案は、シフト因子のない評価の出力の丸めが変わりうる（確かめていない）
  ので採らず、Λ = i の経路を残した（VIII.0、VIII.1 の LRF=7 の 9 評価で出力のバイトが同じ）。格子点を NJOY のように E_λ ± 半値幅に置くと、
  Cu-63 の 69941.6 eV の p 波の準位は 10 eV（半値幅の 2.6 倍）低い所にピークがあり、裾に 2.4 eV の格子の隙間が残って全断面積が
  NJOY と 4.0e-3 違ったので、準位のずれの分だけ動かした
- 2026-11-20: 扱う範囲に、荷電粒子の放出の粒子対（MT600〜849。Coulomb の透過因子）を [ADR-0273](0273-coulomb-penetrability-by-cf2-and-riccati.md)
  で足した。前は「電荷のない中性子の粒子対」だけで、荷電粒子の粒子対はエラー。4 つの版の Cl-35 と Ca-40 の 5 評価が止まっていたため（#0264）
