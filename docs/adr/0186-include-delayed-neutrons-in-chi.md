---
status: accepted
date: 2026-10-22
requires: [ADR-0039]
supersedes: [ADR-0042]
superseded-by:
issues: ["#0133", "#0148", "#0149"]
---

# ADR-0186: χ に遅発中性子を含め、即発と遅発の分布を ν_p と ν_d で重み付けして 1 本にする

## 背景

[ADR-0042](0042-collapse-chi-by-fission-source.md) は、ν を MF1 MT452（即発と遅発の和）、核分裂スペクトル χ を MF5 MT18 から取り、
χ を核分裂源（νσf と重み関数の積）で重み付けして 1 本のベクトルにすると決めた。ENDF/B-VIII.0 の MF5 MT18 は即発中性子の
分布で、遅発中性子（MF1 MT455）には族ごとの別の分布（MF5 MT455）がある。ADR-0042 の χ は、遅発中性子も即発の分布で生まれる
ものとして数えていた。k の固有値の計算では、遅発中性子も核分裂源に入り、それぞれの族の分布で生まれる。

1C のピンセル（新しい燃料。#0133）で、xlc は U-238 の高速の核分裂を生成 1 あたり OpenMC より +3.17% 多く数え、χ は OpenMC の
放出のスペクトルより硬かった（1.353 MeV より上の割合が U-235 で 0.5763 対 0.5721、U-238 で 0.5763 対 0.5681）。OpenMC の
ライブラリ（examples/depletion/openmc-library、ENDF/B-VIII.0 を NJOY で処理したもの）の U-235 の即発と遅発の分布から、熱の
核分裂（0.0253 eV）で生まれる中性子の割合を数えると（2026-10-22）:

| 分布                    | ν       | 1.353 MeV より上 | 10 MeV より上 |
| ----------------------- | ------- | ---------------- | ------------- |
| 即発（MF5 MT18）        | 2.414   | 0.57627          | 0.00139       |
| 遅発（MF5 MT455、6 族） | 0.01585 | 0.04346          | 0             |
| 和（遅発の割合 0.652%） | 2.42985 | 0.57279          | 0.00138       |

即発の割合は xlc の χ（0.57628）と同じで、和から 10 MeV より上を除いて割合を取り直すと 0.57220 になり、OpenMC の tally
（10 MeV より上を数えない）の 0.5721 と合う。U-238（2 MeV の入射）では遅発の割合が 1.68% で、即発の 0.57555 に対し和は
0.56677。χ の差は、遅発中性子を即発の分布で数えたことで説明がつく。群 1（6.07〜10 MeV）の +6.4%（U-235）、+8.4%（U-238）は、
xlc の χ が 10 MeV より上に生まれる中性子を群 1 に入れ（`outgoing_group_integrals`。NJOY GROUPR と同じ）、OpenMC の tally が
数えていなかったことによる比較の方法の差だった（#0133 の比較のスクリプトを直した）。

xlc のライブラリが使う ENDF/B-VIII.0 の核分裂性の核種（U-234、U-235、U-236、U-238、Np-237、Np-239、Pu-238〜Pu-242、
Am-241、Am-242m、Am-243、Cm-242、Cm-244）は、どれも MF1 MT452、MT455、MT456 と MF5 MT18 を持つ。MF5 MT455 は Np-239、
Am-242m、Cm-244 だけにない。OpenMC（NJOY の ACE から読んだもの）は、この 3 核種で遅発中性子を作らず、ν を即発（MT456）
だけにする。

## 検討した案

- 案 C1: ADR-0042 のまま、MT452 の全数を MF5 MT18 の分布で生む。遅発中性子を硬い分布で数え、χ の 1.353 MeV より上を
  U-235 で 0.6%、U-238 で約 1.5% 多くする
- 案 C2: 即発（MF1 MT456 の ν_p と MF5 MT18）と遅発（MF1 MT455 の ν_d と MF5 MT455。族の割合は MF5 の p_k(E)）を、それぞれの
  ν σf w で重み付けして 1 本の χ にする
- 案 C3: 即発と遅発の χ を別々にライブラリに持ち、格子計算で混ぜる。動特性の分岐テーブル（[ADR-0033](0033-store-kinetics-and-cooling-data-in-m4-tables.md)）には
  別々の値が要るが、定常の格子計算は和だけを使う

MF5 MT455 のない評価（Np-239、Am-242m、Cm-244）の扱い:

- 案 D1: エラーにする。燃焼のチェーンの核種のライブラリが作れなくなる
- 案 D2: OpenMC と同じく、遅発中性子を数えない（ν を MT456 にする）。評価が与える遅発中性子の数を捨てる
- 案 D3: 遅発中性子にも MF5 MT18 の分布を使い（案 C1 と同じ数え方）、warn を出す。評価にない分布を作らない

## 決定と理由

案 C2 と案 D3 にする。ν は MF1 MT452 のまま（ADR-0042 と同じ）。χ は、核分裂で生まれる中性子の成分ごとに ν_c(E)σf(E)w(E)
（w は [ADR-0039](0039-weight-with-njoy-iwt4.md) の重み関数）で重み付けし、和が 1 になるように正規化した 1 本のベクトルにする。

- 成分は、MF5 MT455 があれば即発（MT456、MF5 MT18）と遅発（MT455、MF5 MT455）、遅発中性子のデータがなければ全数（MT452、
  MF5 MT18）。MT455 はあるが MF5 MT455 がない評価は全数を MF5 MT18 で数え、warn を出す。MF5 MT455 があるのに MF1 MT455 か
  MT456 がない評価はエラーにする（評価が欠けている）
- 遅発中性子の分布は評価にあり、k の固有値の計算で遅発中性子はその分布で生まれる。案 C2 は評価のとおりに数え、参照解の
  OpenMC とも同じ数え方になる
- 案 C3 は、今の使い方（定常）に要らない形式を先に作る。別々の値は ADR-0033 の項目を作るときに足す
- 案 D3 は、評価が与えるもの（ν の全数と、MF5 MT18 の分布）だけを使う。D2 は ν を評価の全数から変える

## 帰結

- 核分裂性の核種の χ が軟らかくなり、ライブラリを作り直す。examples/depletion/depletion.msgpack は作り直した（1C の k が −34〜−36 pcm
  動く）。examples/vera/vera.msgpack と examples/vera/endfb70 のライブラリも作り直し、VERA の格子と 3A、3B の記録を更新した
  （5-2D の炉心の記録と、同じ vera.msgpack を使う examples/colorset と examples/split は、まだ解き直していない）。examples/depletion/p2 の
  ライブラリと 2A の記録は、まだ作り直していない
- Np-239、Am-242m、Cm-244 は、ν が OpenMC と遅発中性子の分（Np-239 で 0.8%）だけ違う
- NJOY GROUPR との比較（`compare_njoy_gendf` の例）は、GROUPR の MF6 MT18 が MF5 MT18 の分布なので、即発の分布だけで比較する
- ENDF/B-VIII.0 の MF5 MT455 は族ごとに LF=5（一般の蒸発スペクトル。θ = 1 で g(x) がヒストグラムの表）なので、MF5 の読み込みと
  群への振り分けに LF=5 を足した
- 1C のピンセルの U-238 の高速の核分裂の差は、χ を直しても +2.08% 残る。高速の群の群定数（散乱の移行と U-238 の群定数）にあり、
  候補は群の中の重み（確かめていない。#0148）

## 確認

- 1C のピンセルの χ（`examples/depletion/issue-0133/openmc_chi.py`。OpenMC は 10 MeV より上を群 1 に数えた tally、Slurm 448）: 1.353 MeV
  より上の割合は U-235 で 0.57280（OpenMC 0.57268）、U-238 で 0.56840（0.56887）。ライブラリは examples/depletion/depletion.msgpack
  （md5 ea77045e。Slurm 452 で作り、レビューのあとの 457 でも同じ）
- 1C のピンセルの新しい燃料で、U-238 の核分裂（生成 1 あたり）の OpenMC との差が +2.74% から +2.08%、U-235 の核分裂が −0.23% から −0.17%、
  k の差が +2 から −34 pcm（OpenMC の σ 18 pcm）。1C の燃焼で U-235 の消えた量の差が −0.22〜−0.24% から −0.17〜−0.19%（1 リング）。
  数値の出どころは #0133 の経緯と examples/depletion/verification.md の「格子計算の燃焼」
- `each_part_of_the_neutrons_enters_chi_by_its_nubar`、`delayed_neutrons_take_their_own_spectrum_beside_the_prompt_ones`、
  `delayed_nubar_without_its_spectrum_leaves_all_neutrons_to_mt18`、`delayed_spectrum_without_prompt_nubar_is_an_error`、
  `delayed_spectrum_without_delayed_nubar_is_an_error`、`reads_delayed_nubar_after_the_decay_constants`、
  `reads_general_evaporation`、`general_evaporation_splits_the_tabulated_shape_by_group`、
  `general_evaporation_scales_the_outgoing_energy_by_theta`、`general_evaporation_is_cut_at_e_minus_u`
