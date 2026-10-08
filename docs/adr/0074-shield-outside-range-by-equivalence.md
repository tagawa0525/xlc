---
status: accepted
date: 2026-10-07
supersedes: [ADR-0023]
superseded-by:
issues: ["#0001", "#0006"]
---

# ADR-0074: 減速計算の範囲の外は、ライブラリの σ0 の表を等価原理の σ0 で引き、燃料だけを遮蔽する

## 背景

共鳴の自己遮蔽は、はじめ全範囲を等価原理で扱っていた（[ADR-0023](0023-shield-by-bondarenko-and-equivalence.md)）。分離共鳴の範囲は、ピンごとの減速計算に置き換えた
（[ADR-0068](0068-solve-pinwise-slowing-down.md)）。その範囲（[ADR-0069](0069-cap-cell-slowing-down-at-20-kev.md)）の外にも共鳴がある（分離共鳴の上端より上の非分離共鳴、4 eV より下の共鳴など）。
この範囲の共鳴は狭く、NR 近似と等価原理がよく成り立つ。

## 検討した案

計画には、等価原理以外の方式を比べた記録はない（範囲の中の方式の比較は [ADR-0068](0068-solve-pinwise-slowing-down.md)）。Dancoff 係数の求め方は次を比べた。

- ピンセルの衝突確率法（全体計画の最初の版）: 無限に並んだ同じピンを仮定するので、案内管の隣や集合体の端での
  ピンごとの違いが出ない
- 燃料を黒体にした集合体全体の MOC（採用。集合体の MOC はすでにあるため）

## 決定と理由

UFG の範囲の外は、ライブラリの σ0 ごとの群定数を、等価原理の σ0 で補間して使う（[ADR-0023](0023-shield-by-bondarenko-and-equivalence.md) の方式のうち、範囲の外で
有効な部分をここにまとめ直した）。

- Dancoff 係数は、燃料を黒体にした固定源の MOC を集合体全体で解き、ピンごと・群ごとに求める。被覆管と減速材には
  その群の全断面積（[ADR-0087](0087-use-dilute-totals-in-shadows.md)）と、無限に広がれば中性子束が 1 になる一様な等方線源（散乱なし）を置く。孤立した
  ピンへの入射の流れ（周長 × 1/4）に対する、燃料の吸収率の不足分を C とする
- 燃料の脱出確率を、Dancoff 係数 C を入れた Carlvik の2項有理近似 P(x) = Σ_i β_i x / (x + α_i)（x = Σt l̄）とする。
  A = (1 − C)/C から α₁,₂ = [(5A + 6) ∓ √(A² + 36A + 36)] / (2(A + 1))、β₁ = ((4A + 6)/(A + 1) − α₁)/(α₂ − α₁)、
  β₂ = 1 − β₁（C = 0 で Carlvik の α = 2, 3、β = 2, −1）。核種 k の背景断面積は
  σ0,i = (Σ_{j≠k} N_j σt,j + α_i / l̄) / N_k（l̄ はペレットの平均弦長 2r）。2つの σ0 での値を合わせる重みは [ADR-0089](0089-weight-two-terms-with-potential-scattering.md)
- σ0 の格子点の間は ln σ0 に線形に補間し、格子の外はエラーにする
- 核種の間の共鳴の干渉は、同じ物質の他の核種の遮蔽した全断面積を σ0 に入れ、収束するまで反復して近似する（Bondarenko
  の方法の反復。MPACT も同じ。MPACT Theory Manual v4.3, ORNL/SPR-2021/2330, 9.3 節）。UFG の範囲の中は、ポイントワイズの断面積で
  干渉をそのまま扱う
- 等価なセルの中で等価原理で遮蔽するのは燃料だけで（等価なセルの外の構造材の塊は [ADR-0086](0086-shield-structural-lumps.md) で別に遮蔽する）、等価な
  セルの被覆管と減速材は範囲の外の群で無限希釈のままとし、SPH 因子も掛けない。減速材は分離共鳴を持つ核種を含まず
  （O-16 の共鳴は MeV の領域の幅の広いもので、ライブラリでは分離共鳴として扱わない）、被覆管の Zr の分離共鳴（上端
  26〜200 keV）で範囲（4 eV〜15.03 keV）の外に残るのは上の一部だけである。SPH 因子は UFG の中性子束を参照にして求める
  ので、参照のない範囲の外には掛けない
- 1本のピンの中では σ0 を一様とする。ペレット内の半径方向の分布（リム効果）は燃焼を扱う M4 で検討する

1B の OpenMC と比べて、15〜110 keV の群の被覆管の吸収断面積は +2〜8% 大きいが、被覆管の吸収の約 6% の範囲なので
被覆管全体の吸収断面積には効かない（0.009480 /cm、OpenMC 0.009490 /cm。examples/vera/README.md、
docs/plans/m2-assembly.md、256bc5c、2026-10-08）。

## 確認

`interpolates_shielded_cross_sections_in_ln_sigma0`、`rejects_sigma0_outside_the_table`、`isolated_rod_gives_carlvik_coefficients`、
`coefficients_are_ordered_and_betas_sum_to_one`、`black_lattice_has_no_escape`、`isolated_pin_is_not_shadowed`、
`shields_fuel_and_keeps_moderator_dilute`。被覆管の吸収の OpenMC との比較（上記）。

## 結果

- 4 eV より下の共鳴の遮蔽の誤差は #0001。構造材の塊への Carlvik の係数の適用は #0006
