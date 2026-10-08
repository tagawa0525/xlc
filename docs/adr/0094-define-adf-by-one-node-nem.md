---
status: accepted
date: 2026-10-09
supersedes: [ADR-0093]
superseded-by:
issues: ["#0020", "#0026", "#0027"]
---

# ADR-0094: 非反射の格子計算の ADF は、一般化等価理論により、core と同じ NEM で各ノードを1つずつ解いた均質な中性子束から作る

## 背景

提案の時点（2026-10-08）では、ADF は外周がすべて鏡面反射の格子だけで作っていた（[ADR-0093](0093-choose-lattice-boundary-conditions.md)、`adf_needs_reflective_sides`）。均質な解が平らになり、
体積平均の中性子束を均質な面の中性子束としてよいためである。M3 では、反射体の定数（[ADR-0084](0084-make-reflector-constants-with-fuel.md)）とカラーセット
（4-2D）のように、均質な解が平らにならない格子から ADF を作る必要がある。core の NEM は ADF を入れる計画である
（[ADR-0027](0027-solve-core-with-nem.md)、[ADR-0045](0045-expand-nem-to-fourth-order.md)）が、いまはまだ使わない。`crates/core/src/nem.rs` はノードを部分中性子流でつなぎ、
`FewGroup::adf` を読まない。NEM に ADF を入れるのは M3 の進め方の 6 である。core の NEM の外周の境界は、正味の流れ 0
（`Reflective`）と中性子束 0（`ZeroFlux`）だけで、格子計算の真空（入射 0）に当たる境界はない。

- 定数を作った体系では、炉心計算が参照（格子計算の MOC）のノード平均の反応率、面の正味の流れ、k を再現すること
- core の NEM と ADF の定義がそろうこと
- 段階のクレートの独立（[ADR-0002](0002-keep-stage-crates-independent.md)）と、共通の計算の切り出し（[ADR-0003](0003-extract-shared-computation-on-second-use.md)）
- 原理で決まり、手で合わせる係数を持たないこと

## 検討した案

- 体積平均の中性子束で割る（いまの反射境界の定義をそのまま使う）: 非反射の格子では均質な解が平らでないので、一般化
  等価理論（GET）の ADF にならない。反射体（核分裂がない）では意味を持たない
- 均質化した定数で、同じ体系を細かいメッシュの拡散計算で解き、その面の中性子束を使う（GET の均質な解を、空間の離散化
  の誤差なく求める）: core の NEM と解き方が違うので、定数を作った体系でも、NEM の誤差の分だけ参照を再現しない。
  拡散のソルバを lattice に別に持つことになる
- 各ノードを、そのノードの均質化した定数、参照の k、参照（MOC）の面の正味の流れを境界条件として、core と同じ NEM で
  1ノードずつ解き、面の均質な中性子束を求める（採用）。Smith の GET の定義（K. S. Smith, "Assembly Homogenization
  Techniques for Light Water Reactor Analysis," Prog. Nucl. Energy 17(3), 303–335 (1986)。nTRACER-RENUS も ADF の出典に
  挙げる）。MPACT の節点の定数の生成は同じ方式で、1ノードの NEM（4次）または SENM の核を、輸送計算の面の正味の流れと
  ノード平均の定数で解き、ADF = 不均質な面の中性子束 / 均質な面の中性子束とする（MPACT 4.4 Theory Manual、
  ORNL/SPR-2023/3090、13.1.4 節）。オープンソースの
  格子コード Scarabee も、反射体の ADF を NEM の1ノードの解で作る（https://scarabee.readthedocs.io/en/latest/_modules/scarabee/reseau/reflector.html ）
- 同じ量を core で計算する（lattice は不均質な面の中性子束と正味の流れを出し、core が自分の NEM で ADF にする）:
  NEM のコードを動かさずに済むが、段階の間のファイルが ADF でなく途中の量になり、状態の分岐のテーブルで何を
  内挿するかがあいまいになる
- SPH 因子（DRAGON/DONJON が主に使う。A. Hébert, Applied Reactor Physics, 3rd ed., Presses inter Polytechnique (2020)、
  4 章。未確認）: 面の中性子束の不連続ではなく断面積を補正するので、ADF を入れる NEM（[ADR-0027](0027-solve-core-with-nem.md)）と形が合わない
- 真空の境界の近くで ADF が大きく外れたら、横方向の漏れと除去断面積に補正を足し、ADF が 0.8〜1.2 に入るまで Newton 法で
  合わせる（MPACT の adaptive discontinuity factors、同 13.2.3 節）: 目標の範囲が手で決めた値で、原理で決まらない

角の不連続因子（CDF）も比べた。PARCS は、単一集合体の計算での角の中性子束 / 平均の中性子束を CDF とし、ピン出力の
再構成で角の中性子束の連続に使う（PARCS v3.0 Theory Manual、5.1.2.3 節、https://www.nrc.gov/docs/ML1016/ML101610117.pdf ）。
STREAM も CDF をピン出力の再構成のために作る（S. Choi et al., "Recent Development Status of Neutron
Transport Code STREAM," Trans. KNS Spring Meeting (2019)、2.5 節）。MPACT は角の
中性子束を求められないとし、角に接する細かいメッシュの中性子束の体積平均で近似する（同 13.1.3 節）。

## 決定と理由

ADF は、GET の定義（面の不均質な中性子束 / 同じ面の均質な中性子束）で作り、均質な中性子束は、core と同じ NEM で各ノード
を1つずつ解いて求める。core が同じ NEM を使えば、定数を作った体系では参照のノード平均と面の正味の流れを NEM の近似に
よらず再現するためである。細かいメッシュの拡散で解く案は GET の元の定義に近いが、core の NEM の誤差が残る。

- 1ノードの問題は、ノード平均の定数、参照の k、参照のノード平均の中性子束、MOC の辺ごとの正味の流れ（MOC のスイープが
  CMFD のために数える粗メッシュの面の正味の流れを、辺にわたって足して辺の長さで割ったもの。MOC の解の `currents`）を
  与え、方向ごとに core と同じ NEM の式（4次の展開、f1, f2 で重み付けした残差、2次の横方向の漏れ）で解く。ノード平均は
  参照の値に固定する。参照の k と定数ではノードのバランス式を参照の値が満たすので、バランス式から求め直さない
  （k が k∞ と違うと、ノード平均をバランス式で更新する反復は k∞/k の比で増えて発散する）。鏡面反射の辺の正味の流れも、
  0 とせず MOC の値を使う（収束の許容差の程度で 0）
- core の NEM の計算は、入射の部分中性子流から出射を求める。1ノードの問題の面の条件（正味の流れが参照の値）は、
  入射の形で β = 1、外からの入射 −J と書けるので、1ノードの問題は core と同じ `xlc_numeric::nem::outflow` で解き、
  ノード平均を参照の値に固定する。`outflow` は、ノードの中で面の条件、a1, a2 と、3次・4次の係数 a3, a4 の残差の式を
  同時に解き、群の間の結合と、隣のノード（1ノードの問題では真空の辺の横方向の漏れ）だけを前の反復の値とする。
  a3, a4 を前の反復の a1, a2 から求めると、両方の面の正味の流れが決まるとき反復が発散する（正味の流れが a1 + a3/2 と
  3a2 + a4/5 を決め、D/h² ≪ Σr のノードでは a3 ≈ 10 a1、a4 ≈ 35 a2 なので、誤差が反復ごとに約 −5 倍と −2.3 倍になる。
  core でも両側が反射の1ノードの方向で起きていた。#0025）
- 横方向の漏れの2次の形の外側は、格子計算の境界に当たる core の境界の扱いを使う: 鏡面反射は core の `Reflective`
  （外に鏡像のノード）、真空は core の真空の面（面の漏れ = ノードの漏れ × 面の中性子束 / ノード平均。
  [ADR-0101](0101-impose-transport-current-ratio-at-vacuum-faces.md)）。中性子束 0 の境界（`ZeroFlux`）は真空と違う条件なので、
  等価性のためには使わない。アルベドの辺は core に対応する境界がないので、ADF を作る格子ではエラーにする。格子は
  全体を1つのノードとして均質化するので、格子の中の隣のノードはない（カラーセットをノードごとに均質化するときに、
  隣のノードの参照の漏れから2次式を作る）
- 1ノードの問題は 2次元（x, y）で、z の漏れは 0（2次元の格子計算は軸方向に一様）
- 外周がすべて鏡面反射の単一の格子では、正味の流れが 0 で均質な解は平らになり、いまの定義（体積平均の中性子束で割る）に
  一致する
- NEM の1ノードの計算は lattice と core の両方で使うので、共通の計算のクレート（`xlc-numeric`）に移す（[ADR-0003](0003-extract-shared-computation-on-second-use.md)。lattice は
  core に依存しない。[ADR-0002](0002-keep-stage-crates-independent.md)）。1ノードの問題の解き方（上）も同じクレートに置き、
  core の等価性のテストからも使う
- ADF は NEM の展開の次数、横方向の漏れの近似、拡散係数の定義に結びつく。これらのどれかを変えたら、テーブルを作り直す
- CDF も出す。単一集合体（外周が鏡面反射）の格子では、均質な解が平らなので、角の不均質な中性子束 / ノード平均の
  中性子束とする。角の不均質な中性子束は、集合体の外周の余白（薄い行と列。[ADR-0083](0083-generalize-to-rectilinear-lattice.md)）の角の小さな領域（VERA の
  17×17 では一辺 0.04 cm）の中性子束とする。平均自由行程（〜1 cm）に比べて小さく、点の値の近似の誤差は小さいと
  見込む。CDF はピン出力の再構成でだけ使う
- 均質な解が平らでない格子（反射体、カラーセット）の CDF は、角の均質な中性子束をノードの中の 2D の均質な解から求める
  必要があり、決めていない（#0020）。それまで、CDF は外周がすべて鏡面反射の格子でだけ作り、ほかの格子ではエラーにする
- MPACT のように ADF を範囲に収める補正はしない。1ノードの解の面の中性子束が正でないなど、ADF が作れないときはエラーに
  する

### ADR-0093 との関係

[ADR-0093](0093-choose-lattice-boundary-conditions.md) は「ADF は鏡面反射の辺だけで作る」と決めていた（`adf_needs_reflective_sides`）。この ADR は ADR-0093 を
置き換える。ADR-0093 の次の決定はそのまま有効で、ここに書き直して含む（ADR の運用。docs/adr/README.md）。

- 外周の境界条件は、辺ごとに鏡面反射、真空、アルベドから選ぶ。鏡面反射は cyclic tracking（[ADR-0060](0060-track-whole-assembly-cyclically.md)）でトラックを
  つなぎ、真空は入射を 0 とする
- アルベドの辺は鏡面反射と同じくトラックをつなぎ、出ていく角度中性子束を α（戻る分）と 1 − α（漏れ）に分ける。
  0 ≤ α ≤ 1 を確かめ、α = 1 と 0 は鏡面反射と真空として持つ
- 入力では、辺ごとに `"reflective"`、`"vacuum"`、`{ albedo = α }` のどれかを書く

変えるのは ADF の作り方だけで、鏡面反射と真空の辺の格子で作る（上のとおり。アルベドの辺はエラー）。

## 確認

- 1ノードの問題（`xlc_numeric::nem::surface_fluxes`）: 正味の流れが 0 で外周が反射なら面の中性子束はノード平均に
  相対 1e-14 で一致（`zero_currents_give_the_node_average`）、外へ流れ出る面で下がる（`outward_current_lowers_the_face_flux`）、
  入力の誤り（`invalid_one_node_problems_are_errors`）と正でない均質な面の中性子束（`non_positive_surface_flux_is_error`）は
  エラー
- 等価性（core の NEM と、細かいノードの NEM の解を参照にした GET の ADF）: 細かいノード（2.5 cm × 16、3種類の物質）の
  平板を 20 cm の粗いノード2つにまとめ、粗いノードの間の面と真空の面を通るもの（`get_adfs_reproduce_a_two_node_slab`）で
  k の相対差 −1.6e-9（ADF 1 では +6.0e-3）、5 cm × 4 × 4 の 2D を 20 cm の1ノードにまとめ、2つの真空の辺を持つもの
  （`get_adfs_reproduce_a_node_with_vacuum_sides`）で −1.7e-8（ADF 1 では −3.1e-2）。どちらもノード平均の中性子束は
  相対 1e-6 以内（参照の真空の面も core の Marshak の条件なので、P1 の関係による差はない）
- lattice: 外周が鏡面反射なら、辺の上の中性子束 / 体積平均（ADR-0093 までの定義）と相対 1e-6 で一致
  （`reflective_adfs_divide_by_the_average_flux`。`homogeneous_cell_has_unit_adf`、`fuel_pin_raises_the_thermal_adf` も
  引き継ぐ）、真空の辺のある格子は1ノードの問題の均質な面の中性子束で割り、辺の対応は（左、右、上、下）=（-x、+x、+y、
  -y）（`open_lattice_adfs_divide_by_the_one_node_flux`）、アルベドの辺はエラー（`adf_rejects_albedo_sides`）。MOC の
  辺の正味の流れは、真空の辺で漏れと相対 1e-6 で合い、鏡面反射の辺で 0（`side_currents_add_up_to_the_leakage`）
- MOC との比較（examples/vera/verification.md の「真空の辺のある格子の ADF と core の1ノードの NEM」、6c2b970、2026-10-09）:
  VERA 2B の燃料ピンの n × n の格子（右の辺が真空）を core の1ノードで解くと、k は MOC から n = 3、5、9 で +3724、
  +3803、+1510 pcm ずれる（平らな ADF では +25170〜+37875 pcm）。真空の辺の ADF を 2J / 均質な面の中性子束に置き換える
  と 10⁻⁶ で一致するので、差はすべて core の真空の面の条件の P1 の関係 φ_het = 2J による（MOC では φ_het/(2J) =
  0.87〜0.95。#0026）。core の真空の面に格子計算の比 r = J/φ_het を課すと（[ADR-0101](0101-impose-transport-current-ratio-at-vacuum-faces.md)）、
  k は MOC と 0.005 pcm 以内で一致する（同じファイルの「真空の辺の比 J/φ を入れた core の1ノードの NEM」、43ba22d、
  2026-10-09）。17×17 の集合体では熱群の均質な面の中性子束が負になり、ADF を作れない（#0027）
- 4-2D（カラーセット）を単一集合体の ADF を入れた core の NEM で解き、同じ多群定数の xlc の MOC と比べる（判定の閾値の
  案: k ±50 pcm、集合体出力の差の最大 1%。ユーザーの決定が要る）のは、カラーセットの入力と、格子のノードごとの均質化が
  できてから。ADF そのものの OpenMC との比較は #0010

## 結果

- 真空の辺のある格子の等価性は、core の真空の面の条件を格子計算の比 r = J/φ_het にして保つ（ADR-0101。ADR-0099 の
  Marshak の条件では、P1 の関係 φ_het = 2J の分だけ崩れ、上の比較で最大 3800 pcm）。燃料が真空に面する 17×17 の
  集合体では ADF を作れない。扱いは #0027 で決める
- 1ノードの問題と core は同じ `outflow` を使う（a3, a4 をノードの中で a1, a2 と同時に解く。#0025）
- 格子は全体を1つのノードとして均質化する。カラーセットや反射体の格子をノードごとに均質化するときに、内側の面の
  正味の流れと不均質な中性子束、隣のノードの漏れを使う横方向の漏れの2次式を加える
- CDF はまだ出していない（少数群定数の形式にない）。外周が鏡面反射の格子の CDF は、ピン出力の再構成（M3 の進め方の 7）
  とあわせて加える。均質な解が平らでない格子の CDF は #0020
- 反射体の ADF、CDF の使い方（ピン出力の再構成）、分岐のテーブルに入れる量は、この ADR を前提に別に決める
- core の NEM の ADF の面の条件と真空の外周の境界は [ADR-0101](0101-impose-transport-current-ratio-at-vacuum-faces.md)（[ADR-0099](0099-couple-core-nodes-by-get-with-vacuum-boundary.md) を置き換え）
- 単一集合体の拡散係数は B1 の臨界スペクトルから、真空の辺のある格子では outflow の輸送断面積から作る
  （[ADR-0098](0098-derive-diffusion-coefficients-from-b1-and-currents.md)）。ADF は定数を作った体系では D の誤差も吸収するが、
  別の体系では吸収しない
- 細部は、xlc の格子計算（MOC）の結果との比較を通じて段階的に詰め、新しい ADR か issue に記録する（2026-10-08 の
  ユーザーの方針）
