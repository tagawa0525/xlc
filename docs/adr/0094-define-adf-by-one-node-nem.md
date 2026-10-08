---
status: proposed
date: 2026-10-08
supersedes: []
superseded-by:
issues: []
---

# ADR-0094: 非反射の格子計算の ADF は、一般化等価理論により、core と同じ NEM で各ノードを1つずつ解いた均質な中性子束から作る

## 背景

いまの ADF は、外周がすべて鏡面反射の格子だけで作る（[ADR-0093](0093-choose-lattice-boundary-conditions.md)、`adf_needs_reflective_sides`）。均質な解が平らになり、
体積平均の中性子束を均質な面の中性子束としてよいためである。M3 では、反射体の定数（[ADR-0084](0084-make-reflector-constants-with-fuel.md)）とカラーセット
（4-2D）のように、均質な解が平らにならない格子から ADF を作る必要がある。core は NEM に ADF を入れる
（[ADR-0027](0027-solve-core-with-nem.md)、[ADR-0045](0045-expand-nem-to-fourth-order.md)）。

- 定数を作った体系では、炉心計算が参照（格子計算の MOC）のノード平均の反応率、面の正味の中性子流、k を再現すること
- core の NEM と ADF の定義がそろうこと
- 段階のクレートの独立（[ADR-0002](0002-keep-stage-crates-independent.md)）と、共通の計算の切り出し（[ADR-0003](0003-extract-shared-computation-on-second-use.md)）
- 原理で決まり、手で合わせる係数を持たないこと

## 検討した案

- 体積平均の中性子束で割る（いまの反射境界の定義をそのまま使う）: 非反射の格子では均質な解が平らでないので、一般化
  等価理論（GET）の ADF にならない。反射体（核分裂がない）では意味を持たない
- 均質化した定数で、同じ体系を細かいメッシュの拡散計算で解き、その面の中性子束を使う（GET の均質な解を、空間の離散化
  の誤差なく求める）: core の NEM と解き方が違うので、定数を作った体系でも、NEM の誤差の分だけ参照を再現しない。
  拡散のソルバを lattice に別に持つことになる
- 各ノードを、そのノードの均質化した定数、参照の k、参照（MOC）の面の正味の中性子流を境界条件として、core と同じ NEM で
  1ノードずつ解き、面の均質な中性子束を求める（採用）。Smith の GET の定義（K. S. Smith, "Assembly Homogenization
  Techniques for Light Water Reactor Analysis," Prog. Nucl. Energy 17(3), 303–335 (1986)。nTRACER-RENUS も ADF の出典に
  挙げる）。MPACT の節点の定数の生成は同じ方式で、1ノードの NEM（4次）または SENM の核を、輸送計算の面の正味の中性子流と
  ノード平均の定数で解き、ADF = 不均質な面の中性子束 / 均質な面の中性子束とする（MPACT 4.4 Theory Manual、
  ORNL/SPR-2023/3090、13.1.4 節）。オープンソースの
  格子コード Scarabee も、反射体の ADF を NEM の1ノードの解で作る（https://scarabee.readthedocs.io/en/latest/_modules/scarabee/reseau/reflector.html ）
- 同じ量を core で計算する（lattice は不均質な面の中性子束と正味の中性子流を出し、core が自分の NEM で ADF にする）:
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
を1つずつ解いて求める。core が同じ NEM を使えば、定数を作った体系では参照のノード平均と面の正味の中性子流を NEM の近似に
よらず再現するためである。細かいメッシュの拡散で解く案は GET の元の定義に近いが、core の NEM の誤差が残る。

- 1ノードの問題は、ノード平均の定数、参照の k、MOC の面ごとの正味の中性子流（CMFD の粗メッシュの面の流れを、ノードの面に
  わたって足したもの）を与え、方向ごとに NEM（4次の展開、2次の横方向の漏れ）で解く。横方向の漏れの2次の形は、格子の中の
  隣のノードの参照の漏れから作り、格子の外周では core が炉心の外周で使うのと同じ扱いにする
- 外周がすべて鏡面反射の単一の格子では、正味の中性子流が 0 で均質な解は平らになり、いまの定義（体積平均の中性子束で割る）に
  一致する
- NEM の1ノードの計算は lattice と core の両方で使うので、共通の計算のクレート（`xlc-numeric`）に移す（[ADR-0003](0003-extract-shared-computation-on-second-use.md)。lattice は
  core に依存しない。[ADR-0002](0002-keep-stage-crates-independent.md)）
- ADF は NEM の展開の次数、横方向の漏れの近似、拡散係数の定義に結びつく。これらのどれかを変えたら、テーブルを作り直す
- CDF も出す。単一集合体（外周が鏡面反射）の格子で、角の不均質な中性子束 / ノード平均の中性子束とする。角の不均質な
  中性子束は、集合体の外周の余白（薄い行と列。[ADR-0083](0083-generalize-to-rectilinear-lattice.md)）の角の小さな領域（VERA の 17×17 では一辺 0.04 cm）
  の中性子束とする。平均自由行程（〜1 cm）に比べて小さく、点の値の近似の誤差は小さいと見込む。CDF はピン出力の再構成
  でだけ使う
- MPACT のように ADF を範囲に収める補正はしない。1ノードの解の面の中性子束が正でないなど、ADF が作れないときはエラーに
  する

## 確認

まだ実装していない（M3 の進め方の 5）。次で確かめる。

- 外周が鏡面反射の単一の格子で、いまの定義と同じ値になること（`homogeneous_cell_has_unit_adf` などを引き継ぐ）
- 定数を作った体系（2 ノード以上の格子）を core の NEM で解くと、MOC のノード平均の反応率と k を、反復の収束の許容差の
  範囲で再現すること（等価性のテスト）
- 4-2D（カラーセット）を、単一集合体の ADF を入れた core の NEM で解き、同じ多群定数の xlc の MOC で同じ問題を解いた
  結果と比べる（2段階の手法の誤差だけを測る）。あわせて OpenMC（2D）と KENO の値と比べる。判定の閾値（案）: xlc の
  MOC に対して k ±50 pcm、集合体出力の差の最大 1%（ユーザーの決定が要る）
- ADF そのものの OpenMC との比較は #0010

## 結果

- 反射体の ADF、CDF の使い方（ピン出力の再構成）、分岐のテーブルに入れる量は、この ADR を前提に別に決める
- 拡散係数は D = 1/(3Σt)（輸送補正なし。[ADR-0044](0044-collapse-to-two-groups.md)）で、見直しを [ADR-0098](0098-derive-diffusion-coefficients-from-b1-and-currents.md) で提案している。ADF は
  定数を作った体系では D の誤差も吸収するが、別の体系では吸収しない
- 1ノードの解の ADF が真空の境界の近くで大きく外れたら（MPACT が報告する症状）、原因を調べてから扱いを決める
- 細部は、xlc の格子計算（MOC）の結果との比較を通じて段階的に詰め、新しい ADR か issue に記録する（2026-10-08 の
  ユーザーの方針）
