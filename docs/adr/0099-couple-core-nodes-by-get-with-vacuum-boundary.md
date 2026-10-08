---
status: accepted
date: 2026-10-09
supersedes: []
superseded-by:
issues: []
---

# ADR-0099: 炉心の NEM は ADF を GET の面の条件として部分中性子流に入れ、外周に Marshak の真空の境界を加える

## 背景

炉心の NEM で ADF を使うことは [ADR-0027](0027-solve-core-with-nem.md) で決めている。ここでは、その面の条件をノード間の
部分中性子流の結合にどう入れるかと、外周の境界条件を決める。外周は反射と中性子束 0 だけで、中性子束 0 は外挿距離を
無視するので、反射体の外側などの真空の面を表せない。提案中の [ADR-0094](0094-define-adf-by-one-node-nem.md)（1ノードの
NEM で ADF を作る）と [ADR-0095](0095-make-reflector-constants-in-2d-with-normalized-adf.md)（反射体の定数）は、core に
入射の部分中性子流が 0 の境界と、その境界での横方向の漏れの扱いを求めている。

## 検討した案

- 真空の面: 中性子束 0 の面を外挿距離だけ外に置く。外挿距離を入力で決める必要があり、ノードの幅と面の位置がずれる
- 真空の面: 均質な中性子束 φ に Marshak の条件 J- = φ/4 − J/2 = 0 を課し、外周の面の ADF を使わない。ブランチ
  `feat/core-adf` で試した（e216c37、2026-10-08）。格子計算の真空の境界は不均質な解に課すので、外周の面の ADF が 1 で
  なければ GET の等価性が崩れる
- 真空の面: 不均質な中性子束 f φ に Marshak の条件を課す（採用）

## 決定と理由

- 面の条件は一般化等価理論（GET。K. S. Smith, Prog. Nucl. Energy 17, 303 (1986)）のとおり、ADF f をかけた不均質な
  中性子束 f φ と正味の中性子流を連続にする。自分の側の ADF を f、隣を f'、出ていく部分中性子流を o、o' とすると、
  入射は i = ((f' − f) o + 2 f' o')/(f + f')。外周の境界と同じ i = β o + 外からの入射 の形なので、ノードのバランス式と
  同時に解く既存の計算をそのまま使える
- 径方向の ADF は少数群定数の辺ごとの値（−x が左、+x が右、−y が下、+y が上）とする。軸方向の面は 1 とする。2次元の
  格子計算は軸方向に一様で、軸方向の ADF を出さないため。ADF は少数群定数の必須の値なので、欠けたり正でなかったりすれば
  エラーにする
- 外周の面ごとに `vacuum` を選べるようにし、不均質な入射の部分中性子流を 0 とする（Marshak の条件を GET の不均質な
  中性子束に課したもの）。面の ADF を f、正味の中性子流を J = o − i とすると、f φ/4 − J/2 = (f(o + i) − (o − i))/2 = 0
  から i = β o、β = (1 − f)/(1 + f)。隣の面の式で f' = 1、o' = 0 としたものに等しい。外挿距離は 2D/f になる
- 反射（J = 0、β = 1）と中性子束 0（f φ = 0 は φ = 0 と同じ、β = −1）は、GET の条件にしても f によらない
- 横方向の漏れの2次式では、真空の面での漏れを、ノードの漏れ × 面の中性子束 / ノード平均とする。中性子束 0 の面で漏れを
  0 とする既存の条件と同じく、漏れが中性子束に比例するとみなすものである

## 確認

- ADF: すべて 1 なら 47c63a7 の結果と一致（`unit_adfs_leave_results_unchanged`）、内側の面で f φ と正味の中性子流が連続
  （`adfs_make_homogeneous_flux_jump_at_faces`）、同じ断面積で ADF の違う2領域の平板の解析解
  （f_A tan(Ba) + f_B tan(Bb) = 0）と相対 1e-8 で一致（`two_region_slab_matches_analytic_jump`）、すべての ADF を群ごとに
  同じ倍率で変えても、真空の面がなければ k は変わらない（`uniform_adf_scaling_keeps_k_without_vacuum`）、x、y、z の
  各方向で辺の取り違えがない（`side_mapping_matches_two_region_analytic`。境の2ノードだけに跳びを置き、上と下を取り違える
  と k が 3e-3 ずれる）
- 真空: 不均質な入射が 0（`vacuum_face_has_no_incoming_heterogeneous_current`）、k が中性子束 0 と反射の間
  （`vacuum_k_lies_between_zero_flux_and_reflective`）、1群の裸の平板の解析解（tan(Bl/2) = f/(2DB)）と f = 0.7〜2 で
  相対 1e-8（`bare_slab_matches_marshak_analytic`、`outer_adf_at_vacuum_matches_analytic`）、各方向に分離できる裸の直方体の解析解と 8³ ノードで相対 1e-4
  （`bare_cuboid_matches_marshak_analytic`）、物質の境が真空の角にある分離できない 2D で 32 × 32 ノードの解に近づく
  （`vacuum_corner_converges_with_refinement`。4 × 4 で 9.6e-5、8 × 8 で 1.1e-5）、真空の面のノード平均が正でなければ
  エラー（`no_flux_next_to_a_vacuum_face_is_error`）、入力の `parses_vacuum_boundary`

## 結果

- ADR-0094 と ADR-0095 が求める、入射の部分中性子流が 0 の core の境界と、その境界での横方向の漏れの扱いはこの ADR の
  とおり。ADR-0094 の1ノードの NEM は、格子計算の参照の正味の流れから均質な解を作り、真空の辺の ADF を不均質な面の
  中性子束 / 均質な面の中性子束とする。core は真空の面の条件をその ADF をかけた f φ に課すので、格子計算の真空の境界
  （不均質な解の入射 0）と同じ条件になる。1ノードの NEM の横方向の漏れも、真空の面ではこの ADR と同じ扱いを使う
- アルベドの境界は core にない。必要になったら決める
