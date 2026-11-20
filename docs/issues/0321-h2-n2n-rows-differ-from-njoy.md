---
type: Issue
status: open
labels: [xs]
opened: 2026-11-20
closed:
milestone:
blocked-by: []
resolution:
---

# H-2 の (n,2n) の散乱行列のしきい値の群の出射の分け方が NJOY と 6.95e-1 違う原因が候補のまま

## 問題

- ENDF/B-VIII.0 の H-2 の MT16（MF6 の LAW=6、3 体の位相空間分布）の P0 の散乱行列が、WIMS 69 群のしきい値を含む 2 群（2.231〜3.679 MeV）の行で、
  NJOY2016 GROUPR と 6.95e-1（行の差の絶対値の和 / P0 の行の和）違う。行の和（断面積）は 6.4e-4 で一致する（examples/depletion/verification.md の
  「MF6 の LAW=6、7 の散乱行列（NJOY との比較）」）。群の番号は `compare_njoy_gendf` の 0 始まり（WIMS の 1 始まりの番号は 1 を足したもの）
- 確かめたこと:
  - xs の群ごとの割合は、入射エネルギーごとに重心系の細かい中点則と 2e-7 で一致し（テスト `phase_space_spreads_by_the_centre_of_mass_motion`）、
    しきい値の群の行は入射エネルギーの細かい中点則と 1e-5 で一致する（`phase_space_rows_follow_the_kinks_above_threshold`）
  - 2 群の入射では、実験室系の出射エネルギーは (√Ec − √E_max)² ≥ 63 keV なのに、GROUPR はその行の 1.6% を 40.85 keV より下の 11〜14 群に置く
    （運動学で届かない）
  - Li-6 の VIII.1 の MT41（4 体）の同じ差（1.59e-1）は、GROUPR と同じ入射エネルギーの格子（1.2 倍ずつ）の実験室系の表を xs に与えると
    4.07e-3 になる（examples/depletion/issue-0311/tabulate_phase_space.py の `lab`）。H-2 では同じ表で 1.07 になり、GROUPR を再現しない
- 原因の候補（確かめていない）: GROUPR の入射エネルギーの格子と実験室系の単位基底の補間に加えて、しきい値のすぐ上の最初の格子の点
  （E_max ≤ 0 を 1 eV に置き換える `f6psp`）、`cm2lab` の適応的な点の取り方と許容誤差、`getmf6` の正規化（`yld/test`）と、上方散乱を
  自群に移す処理（`jgmax`）

関係する ADR: [ADR-0271](../adr/0271-group-mf6-law6-and-law7-by-endf-102.md)

## 完了の条件

- [ ] 残りの差の原因を GROUPR の処理の段ごとに確かめ、xs の側に誤りがないことを示すか、誤りを直す。戻す条件: H-2 の (n,2n) の出射の
  分け方が、結果（重水を含む体系など）に効くと分かったとき

## 経緯

- 2026-11-20: #0311 で LAW=6 を読み、NJOY と比べて残った（#0311 を閉じるので、ここに移した）
