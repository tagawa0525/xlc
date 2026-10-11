---
status: superseded
date: 2026-10-08
supersedes: []
superseded-by: ADR-0081
issues: []
requires: []
---

# ADR-0041: 群の間の移行は、反応ごとの運動学で作る

## 背景

M0 の縦切り（2026-10-06）で、散乱行列（群の間の移行）を、反応の種類ごとにどう作るかを決める必要があった。

## 検討した案

計画には、ほかの案を比べた記録はない。

## 決定と理由

- 弾性散乱は、4 eV 以上では静止した標的との重心系等方散乱、4 eV 未満では自由ガス模型とする
- 離散準位の非弾性散乱は、重心系等方の運動学で扱う
- 連続的な中性子放出反応（MT91, (n,2n) など）は、MF6 LAW=1 のエネルギー分布と多重度を使う。角度と座標系
  （重心系か実験室系か）は無視する

## 確認

`hydrogen_elastic_matrix`、`discrete_inelastic_uses_q_value`、`free_gas_rows_sum_to_the_broadened_elastic_cross_section`、
`mf6_multiplicity_scales_the_row`。NJOY GROUPR との比較（examples/m1/verification.md）。

## 帰結

M3 で、2体の反応と MF5 の連続的な放出に MF4 の角度分布を入れる決定に置き換えた（[ADR-0081](0081-use-mf4-angular-distributions.md)）。まだ有効な部分
（4 eV を境にした弾性散乱の扱い、離散準位の Q 値の2体反応、MF6 LAW=1 のエネルギー分布と多重度）は [ADR-0081](0081-use-mf4-angular-distributions.md) に
まとめ直した。重心系の MF6 LAW=1 の出射エネルギーは、計画に書かないまま、重心系で等方として実験室系に変換する
実装にしていた（be3d8f8、2026-10-06）。
