---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
---

# ADR-0041: 群の間の移行は、反応ごとの運動学で作る

## 背景

散乱行列（群の間の移行）を、反応の種類ごとにどう作るかを決める必要があった。

## 検討した案

計画には、ほかの案を比べた記録はない。

## 決定と理由

- 弾性散乱は、4 eV 以上では静止した標的との重心系等方散乱、4 eV 未満では自由ガス模型とする
- 離散準位の非弾性散乱は、重心系等方の運動学で扱う
- 連続的な中性子放出反応（MT91, (n,2n) など）は、MF6 LAW=1 のエネルギー分布と多重度を使う。角度と座標系
  （重心系か実験室系か）は無視する

## 確認

`hydrogen_elastic_matrix`、`discrete_inelastic_uses_q_value`、`free_gas_rows_sum_to_the_broadened_elastic_cross_section`、
`mf6_multiplicity_scales_the_row`。NJOY GROUPR との比較（examples/m1/README.md）。

## 結果

- H in H2O は 4 eV 未満で S(α,β) を使う（[ADR-0052](0052-scatter-h-in-h2o-by-sab.md)）
- 散乱行列の Legendre モーメントは [ADR-0076](0076-store-legendre-scattering-moments.md)、MF4 の角度分布は [ADR-0081](0081-use-mf4-angular-distributions.md)、MF6 の2体の法則は [ADR-0056](0056-treat-mf6-law2-law3-as-two-body.md)
