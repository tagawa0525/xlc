---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
---

# ADR-0056: MF6 の LAW=2, 3 の離散準位の非弾性散乱は、Q 値による2体の運動学で扱う

## 背景

Cr と Sn-120 は、離散準位の非弾性散乱（MT51〜90）の出射分布を MF6 の LAW=2, 3 で持つ（Cr は LAW=3、Sn-120 は LAW=2）。

## 検討した案

計画には、ほかの案を比べた記録はない。

## 決定と理由

MF3 だけの離散準位と同じく（MF4 の離散準位と同じく）、Q 値による2体の運動学で扱う。LAW=3 は重心系で等方、LAW=2 は
重心系の角度分布を使う（[ADR-0076](0076-store-legendre-scattering-moments.md)）。

## 確認

`two_body_mf6_level_uses_q_value`、`two_body_neutron_law_is_marked`、`isotropic_two_body_neutron_law_is_marked`、
`lab_frame_two_body_angles_are_error`。

## 結果

MF6 の LAW=1 と LAW=2 の角度のデータを P1 以上に入れるかは [ADR-0082](0082-leave-out-mf6-angular-data.md)。
