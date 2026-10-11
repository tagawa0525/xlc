---
status: superseded
date: 2026-10-08
supersedes: []
superseded-by: ADR-0082
issues: []
requires: []
---

# ADR-0056: MF6 の LAW=2, 3 の離散準位の非弾性散乱は、Q 値による重心系で等方な2体の運動学で扱う

## 背景

Cr と Sn-120 は、離散準位の非弾性散乱（MT51〜90）の出射分布を MF6 の LAW=2, 3 で持つ（Cr は LAW=3、Sn-120 は LAW=2）。
M2（2026-10-06）で扱い方を決める必要があった。

## 検討した案

計画には、ほかの案を比べた記録はない。

## 決定と理由

MF3 だけの離散準位と同じく、Q 値による2体の運動学（重心系で等方）で扱う。LAW=3 は重心系で等方なので厳密で、LAW=2 の
角度分布を無視するのは、散乱を P0 で扱う M0 からの近似と同じである（docs/plans/m2-assembly.md、0d01ad8、2026-10-06）。

## 確認

`two_body_mf6_level_uses_q_value`、`two_body_neutron_law_is_marked`、`isotropic_two_body_neutron_law_is_marked`、
`lab_frame_two_body_angles_are_error`。

## 帰結

M3 で、MF6 の角度のデータを使わないことを P3 の散乱行列のもとで決め直した（[ADR-0082](0082-leave-out-mf6-angular-data.md)）。扱いは同じで、[ADR-0082](0082-leave-out-mf6-angular-data.md) が
理由を改めて示し、この ADR を置き換える。間の M2 の計画で LAW=2 の角度分布を使うことに一時したが、実装していない
（[ADR-0082](0082-leave-out-mf6-angular-data.md) の「検討した案」）。
