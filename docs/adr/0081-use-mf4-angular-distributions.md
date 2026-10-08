---
status: accepted
date: 2026-10-08
supersedes: []
superseded-by:
issues: ["#0009"]
---

# ADR-0081: MF4 の角度分布を読み、散乱行列の P1 以上に入れる

## 背景

M2 では2体の反応を重心系で等方として運動学の分だけを入れていた（[ADR-0076](0076-store-legendre-scattering-moments.md)）。これは k の差（断面積の側）の候補の
一つで、M3 の新しい核種（C、Si、Ni、Mn など）の弾性散乱も強い非等方を持つ。

## 検討した案

計画には、ほかの案を比べた記録はない。

## 決定と理由

MF4 の角度分布（重心系・実験室系、Legendre・表・混合）を xs で読み、散乱行列の P1 以上に入れる。2体の反応に加え、
エネルギー分布を MF5 で持つ連続的な放出（MT16、17、22、28、91 など。Ag、Cd、In、P、C で使われる）も、MF4 の角度分布で
P1 以上を入れる。

## 確認

`elastic_angles_move_the_outgoing_energies`、`tabulated_moments_are_normalized`、`legendre_coefficients_are_interpolated_between_energies`、
`mixed_distribution_switches_at_the_last_legendre_energy`、`continuum_lab_angles_scale_the_spectrum`。NJOY GROUPR（Legendre の次数 3）と
比べ、O-16 と H-1 の弾性散乱と O-16 の第1準位はどの次数も P0 の行の和に対して 1.3e-4 以内で一致する
（examples/vera/README.md）。

## 結果

U-238 の弾性散乱だけ、最も高い群で NJOY と 2e-3 違う（#0009）。MF6 の角度のデータは [ADR-0082](0082-leave-out-mf6-angular-data.md)。
