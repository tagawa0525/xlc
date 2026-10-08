---
status: accepted
date: 2026-10-08
supersedes: [ADR-0041]
superseded-by:
issues: ["#0009"]
---

# ADR-0081: 群の間の移行は、反応ごとの運動学と MF4 の角度分布で作る

## 背景

M0 では、2体の反応を重心系で等方とし、連続的な放出の角度を無視していた（[ADR-0041](0041-transfer-by-two-body-kinematics.md)）。M2 で散乱行列を P0〜P3 の
モーメントで持つようにしたが（[ADR-0076](0076-store-legendre-scattering-moments.md)）、2体の反応は重心系で等方として運動学の分だけを入れていた。これは k の差
（断面積の側）の候補の一つで、M3 の新しい核種（C、Si、Ni、Mn など）の弾性散乱も強い非等方を持つ。

## 検討した案

- 重心系で等方のまま（[ADR-0041](0041-transfer-by-two-body-kinematics.md)）: 上のとおり
- MF4 の角度分布を入れる（採用）

## 決定と理由

MF4 の角度分布（重心系・実験室系、Legendre・表・混合）を xs で読み、散乱行列に入れる。2体の反応では、角度分布は
出射エネルギーも変えるので、P0 の行にも効く。反応ごとの扱いは次のとおり（[ADR-0041](0041-transfer-by-two-body-kinematics.md) のうち有効な部分を含む）。

- 弾性散乱: 4 eV 以上は静止した標的との2体反応で、MF4 の角度分布を使う。4 eV 未満は自由ガス模型（H in H2O は
  S(α,β)。[ADR-0052](0052-scatter-h-in-h2o-by-sab.md)）。MF6 に MT2 があっても、2体反応と自由ガスで扱う
- 離散準位の非弾性散乱: Q 値による2体反応で、MF4 の角度分布を使う。角度分布がなければ重心系で等方とする。MF6 の
  LAW=2, 3 で与えられる準位は [ADR-0082](0082-leave-out-mf6-angular-data.md)
- 連続的な中性子放出反応（MT91, (n,2n) など）: MF6 LAW=1 のエネルギー分布と多重度を使う。重心系で与えられた分布は、
  重心系で等方として実験室系に変換する。MF6 の角度のデータは使わない（[ADR-0082](0082-leave-out-mf6-angular-data.md)）
- エネルギー分布を MF5 で持つ連続的な放出（MT16、17、22、28、91 など。Ag、Cd、In、P、C、Zircaloy-4 の Sn で使われ、
  Sn では MF4 が実験室系で非等方）: MF4 の実験室系の角度分布で P1 以上を入れる。出射エネルギーと角度は独立とする。
  重心系の MF4 は、等方なら P1 以上を 0 とし、非等方ならエラーにする

## 確認

`elastic_angles_move_the_outgoing_energies`、`center_of_mass_distribution_shifts_the_outgoing_energies`、
`tabulated_moments_are_normalized`、`legendre_coefficients_are_interpolated_between_energies`、
`mixed_distribution_switches_at_the_last_legendre_energy`、`continuum_lab_angles_scale_the_spectrum`、
`center_of_mass_continuum_angles_are_error`、`isotropic_center_of_mass_continuum_is_accepted`、
`elastic_uses_free_gas_even_with_mf6_data`、`discrete_inelastic_uses_q_value`、`center_of_mass_energy_spreads_in_laboratory`、
`mf6_multiplicity_scales_the_row`。NJOY GROUPR（Legendre の次数 3）と比べ、O-16 と H-1 の弾性散乱と O-16 の第1準位は
どの次数も P0 の行の和に対して 1.3e-4 以内で一致する（examples/vera/verification.md）。

## 結果

- U-238 の弾性散乱だけ、最も高い群で NJOY と 2e-3 違う（#0009）
- MF6 の角度のデータは [ADR-0082](0082-leave-out-mf6-angular-data.md)、モーメントの持ち方は [ADR-0076](0076-store-legendre-scattering-moments.md)
