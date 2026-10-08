---
status: accepted
date: 2026-10-08
supersedes: []
superseded-by:
issues: []
---

# ADR-0080: 565 K の H in H2O は、550 K と 573.6 K の熱中性子の散乱行列を温度で線形に内挿する

## 背景

HZP は 565 K だが、H in H2O の S(α,β) は 550 K と 573.6 K にしかない。格子計算は温度方向に補間しないとしていた
（[ADR-0047](0047-choose-temperatures-per-nuclide.md)）。

- 参照解（OpenMC）と同じ扱いにできること

## 検討した案

- S(α,β) の表そのものを温度で内挿する: NJOY にも OpenMC にもないので、参照解と同じ扱いにできない
- 群定数を温度について線形に内挿する（採用）

## 決定と理由

xs は両方の温度の熱中性子の散乱行列をライブラリに持ち、lattice が物質の温度で群定数を温度について線形に内挿する。
参照解の OpenMC も温度の内挿のモードで解く。OpenMC はこのモードで、両側の温度の S(α,β) を温度（kT）の比で確率的に
選ぶ（src/thermal.cpp）ので、期待値は同じ線形の内挿になり、扱いがそろう。S(α,β) を持たない核種は、xs も OpenMC の
ライブラリも計算する温度（565 K、等温温度係数の 560 K と 570 K）そのもので作り、内挿するのは S(α,β) だけにする。

KENO の値は、仕様書が S(α,β) の温度を補正した値（550 K と 600 K の内挿、または k の補正）を使う。どちらも 565 K の
S(α,β) を内挿で作ったものに当たるため。

## 確認

`thermal_scattering_is_interpolated_in_temperature`、`other_nuclides_are_not_interpolated`。

## 結果

[ADR-0047](0047-choose-temperatures-per-nuclide.md) の「温度方向に補間しない」の例外は、S(α,β) の熱中性子の散乱行列だけである。
