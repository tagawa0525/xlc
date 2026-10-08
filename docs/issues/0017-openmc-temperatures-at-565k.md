---
status: open
labels: [verification, xs]
opened: 2026-10-08
closed:
milestone: M3
blocked-by: []
resolution:
---

# 565 K の OpenMC の参照解の温度の扱いを xlc とそろえる

## 問題

[ADR-0080](../adr/0080-interpolate-sab-in-temperature.md) は、565 K の H in H2O の群定数を 550 K と 573.6 K の表から温度で内挿し、参照解の OpenMC も温度の内挿の
モードで解くとしている。いまの examples/vera/openmc_library.py は、燃料以外の核種と H in H2O を 600 K だけで作り、
examples/vera/openmc_model.py は最も近い温度のモード（`"nearest"`）で解く。OpenMC の内挿のモードは、温度を1つしか
持たない核種があると全体を最も近い温度のモードに戻す（OpenMC の src/nuclide.cpp）。また、xlc は H-1 の表全体
（自由原子の全断面積とポイントワイズの断面積を含む）を内挿するが、OpenMC は H-1 を 565 K で処理していればその値を使う。

関係する ADR: [ADR-0080](../adr/0080-interpolate-sab-in-temperature.md)

## 完了の条件

565 K（と 560 K、570 K）の OpenMC のライブラリと入力を作り、内挿のモードで解けることを確かめる。H-1 の自由原子の
データの扱いの違いが k と反応率に効く大きさを見積もり、必要なら xlc の内挿の範囲を見直すか、差を記録する。

## 経緯

- 2026-10-08: ADR のレビューで、[ADR-0080](../adr/0080-interpolate-sab-in-temperature.md) が OpenMC 側の準備を前提にしていることが分かった
