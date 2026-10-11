---
status: accepted
date: 2026-10-08
supersedes: []
superseded-by:
issues: ["#0004"]
requires: []
---

# ADR-0088: 枠のあるピンの等価なセルは、減速材のリングの外に枠の物質のリングを置く

## 背景

支持格子の枠のある格子では、格子の Dancoff 係数に枠の影が入る。等価なセル（[ADR-0068](0068-solve-pinwise-slowing-down.md)）に枠がないと、その影を水で作る
ことになり、枠の Zr が減速材の中性子束を変える効果も落ちる。VERA 1B では、枠による k の変化の OpenMC との差
（−170 pcm。2B は −126 pcm）の大部分がこれによる（examples/vera/verification.md、612913c と 9fb742c、2026-10-08）。

## 検討した案

- 枠を減速材に混ぜる: k はほぼ同じになる（15 pcm 以内）が、Zr を薄めて枠の中の自己遮蔽を落とし、減速材のリングの
  断面積から枠の核種を除く手間が要る
- 減速材のリングの外に、枠と同じ面積の枠の物質のリングを置く（採用）

## 決定と理由

枠のあるピンの UFG の等価なセルは、減速材のリングの外に、枠と同じ面積の枠の物質のリングを置く。外径は、枠のリングを
含めたセルで Dancoff 係数を合わせて決める。枠のリングの断面積は格子には使わず、枠は今までどおり塊として遮蔽する
（枠は隣のセルと共有されるため。[ADR-0086](0086-shield-structural-lumps.md)）。枠のリングの SPH 因子は 1 のままにする。

## 確認

`frames_are_the_outer_ring_of_the_cell`、`fixed_rings_keep_one`。枠による k の変化の OpenMC との差は、1B で −170 → +40 pcm、
2B で −126 → +38 pcm になった（examples/vera/verification.md、1d4606f、2026-10-08）。OpenMC の2つの計算の差の統計誤差は
1σ で約 10 pcm。

## 結果

残りの +40 pcm（約 4σ）の原因は分けていない（#0004）。
