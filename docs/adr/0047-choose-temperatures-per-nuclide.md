---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
---

# ADR-0047: ライブラリの温度は核種ごとに指定し、格子計算は材料の温度がライブラリにあることを求める

## 背景

H in H2O の S(α, β) は 283.6〜800 K しかなく、水を 900 K 以上で扱う必要もない。全核種に同じ温度の組を求めると、
水を作れない。

## 検討した案

- 全核種で共通の温度の組: H in H2O を 900 K 以上で作れない
- 核種ごとに温度を指定する（採用）

## 決定と理由

温度は核種ごとに指定する（燃料と被覆管は 293.6, 600, 900, 1200 K、減速材は 293.6, 600 K のように選ぶ）。格子計算は
群定数を温度方向に補間せず、材料の温度がライブラリにあることを求める。例外は H in H2O のような S(α,β) の熱中性子の
散乱行列だけで、それは [ADR-0080](0080-interpolate-sab-in-temperature.md) で決める。

## 確認

`temperature_must_be_in_the_library`、`temperatures_are_looked_up_exactly`。

## 結果

S(α,β) を持たない核種は、計算する温度そのものでライブラリを作る（565 K、560 K、570 K など。[ADR-0080](0080-interpolate-sab-in-temperature.md)）。
