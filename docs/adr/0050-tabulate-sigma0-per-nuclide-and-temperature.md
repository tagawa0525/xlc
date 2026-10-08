---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
---

# ADR-0050: 多群ライブラリは核種・温度ごとに σ0 ごとの群定数を持つ

## 背景

M1 で、温度と自己遮蔽に依存する多群ライブラリを作る。ライブラリの形式と中身、背景断面積 σ0 の格子、散乱行列の
自己遮蔽の扱いを決める必要があった。

## 検討した案

- 散乱行列も σ0 ごとに持つ: 効果が小さいので後回しにする
- 散乱行列は無限希釈の形のまま、行の和を自己遮蔽した弾性散乱に合わせる（採用）

## 決定と理由

- ライブラリは、核種・温度ごとに、σ0 ごとの全断面積・弾性散乱・核分裂・νσf と、χ、無限希釈の散乱行列を持つ。
  温度と σ0 に依存する表を持てるように、形式を変える
- σ0 は 1, 10, 30, 100, 300, 1e3, 1e4, 1e5, 1e10 b とし、全核種で共通にする
- 散乱行列は無限希釈の形のまま、行の和を自己遮蔽した弾性散乱断面積に合わせる（散乱行列の自己遮蔽は効果が小さいので
  後回しにする）
- M2 までは、格子計算は指定した温度の無限希釈の値を使う

## 確認

`processes_each_temperature_and_sigma0`、`constant_total_gives_infinite_dilution_for_every_sigma0`、
`shielded_cross_sections_at_infinite_dilution_match_averages`。

## 結果

- σ0 ごとの重みは [ADR-0051](0051-weight-shielded-data-by-nr.md)（のちに [ADR-0072](0072-weight-shielded-data-by-flux-calculator.md)）
- ファイルの形式は M2 で MessagePack に改め、形式の版も廃止した（[ADR-0011](0011-write-data-files-in-messagepack.md)、[ADR-0012](0012-drop-backward-compatibility.md)）
- ポイントワイズの断面積などの追加は [ADR-0073](0073-store-pointwise-data-in-library.md)、散乱行列の Legendre モーメントは [ADR-0076](0076-store-legendre-scattering-moments.md)
