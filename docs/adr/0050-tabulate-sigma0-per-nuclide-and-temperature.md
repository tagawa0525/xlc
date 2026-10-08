---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
---

# ADR-0050: 多群ライブラリは核種・温度ごとに σ0 ごとの群定数を持つ

## 背景

M1 で、温度と自己遮蔽に依存する多群ライブラリを作る。ライブラリの中身と、背景断面積 σ0 の格子を決める必要があった。

## 検討した案

- 散乱行列も σ0 ごとに持つ: 効果が小さいので後回しにする
- 散乱行列は無限希釈のものだけを持つ（採用。自己遮蔽した核種での行の扱いは [ADR-0092](0092-match-scattering-rows-to-shielded-elastic.md)）

## 決定と理由

- ライブラリは、核種・温度ごとに、σ0 ごとの全断面積・弾性散乱・核分裂・νσf と、χ、無限希釈の散乱行列を持つ
- σ0 は 1, 10, 30, 100, 300, 1e3, 1e4, 1e5, 1e10 b とし、全核種で共通にする
- 散乱行列は σ0 ごとに持たない（散乱行列の自己遮蔽は効果が小さいので後回しにする）

M1 の計画は、散乱行列の行を自己遮蔽に合わせる扱いもこの項目に書いていた。その部分は別の決定として [ADR-0092](0092-match-scattering-rows-to-shielded-elastic.md) に分けた
（のちに [ADR-0076](0076-store-legendre-scattering-moments.md) で置き換え）。

## 確認

`processes_each_temperature_and_sigma0`、`constant_total_gives_infinite_dilution_for_every_sigma0`、
`shielded_cross_sections_at_infinite_dilution_match_averages`。

## 結果

- σ0 ごとの重みは [ADR-0051](0051-weight-shielded-data-by-nr.md)（のちに [ADR-0072](0072-weight-shielded-data-by-flux-calculator.md)）
- M1 の計画では、M2 までは格子計算は指定した温度の無限希釈の値を使うとしていた。格子計算で自己遮蔽を使い始めたのは
  M2（[ADR-0068](0068-solve-pinwise-slowing-down.md)、[ADR-0074](0074-shield-outside-range-by-equivalence.md)）
- ファイルの形式は M2 で MessagePack に改め、形式の版も廃止した（[ADR-0011](0011-write-data-files-in-messagepack.md)、[ADR-0012](0012-drop-backward-compatibility.md)）
- ポイントワイズの断面積などの追加は [ADR-0073](0073-store-pointwise-data-in-library.md)、散乱行列の Legendre モーメントは [ADR-0076](0076-store-legendre-scattering-moments.md)
