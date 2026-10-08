---
status: accepted
date: 2026-10-07
supersedes: []
superseded-by:
issues: []
---

# ADR-0057: 出射分布の入射エネルギーの補間則は ENDF のとおりに扱う

## 背景

Cr の MF6 は、入射エネルギーの間の補間に対応点の補間（INT=12）を使う。補間則の扱いを決める必要があった。

## 検討した案

- ふつうの線形（INT=2）を、NJOY と同じく単位基底に置き換える: NJOY は散乱源を滑らかにするためにそうするが、
  ファイルの指定と違う
- ファイルの指定どおりにする（採用）

## 決定と理由

出射分布の入射エネルギーの補間則は ENDF のとおりに扱う。単位基底（INT=21〜25）は出射エネルギーの範囲を正規化して
補間する。ふつうの線形（INT=2）は、ファイルの指定どおりにする。対応点（INT=11〜15）は、隣り合う分布の点の数が違うと
定義できないので、NJOY と同じく単位基底に切り替える（Cr はどの反応もこれに当たる）。

## 確認

`incident_interpolation_follows_its_law`、`incident_interpolation_schemes`、`unit_base_keeps_the_shape_between_incident_energies`、
`corresponding_points_with_equal_grids_are_not_supported`。

## 結果

INT=2 の分布では、NJOY と補間の仕方が違う。
