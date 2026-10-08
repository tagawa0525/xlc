---
status: open
labels: [xs]
opened: 2026-10-08
closed:
milestone:
blocked-by: []
resolution:
---

# 点の数が等しい対応点の補間（INT=12〜15）を実装していない

## 問題

[ADR-0057](../adr/0057-interpolate-outgoing-distributions-as-endf.md) は、出射分布の入射エネルギーの補間則を ENDF のとおりに扱うと決めている。対応点の補間（INT=11〜15）は、
隣り合う分布の点の数が違うと定義できないので単位基底に切り替えるが、点の数が等しいときの本来の対応点の補間は
実装しておらず、エラーにしている（crates/xs/src/endf/mf5.rs の `incident_laws`。テスト
`corresponding_points_with_equal_grids_are_not_supported`）。いまの対象核種で、これに当たるファイルは見つかっていない。

関係する ADR: [ADR-0057](../adr/0057-interpolate-outgoing-distributions-as-endf.md)

## 完了の条件

点の数が等しい対応点の補間を実装し、テストで確かめる。対象核種で使われないままなら、ADR でエラーにすることを決める。

## 経緯

- 2026-10-08: ADR のレビューで、決定と実装の食い違いとして見つかった
