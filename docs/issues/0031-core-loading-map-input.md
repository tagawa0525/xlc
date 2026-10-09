---
status: closed
labels: [core]
opened: 2026-10-09
closed: 2026-10-09
milestone: M3
blocked-by: []
resolution: ADR-0105
---

# 炉心の装荷図・配置の入力形式（段階 6 で決める。形式はユーザーが決める）

## 問題

core の入力は1つの少数群定数のファイルで均質な直方体だけを書け、種類の違うノード（燃料と反射体、種類の違う
燃料）を並べられない。反射体の定数の等価性（ADR-0095）は、core のライブラリの API で問題を組んで確かめる。

ブランチ feat/reflector-constants で暫定の形を試した（84d96d1 のテストと f901f07 の実装。ブランチからは外した）:
`few-group` を名前つきのファイルの表にし、`layout` に x-y のノードの名前を上の行（+y）から並べる（z は一様、行と列の
数は y と x のノード数に等しい）。core は位置ごとの群の中性子束と核分裂率を出す。

```toml
few-group = { fuel = "fg-fuel.msgpack", reflector = "fg-reflector.msgpack" }
layout = ["fuel reflector"]
```

関係する ADR: ADR-0105

## 完了の条件

M3 の進め方の 6 で、装荷図（集合体の種類と位置、軸方向の層、反射体の向き、状態の分岐のテーブル）の入力形式を、
ほかのコード（PARCS、SIMULATE-3、VERA の入力の `assm_map` など）と比べてユーザーが決め、ADR にして実装する。

## 経緯

- 2026-10-09: 反射体の定数の作業で暫定の形を作り、話題を分けるために外した（ユーザーの方針）
- 2026-10-09: ユーザーと先行事例（VERAin、PARCS）を比べ、集合体の ID で追う形（シャッフリング）に決めた（ADR-0105）
