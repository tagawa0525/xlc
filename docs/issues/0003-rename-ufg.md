---
status: open
labels: [lattice, refactor]
opened: 2026-10-08
closed:
milestone:
blocked-by: []
resolution:
---

# 「UFG」という名前が中身（ポイントワイズ）と合わない

## 問題

lattice のピンごとの減速計算は、断面積をポイントワイズの点の格子で持ち、点の間を線形として解く方式（CENTRM と同じ
系統）だが、コード（`ufg` モジュール、`UfgError` など）と docs では「超微細群（UFG）」と呼んでいる。UFG は普通、
等レサジーの数万群に群平均の断面積を置く方式（AEGIS、nTRACER の PSM）を指すので、誤解を招く。

## 完了の条件

名前を決め（例: `pointwise`、`slowing_down`）、コードと docs の呼び方をそろえる。

## 経緯

- 2026-10-08: ユーザーとの議論で、UFG とポイントワイズの違いが話題になった
