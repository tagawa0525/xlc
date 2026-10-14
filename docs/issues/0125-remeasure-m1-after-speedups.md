---
status: open
labels: [verification]
opened: 2026-10-14
closed:
milestone:
blocked-by: []
resolution:
---

# M1 の記録を高速化の後のコードで測り直す

## 問題

高速化の後の記録の測り直し（ブランチ docs/remeasure-after-speedups、Slurm のジョブ 332、分岐点 d07475d3、2026-10-14）は、
DVC の段階のある例と、examples/vera の 1C〜2Q、examples/mox を測り直した。examples/m1 の記録の数値（examples/m1/verification.md、
1f4fcaaf、2026-10-12）は測り直していない。M1 も ADR-0156 で変わった自己遮蔽のコードを使う。

関係する ADR: なし

## 完了の条件

- examples/m1 を今のコードで走らせ（`run.sh` の手順）、記録の数値を書き直す。前の値との差を書く

## 経緯

- 2026-10-14: 測り直しのレビューで、M1 が抜けていることが分かった
