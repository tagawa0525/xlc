---
status: closed
labels: [verification]
opened: 2026-10-14
closed: 2026-11-19
milestone:
blocked-by: []
resolution: "wontfix: #0224 に含まれる。examples/m1 のステージは今の bin で再実行してあり、残る記録の本文は #0224 の一覧にある"
---

# M1 の記録を高速化の後のコードで再測定する

## 問題

高速化の後の記録の測り直し（ブランチ docs/remeasure-after-speedups、Slurm のジョブ 332、分岐点 c59097f0、2026-10-14）は、
DVC のステージのある例と、examples/vera の 1C〜2Q、examples/mox を再測定した。examples/m1 の記録の数値（examples/m1/verification.md、
8879202c、2026-10-12）は再測定していない。M1 も ADR-0156 で変わった自己遮蔽のコードを使う。

関係する ADR: なし

## 完了の条件

- examples/m1 を今のコードで走らせ（`run.sh` の手順）、記録の数値を書き直す。前の値との差を書く

## 経緯

- 2026-10-14: 再測定のレビューで、M1 が抜けていることが分かった
- 2026-11-19: issue の棚卸しで閉じた。examples/m1 のステージは新しい既定の後の再実行で再実行した。記録の k と表の書き直しは #0224 の一覧が受け持つ（前の値との差は、記録は今の結果だけという規則で書かない）
