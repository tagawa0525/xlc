---
status: closed
labels: [lattice, core, depletion]
opened: 2026-10-23
closed: 2026-11-07
milestone: M4
blocked-by: []
resolution: ADR-0190
---

# core のノードの出力を核分裂のエネルギーで重みづける

## 問題

なし（ADR-0190 の「## 変更」の 2026-11-07 の項で解決した。熱水力のノードの出力の重みは #0266）

関係する ADR: ADR-0190

## 完了の条件

lattice がテーブルの状態ごとに κΣf（Q で重みをつけた核分裂の断面積）を書き、core がノードの出力の割合をそれで求める。差の大きさを
小さな炉心（燃焼度の違う集合体を並べたもの）で測り、検証の記録に書く。Problem 9 の燃焼（docs/plans/m4-depletion.md の「進め方」9）の前に行う。

## 経緯

- 2026-10-23: ADR-0190 の最初の縦切りで、テーブルに κΣf がないため核分裂率の割合で進め、この issue に残した
- 2026-11-07: テーブルの状態ごとの κΣf と、core のノードの出力の割合を κΣf φ V で求める形を入れ、小さな炉心で差を測って閉じた（ADR-0190 の「## 変更」、examples/depletion/verification.md の「ノードの出力の重み（小さな炉心）」、c899bf3d9 の版、2026-11-07）
