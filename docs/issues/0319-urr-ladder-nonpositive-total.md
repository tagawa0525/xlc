---
type: Issue
status: closed
labels: [xs, depletion]
opened: 2026-11-20
closed: 2026-11-20
milestone:
blocked-by: []
resolution: ADR-0274
---

# 非分離共鳴のラダーのサンプルの全断面積が正でなく、xs が止まる（ENDF/B-VII.0 の Ta-181、Es-253）

## 問題

- ENDF/B-VII.0 の Ta-181 と Es-253 の `xs nuclide`（565 K の 1 点）は、非分離共鳴の自己遮蔽のラダーのサンプリング（`resonance/ladder.rs` の
  `bondarenko` の和、`add_samples`）で、サンプルの点の全断面積が 0 以下になり、`UrrError::NonPositiveTotal`（「the sampled total cross section
  is not positive … eV from … eV」）で止まる。Ta-181 は 330 eV（非分離共鳴の下端）の窓の −251.8 eV の点、Es-253 は 8500 eV の窓の +5.06 eV の点
  （止まった文と一覧の出どころは examples/depletion/verification.md の「中性子の評価のすべての `xs nuclide`」）。`xs check` は群定数を作らない
  ので、読み込みの走査（ADR-0266）では止まらない
- 4 つの版でこの形で止まるのはこの 2 つだけ。VII.1、VIII.0、VIII.1 の同じ核種は通る
- 候補（確かめていない）: 2 つとも、非分離共鳴の範囲の MF3 のバックグラウンドが負。Ta-181 は非分離共鳴の範囲（330〜5000 eV、LRF=2）の
  MF3 の 7 点の MT1、MT2 がどれも負で、MT102 は 330 eV で −1.61 b、5000 eV で −0.113 b。Es-253 は範囲（101.2〜10000 eV）の上端の 10 keV で
  MT1 が −6.82 b、MT102 が −0.698 b（ファイルの MF3 の TAB1 を読んだ値、2026-11-20）。ラダーは弾性散乱が負の点だけを 1 μb にし（NJOY PURR と
  同じ。ADR-0177）、捕獲のバックグラウンドが負の点は全断面積が 0 以下になりうる
- NJOY PURR が捕獲のバックグラウンドの負の点をどう扱うかは確かめていない
- 再現: examples/depletion/process-survey/VII.0 で `../../../../bin/xs --threads 2 nuclide xs.toml xs-range.toml Ta-181 <出力>`（Es-253 も同じ）

関係する ADR: [ADR-0274](../adr/0274-scale-sampled-reactions-with-negative-unresolved-backgrounds.md)

## 完了の条件

- [x] 負の MF3 のバックグラウンドで全断面積が 0 以下になる点の扱いを決め（NJOY PURR の扱いを確かめてから）、VII.0 の Ta-181 と Es-253 が
  `xs nuclide` を通る

## 経緯

- 2026-11-20: 4 つの版の中性子の評価のすべてを `xs nuclide` に通し（#0292）、VII.0 の 2 つが止まって起こした
- 2026-11-20: 原因は MF3 の捕獲のバックグラウンドが負で、各サンプルに定数で足すと共鳴の間の点の捕獲が負になること（Ta-181 は 330 eV で
  −1.61 b、Es-253 は 10 keV で −0.698 b）と分かった。NJOY PURR は止まらず、確率テーブルに負の値を残し、Ta-181 の σ0 = 1 b の直接のサンプルの
  全断面積は −195 b（references/njoy-runs/issue-0319/）。同じ状態の評価を 4 つの版で数え（25 評価）、既定はバックグラウンドの値で止め、
  宣言 `unresolved-negative-background = "scale"` でその反応のサンプルを (平均 + バックグラウンド)/平均 倍にすると決めた（ADR-0274）。宣言した
  Ta-181 と Es-253 は `xs nuclide` を通り、閉じた
