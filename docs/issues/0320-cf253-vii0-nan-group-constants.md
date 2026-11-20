---
type: Issue
status: open
labels: [xs, depletion]
opened: 2026-11-20
closed:
milestone:
blocked-by: []
resolution:
---

# ENDF/B-VII.0 の Cf-253 の σ0 が 10 b 以下の群定数に NaN が出て、xs が止まる

## 問題

- ENDF/B-VII.0 の Cf-253 の `xs nuclide`（565 K の 1 点）は、ライブラリを書く前の検査（crates/format の `NuclideData` の温度ごとの値の検査）で
  「Cf-253 at 565 K: non-finite value NaN」で止まる（止まった文と一覧の出どころは examples/depletion/verification.md の「中性子の評価のすべての
  `xs nuclide`」）。VII.1 と VIII.0、VIII.1 の Cf-253 は通る
- 入力を変えた切り分け（xlc 49e3045c7 の bin/xs、手元、2026-11-20）: σ0 を [1e10] だけ、[1e10, 1e5]、[1e10, 1e3] にすると通り、[1e10, 10] と
  [1e10, 1] では止まる。温度（1 K、293.6 K、1200 K）とポイントワイズの範囲（[4, 1.353e6] と [4, 9.877] eV。後者は 293.6 K だけ）によらず止まる。NaN の入る量
  （χ、散乱行列、σ0 ごとの遮蔽した断面積、反応）と群は確かめていない
- VII.0 の Cf-253 の評価は MF3 が MT1、2、18、102 だけで、核分裂（MT18）は 1e4〜1.1e4 eV の外で 0（MF3 の TAB1、2026-11-20）。分離共鳴は
  SLBW で、J が標的のスピンと合わない（xs の警告「SLBW resonances with L=0 have J [3.5], not allowed for the target spin 3.5」）
- 候補（確かめていない）: σ0 が小さいときの中性子束の重み（減速計算 `multigroup/slowing.rs` の分母 σt + σ0 − …、範囲の上の 1/(σt + σ0)）が
  0 か負になる点、核分裂の 0 の群での 0/0
- 再現: examples/depletion/process-survey/VII.0 で `../../../../bin/xs --threads 2 nuclide xs.toml xs-range.toml Cf-253 <出力>`

関係する ADR: [ADR-0274](../adr/0274-scale-sampled-reactions-with-negative-unresolved-backgrounds.md)

## 完了の条件

- [ ] NaN を作る計算と原因を特定し、VII.0 の Cf-253 が `xs nuclide` を通るか、評価の側の誤りとして扱いを決める

## 経緯

- 2026-11-20: 4 つの版の中性子の評価のすべてを `xs nuclide` に通し（#0292）、VII.0 の Cf-253 が止まって起こした
