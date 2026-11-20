---
type: Issue
status: open
labels: [xs, vera]
opened: 2026-11-21
closed:
milestone:
blocked-by: []
resolution:
---

# ADR-0274 の宣言で変わる Cd-106、Cd-111 の非分離共鳴の群定数が、VERA の AIC の制御棒の結果に効く大きさを確かめていない

## 問題

- ENDF/B-VII.0〜VIII.1 の Cd-106 と Cd-111 は、非分離共鳴（LSSF=0）の MF3 の捕獲のバックグラウンドが負で（MT102 の最小 −0.0647 b、−0.146 b）、
  xs は入力の宣言 `unresolved-negative-background = "scale"` なしで止まる（ADR-0274）。AIC の Cd を輸送する入力（examples/vera/xs.toml、
  examples/vera/endfb70/xs.toml、examples/depletion/p2/xs.toml、examples/core/xs.toml、examples/m3-nuclides/xs.toml）に宣言を足したが、ライブラリは
  作り直していない（ステージは changed）
- 宣言で変わる遮蔽した捕獲（VIII.0、565 K、範囲にかかる群の最大）は、Cd-111 が σ0 = 1e3 b で +0.013 %、1 b で +1.67 %、Cd-106 が 1e3 b で +0.086 %、1 b で
  +13.2 %（examples/depletion/verification.md の「版ごとのライブラリの非分離共鳴の負のバックグラウンド」）。変わるのは非分離共鳴の範囲（Cd-106 は
  6〜100 keV、Cd-111 は 2.3〜100 keV）にかかる群だけ
- AIC の Cd の吸収のうち非分離共鳴の範囲が占める割合と、AIC の中の Cd の σ0 は、既にある出力から出していない。KENO との比較（VERA の
  Problem 2〜5 の AIC の状態）と制御棒の価値に効く大きさかは確かめていない

関係する ADR: [ADR-0274](../adr/0274-scale-sampled-reactions-with-negative-unresolved-backgrounds.md)

## 完了の条件

- [ ] VERA のライブラリを作り直したとき（ADR-0258 の作り直しの残高）に、AIC の状態の k の、宣言の前後の差を記録に書く

## 経緯

- 2026-11-21: #0319 のブランチで、AIC の Cd を輸送する入力に宣言を足し、作り直しをしないまま残すことにして起こした
