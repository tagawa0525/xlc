---
status: closed
labels: [lattice, depletion, performance]
opened: 2026-10-15
closed: 2026-10-24
milestone: M4
blocked-by: []
resolution: ADR-0197
---

# 燃焼の中点で自己遮蔽を解き直さない場合の誤差と時間

## 問題

ADR-0173 は、CE/CM の始まりと中点の両方の輸送で自己遮蔽を解き直す（案 F1）。格子計算の燃焼の時間の大部分は自己遮蔽で
（VERA 1C のピンセルで 95%、2A の 1 ステップで 73%。examples/depletion/verification.md の「格子計算の燃焼」）、中点で始まりの
微視的断面積を使えば（案 F2）、自己遮蔽の回数はほぼ半分になる。誤差は、中点までの組成の変化（初めの Xe、Sm の平衡、Gd）で
決まり、測っていない。

関係する ADR: ADR-0173、ADR-0197

## 完了の条件

集合体（毒物なし、Gd）の燃焼で、F1 と F2 の k∞ と数密度の差と時間を測り、ADR で頻度を決める。

## 経緯

- 2026-10-15: ADR-0173 の「帰結」から
- 2026-10-24: `[depletion]` の `shielding = "step-start"` で F2 を選べるようにし、2A と 2O の 0〜60 MWd/kg で F1 と比べた
  （examples/depletion/verification.md の「自己遮蔽の頻度（#0134）」）。k∞ の差は最初のステップで +7〜+9 pcm、長いステップで最大
  −17.9 pcm、数密度の差は 0.5% 以下で、自己遮蔽の時間はほぼ半分。既定は F1 のまま（ADR-0197）
