---
status: open
labels: [core, lattice, format]
opened: 2026-10-30
closed:
milestone: M4
blocked-by: []
resolution:
---

# 履歴の反応率の補正（ADR-0215 の案 S2）のレビューで残したもの

## 問題

ブランチ feat/0190-history-micro-sigma（#0190）のローカルのレビュー（1 周目、HEAD 23863f7f）で、再現できる誤りでない指摘を残した。
どれも出力の数値は変えない。

関係する ADR: ADR-0215

- 型: crates/core/src/case.rs の `node_micro` は、線形の補正の係数を読む `correct` と同じ型 `MicroAt` を返し、`fission`（Σf の係数でもある）と
  収率の境界より下の割合だけを補正する。補正した値が線形の補正に渡ると、ADR-0215 の案 S2 の「係数は補正しない」が崩れ、履歴の分岐の点に
  (cᵢ^j − cᵢ(j)) ΔN̄ᵢ が残る。今は `node_constants` が `interpolate_micro` を、反応率が `node_micro` を別に呼ぶので起きないが、型が止めない
  （反応率の断面積を別の型にする候補）
- テスト: crates/core/src/depletion/tests.rs の `micro_depletion_rates_take_the_weighted_rate_residuals_of_the_history` の w = 0 の半分は、
  Δσ を 0 にした `histories` を持つテーブルと比べ、`histories` を消したテーブルとは比べない（どちらも w = 0 で同じ経路を通る）。注釈は
  「補正のないテーブル」と書く。case の `nominal_histories_leave_the_rate_cross_sections_bit_identical` は `histories` を消して比べている
- エラーの文: 核分裂する追う核種の σf の群の和が内挿で 0 以下になると、`interpolate_micro` は割合を `None` にし、w が 0 でないノードでは
  `add_rate_residual` が「核分裂の有無が違う」と報告する（本当の原因は σf の和）
- 重複: crates/lattice/src/residual.rs の `residuals` で、Δσ と R が同じ燃焼度の点の走査と条件の点の引き直しを 2 回する（1 回の走査で組を返せる）。
  crates/core/src/state.rs の `scaled_rate`、`add_rate`、`add_rate_residual` の反応の経路の足し方が `scaled_micro`、`accumulate_micro` と同じ形で、
  (MT、経路の数) の形の比べ方が lattice の `rate_difference`、形式の `validate_rate_residuals`、`validate_micro` にある
- 効率: crates/core/src/depletion.rs の `Solve` が解ごとにノードの履歴の変数のコピーを持つ（平衡 Xe の反復でも解ごと）。反応率を作るときに
  読むだけなので、`n` と同じく引数で渡せる。炉心で効くかは測っていない

## 完了の条件

- 上の項目ごとに、直すか、直さない理由を書く

## 経緯

- 2026-10-30: feat/0190-history-micro-sigma のレビューで起こした
