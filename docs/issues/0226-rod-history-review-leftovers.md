---
status: open
labels: [core, lattice, format]
opened: 2026-11-04
closed:
milestone: M4
blocked-by: []
resolution:
---

# 制御棒の履歴（H_r）と変数の間の分け方のレビューで残した指摘

## 問題

ブランチ feat/0190-rod-history（#0190 の制御棒の履歴 H_r と ADR-0228）のローカルのレビューで、再現を書けない（契約の外の入力、
効率、重複、文言）として直さなかった指摘。

- core の `linearization_weights`（crates/core/src/case.rs）は Σ_k w_k P_k がちょうど 0 のときだけ止め、打ち消し合って 0 に近いときは
  a_j が大きく（符号も逆に）なりうる（ADR-0228 の帰結。減速材の密度の履歴が公称の逆の側（w_m < 0）で制御棒を入れたノード。相対の
  閾値は決めていない。P / P_j の上限と同じく #0190 の「残ったレビューの指摘」）
- core の制御棒の挿入の割合は、棒のない区間でない区間をすべて数え、制御棒の値（aic、b4c など）を区別しない。履歴の分岐は 1 本
  （Hᵇ = 1）で、どの吸収体で燃焼させたかをテーブルが持たない
- lattice の `residuals`（crates/lattice/src/residual.rs）は `rods` の Hᵇ = 1 を format の `rod_history_condition` と別に持つ。format の
  `rod_history_condition` と `moderator_history_condition` は条件の格子点の組み立てを重ねて持つ
- lattice の条件の食い違いのエラーの文は、`rods` でも「the moderator at the density 1」と書く
- lattice の `rod_value` は制御棒の軸のないテーブルも複製し、`run` は同じファイルを制御棒の値ごとに読み直す。履歴の残りを持つ
  テーブルを分岐に渡すと、選んだ値の軸に無い `out` の残りで検査が誤った文で止まる
- core の `node_micro` の Σ_j a_j ΔN̄ᵢ^j は最初の変数の長さで zip する（長さはテーブルの検査が揃える）。`expected` の名前が
  `linearization_weights` の Σ w_k P_k と重なる

関係する ADR: ADR-0228

## 完了の条件

上の指摘を、直すか、直さない理由を書いて閉じる。

## 経緯

- 2026-11-04: feat/0190-rod-history のレビュー 1 周目で残した
