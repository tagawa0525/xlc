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
- 2 周目: examples/depletion/history/compare_rods_switch.py の `compare()` と `report()` を呼ぶテストが無い（`continued` の読み方の直しも
  テストなし）。比較のツールのテストと実装は 1 つのコミット（70f91dde3）。`compare()` は core が参照の `continued` の燃焼度で切り替えたこと、
  20 MWd/kg まで H_r が 1、減速材の密度の重みが 0 であることを照らさず、JSON の `density_weight` は w_m でなく H_m。`nodes`、`at_burnup`、
  `production`、`scaled` は compare_core_history.py などの複製。燃焼度の軸を最初の軸と置く。切り替えの後の点が無いと `report()` が空の列で
  止まる。`test_the_overcount_adds_the_extra_weights_on_the_linearizations` は値を見ず符号だけを見る
- 2 周目: 2 つの変数のテーブル（residual-rods.toml）は減速材の密度の残りを制御棒の値 out だけで持つので、制御棒を入れたノードで w_m が 0 でないと
  core が止まる（ADR-0228 の分け方を制御棒を入れたノードで使えない。減速材の密度の履歴の分岐を制御棒を入れた値でも解けば足せる）
- 3 周目: 記録の「制御棒を途中まで入れて抜く燃焼（2A）」の表は Pu-239〜241 だけを描き、本文の 5 核種の最大（Xe-135 の 0.459%）は JSON にだけある。
  core の `add_histories` と `node_micro` は同じノードの a_j を 2 度作る

関係する ADR: ADR-0228

## 完了の条件

上の指摘を、直すか、直さない理由を書いて閉じる。

## 経緯

- 2026-11-04: feat/0190-rod-history のレビュー 1〜3 周目で残した
- 2026-11-05: ADR-0215 の「出力の形」の 2 文を直した後（4bf10676e、c0938d5f5）にブランチ全体へ回した周（high）で残した文言の指摘（どれも
  non-blocking。ユーザーは 2026-11-05 00:44 に文言の整合で周を重ねないと決めた）: 条件 j を公称のテーブルの格子点に引く文の主語が両方の
  変数になっているが `rods` は制御棒の値の公称の格子点を固定で使う（`rod_history_condition`）、制御棒の軸を持たないテーブルは値 `out` にしか
  渡せない規則（residual.rs の 247〜249 行）が書かれていない、ハイブリッドの「同じ軸と点」は `rod_value` で値を選んだ後の照合である
  （`check_hybrid`）、「`rods` は制御棒の軸が要る」と「H^0 = 0、H^b = 1」が同じ段落に 2 度ある、「## 変更」の 2026-11-05 の項の
  「変えた部分」が検査の規則の変更に読める（実際は `moderator-density` だけを書いていた括弧を変数ごとにした文言の直し）
