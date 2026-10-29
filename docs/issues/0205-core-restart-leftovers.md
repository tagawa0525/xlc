---
status: open
labels: [core, depletion]
opened: 2026-10-30
closed:
milestone: M4
blocked-by: []
resolution:
---

# core の restart と装荷替えの残り（戻した集合体、燃焼しないケースの restart、部分ノードの回転）

## 問題

ADR-0220 の restart は、前のサイクルの燃焼の履歴 1 本の最後の点から、集合体の ID で状態を引き継ぐ。次は入れていない（ADR-0220 の帰結）:

- 2 サイクル以上前に取り出してプールに置き、戻した集合体: 前のサイクルの履歴になく、`fresh` でもないのでエラーで止まる。その集合体の最後の
  サイクルの履歴を読む形（restart を複数の履歴にする、集合体ごとの状態のファイル）が要る。プールの間の崩壊（取り出してから戻すまでの日数）も要る
- 燃焼しないケースの restart: Layout が restart のサイクルと燃焼しないケースの組み合わせをエラーにする。Problem 10 の BOC の HZP の解は、restart の
  燃焼のケースの燃焼点 0 で得られるが、燃焼のない状態のケース（ピン出力、複数の状態の並列の解。ADR-0146）として解く形が要るかは、Problem 10 の入力を
  作るときに決める
- `radial-nodes` が 1 でない燃焼: 部分ノードの状態を集合体の自分の座標系で持ち、回転で並べ替える（ADR-0036）。燃焼そのものが `radial-nodes = 1`
  だけ（ADR-0190）
- 挿入物（`inserts`。Pyrex を抜く集合体）と設計の違う集合体（IFBA、WABA、ブランケット）のテーブルの履歴の扱い（ADR-0215、#0190）

レビュー（ブランチ feat/core-restart-reload の 1 周目）で挙がり、結果を変えないので直さなかった整理:

- examples/verification_tables.py の `difference` は `signed_exponent` と桁の数だけ違う複製で、負にならない量（燃焼度と数密度の差の最大）にも
  `+` を付ける。桁の引数と符号なしの形を `signed_exponent` に足せばまとまる
- `rho_pcm` が examples の比較のツールに 4 つ目の複製（compare_core_restart.py）。examples/depletion/core_history.py にまとめられる
- crates/core/src/micro.rs の `MicroChain::auxiliary_columns` は呼ぶたびに作り直す（ノードごとに 2 回）。`MicroChain::new` で 1 度作って持てる
- crates/core/src/depletion.rs の restart の組成の置き方は、`placed` を 2 度引き、`inherited` から移せる組成を複製する

関係する ADR: ADR-0220

## 完了の条件

Problem 10 の入力（仕様書の付録 K の装荷図）を書いて、上のそれぞれが要るかを決め、要るものを入れる。

## 経緯

- 2026-10-30: ADR-0220 の実装（ブランチ feat/core-restart-reload）の残りとして起こした
