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

2 周目のレビューで挙がり、直さなかったもの:

- 引き継ぐ集合体の重金属の質量は、最初の状態の並びのノードのファイル（部分挿入のノードは最初の区間のファイル）のテーブルから求める。制御棒の区間を
  別のファイルのテーブルにした炉心で、前のサイクルと違う区間のファイルになる位置へ動かすと、同じ燃料でも最後のビットの違いで止まりうる（今の炉心は
  制御棒の軸を持つ 1 つのテーブルなので起きない。確かめていない）
- `inherit` のミクロ燃焼の有無の照合の最後の分岐は、`nuclides` があって `auxiliary` のない履歴（コードで組んだもの）に、ミクロ燃焼でないという
  誤った文を出す（ファイルからは `CoreHistory::from_msgpack` が先に止める）
- テスト `reloads_carry_the_states_by_assembly_id` の補助の核種の確かめは、新しい集合体の最後の点の値が 0 でないことだけで、引き継いだ集合体の補助の
  核種と新しい集合体の燃焼点 0 の補助の核種は確かめていない（同じ装荷の restart は `restarts_continue_the_depletion_of_one_case` と
  `restarts_without_the_auxiliary_nuclides_would_not_continue` が確かめる）。`restarts_continue_the_depletion_of_one_case` は点の数を EFPD の列で確かめ、
  比べる繰り返しは点の数を確かめない
- 記録の表 `depletion-small-core-restart` は切る点を EFPD だけで示し、名前と日数を出さない（b の 182.2442 EFPD は停止の後）
- 前のサイクルの最後の点と次のサイクルの最初のステップの状態か出力が違う restart（サイクルの途中で切る）は、続けて解いた燃焼と違う燃焼点 0 になる
  （ADR-0220 の帰結）。エラーにも警告にもしていない

3 周目のレビューで挙がり、直さなかったもの:

- テスト `restart_mismatches_are_errors` の「数密度を持つノードが違う」の組は、追う核種の行だけを空にするので `MicroChain::compose` の数の違いで止まり、
  `deplete` のノードの照合（ミクロ燃焼するノードと履歴の数密度を持つノードの違い）の分岐に届かない
- examples/depletion/compare_core_restart.py の `if not count` は届かない分岐で、`relative_density_differences` は 1 か所からだけ呼ぶラッパー
- small_core_cases.py が書く b の前と後の半分のコメントは「停止 の点で切った」で、18 日の停止の後の点と読めない

関係する ADR: ADR-0220

## 完了の条件

Problem 10 の入力（仕様書の付録 K の装荷図）を書いて、上のそれぞれが要るかを決め、要るものを入れる。

## 経緯

- 2026-10-30: ADR-0220 の実装（ブランチ feat/core-restart-reload）の残りとして起こした
