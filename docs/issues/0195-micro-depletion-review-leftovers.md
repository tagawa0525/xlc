---
status: open
labels: [core, depletion]
opened: 2026-10-29
closed:
milestone: M4
blocked-by: []
resolution:
---

# core のミクロ燃焼のレビューで残した指摘

## 問題

core のミクロ燃焼（ADR-0216、ブランチ feat/0153-core-micro-depletion）の /code-review を 3 周で打ち切り（ユーザーの決定 2026-10-29）、
次の指摘を残した。どれも 2A の集合体 1 体の結果を変えない。

- 収率の境界より下の核分裂の割合は、割合 × 群の和の σf と σf を同じ重みで内挿した比にした（c82f3a84 の前の ac60d199）。群の和の σf は
  ノードの中性子束の重みでなく、状態ごとの割合を作った格子の中性子束の重みとも違う。境界より下の核分裂の断面積を群ごとにテーブルに持てば、
  ノードの中性子束で割合を作れる（形式の項目を足すので、本番のテーブルの前に決める。ADR-0033、#0153）
- 重複した計算: `deplete_nodes`（crates/core/src/depletion.rs）は、解の中で線形の補正に使った微視的断面積を、同じ状態と燃焼度でノードごとに
  もう一度内挿する。`MicroNodes::new` は、テーブルのファイルごとに `MicroChain`（チェーン全体の崩壊の逆向きのクロージャを含む）を作り、
  最初の 1 つだけを使う（ほかは反応の並びの検査のため）。状態の重みの組み立て（`state_weights`）は、少数群定数の内挿、微視的断面積の内挿、
  `deplete_nodes` の微視的断面積の内挿で、同じ状態と燃焼度について 3 回行う。核分裂のあるテーブルかの判定（`fissile`）を位置ごとに評価する
- 格子計算の燃焼との複製: 収率の点の選び方（`yield_weights`）、熱と高速の点、境界、`ELECTRON_VOLT`、`BARN`、負の数密度の許容は
  crates/lattice/src/depletion.rs と同じ値を core に複製した（core のブランチは lattice を変えない範囲）。片方だけを直すと、公称の条件の
  再現（ADR-0216 の確認 (b)）が崩れる。xlc-depletion に置けば 1 つになる（lattice の実行ファイルが変わるので、ステージの付け直しの費用を数える）
- コミット ca8d6eff の本文は、2A のミクロ燃焼の履歴が変わらないと書いたが、丸めの水準で変わった（訂正は 126c3215 の本文と記録）

関係する ADR: ADR-0216

## 完了の条件

- 各項目を直すか、直さない理由を書いて閉じる

## 経緯

- 2026-10-29: レビューの 3 周目の後に残した
- 2026-10-29: 親の 2 周目（head 08f8aca79、/code-review high）で残した指摘（どれも Medium と判断し、周の上限で直さず残す）:
  - `solve_compositions`（pub）は、テーブルの `micro` の核種と `Compositions` のチェーンの核種の並びが合うことを確かめず、ずれれば
    位置で zip され、短ければ `tracked()` が index で panic する。並びと数を確かめてエラーにする
  - `interpolate_micro` の負の重みの検査が散乱行列を含まない（少数群定数の `non_negative()` は含む）。負に内挿した σ_s が補正に入る
  - `CoreHistory::from_msgpack` が、`micro` の核種が空か重複し、数密度が全ノードで空の履歴を受ける
  - `compare_core_micro.py` の `summarize` は、ミクロ燃焼の無いノードで KeyError になる（SystemExit で原因を言うべき）
  - `rates()` は `fission` の群数しか中性子束と比べず、`total`、`nu_fission`、反応の経路は zip で黙って切れる
  - σ_f の群の和が 0 以下の核分裂性の核種は `thermal_fission` を nil にし、`rates()` が核分裂しない核種として扱う（黙ったフォールバック）
  - `deplete_nodes` がノードごとに `solve_with` と同じ σᵢ の内挿を繰り返す（予測子と修正子で 2 回。上の重複の項と同じ）
  - `MicroNodes::new` が核分裂のあるテーブルのファイルごとに `MicroChain` 全体（3,820 核種の親のクロージャ）を作り、最初のものだけ使う
