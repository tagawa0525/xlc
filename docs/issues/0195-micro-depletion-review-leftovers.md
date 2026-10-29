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

- 収率の境界より下の核分裂の割合 `thermal_fission` を、状態と燃焼度の重み（負の重みを含む）で割合のまま内挿し、0..1 を外れたらエラーにする
  （crates/core/src/state.rs の `interpolate_micro`）。割合は核分裂率で重みを付けた量ではなく、1 に近い核種は 3 次 Hermite や cut-HDMR の
  負の重みで 1 を超えうる。境界より下の核分裂の断面積をテーブルに持ち、ノードの中性子束で割合を作れば範囲の外に出ない（形式の項目を足すので、
  本番のテーブルの前に決める。ADR-0033、#0153）
- 重複した計算: `deplete_nodes`（crates/core/src/depletion.rs）は、解の中で線形の補正に使った微視的断面積を、同じ状態と燃焼度でノードごとに
  もう一度内挿する。`MicroNodes::new` は、テーブルのファイルごとに `MicroChain`（チェーン全体の崩壊の逆向きのクロージャを含む）を作り、
  最初の 1 つだけを使う（ほかは反応の並びの検査のため）。ノードの 1 回の燃焼の約 0.7 ms（ADR-0216 の確認 (d)）に内挿の 2 回目が含まれる
- 格子計算の燃焼との複製: 収率の点の選び方（`yield_weights`）、熱と高速の点、境界、`ELECTRON_VOLT`、`BARN`、負の数密度の許容は
  crates/lattice/src/depletion.rs と同じ値を core に複製した（core のブランチは lattice を変えない範囲）。片方だけを直すと、公称の条件の
  再現（ADR-0216 の確認 (b)）が崩れる。xlc-depletion に置けば 1 つになる（lattice の実行ファイルが変わるので、ステージの付け直しの費用を数える）
- コミット ca8d6eff の本文は、2A のミクロ燃焼の履歴が変わらないと書いたが、丸めの水準で変わった（訂正は 126c3215 の本文と記録）

関係する ADR: ADR-0216

## 完了の条件

- 各項目を直すか、直さない理由を書いて閉じる

## 経緯

- 2026-10-29: レビューの 3 周目の後に残した
