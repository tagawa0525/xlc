---
status: open
labels: [core, lattice, depletion]
opened: 2026-11-07
closed:
milestone: M4
blocked-by: []
resolution:
---

# ノードの出力の重み（#0154）のレビューの残り

## 問題

ブランチ feat/0154-node-power-by-fission-energy のレビューの、直さなかった指摘（再現を書けないもの、重複と効率）:

- crates/core/src/case.rs の `node_constants`: 内挿した Σf が 0 以下の群は核分裂あたりのエネルギーを 0 とし、補正で Σf が正になっても κΣf は 0 の
  まま。補正した κΣf の符号も確かめない（通常の入力で内挿した Σf が 0 以下になる例は見つけていない）
- crates/lattice/src/history.rs の `history_branches`: チェーンと履歴の照合は核種の名前と順だけで、同じ核種の並びで Q の違うチェーン（作り直した
  チェーン）を通す。履歴にチェーンの Q か識別を持たせれば止められる
- crates/lattice/src/history.rs: `kappa_fission` は `micro` と同じ SPH 因子と少数群の中性子束を状態ごとに求め直す（重複と計算の繰り返し）
- crates/core/src/nem.rs: `fission_powers` は `fission_rates` のループのコピー（断面積とエラーだけが違う）
- crates/core/src/rods.rs の `mix`: κΣf の有無を数えてから NaN を返す取り出しで混ぜる（すべてか、どれも持たないかの検査を先に書けば NaN が要らない）
- crates/core/src/depletion.rs: 燃焼の前の κΣf の検査がすべての状態を見る（テーブルの検査が、すべての状態が持つか、どれも持たないかを保証する）
- examples/depletion/compare_core_energy.py: 集合体の燃焼度を重金属の質量から求め直す（履歴の点の `assembly_burnups` と同じ定義）

関係する ADR: なし

## 完了の条件

点ごとに、直すか直さないかを決めて、直すものを直す。

## 経緯

- 2026-11-07: #0154 のブランチの 1 周目のレビューで出た、直さなかった指摘を集めた
