---
status: open
labels: [examples, vera]
opened: 2026-11-05
closed:
milestone: M3
blocked-by: []
resolution:
---

# p5 の棒のない状態の照合のレビューの残り（棒のない格子の常設の費用、比較のツールの重複と入力の端）

## 問題

ブランチ chore/0213-rerun-p5-hzp（#0225 の照合。examples/vera/p5 の `hzp-rodless-*` と `hzp-compare-rodless`、`compare_rodless.py`）の
レビュー（/code-review high の 3 周、2026-11-05）で出た、直さずに残した指摘:

- 棒のない格子 8 本は p5 の常設のステージで、`bin/lattice` か `vera.msgpack` が変わるたびに、p5 の再実行に約 30 分（棒のない 8 本を
  解いた Slurm 694 が 29 分 28 秒）を足す。炉心はその出力を読まない。1 周目のレビューは、`dvc add` で固めた棒のない出力を比べる形
  （照合の先が古い bin のまま残り、分岐の格子点を変えると軸の違いで `dvc repro` が止まる）を誤りとし、ステージにした。既定の
  `dvc repro -P` から外す形か、照合を 1 度だけ記録して Rust のテスト（`nodes_without_replaced_cells_hold_the_rodless_states`）に任せる形かは
  決めていない
- `compare_rodless.py` の `md5_inputs` と `leaves` は、examples/vera/compare_p5.py の `md5_inputs` と
  examples/depletion/warm-start/compare_warm_start.py の `leaves` と同じ形の複製。入力のファイルを `read` と `md5_inputs` で 2 回読む
- `compare_rodless.py` の入力の端: 2 つのファイルがどちらも `format` を持たないと KeyError で止まる（`unknown format` の SystemExit に
  ならない）。k∞ が 0 の状態で 0 で割る。文字列に直すと同じになるキー（整数の 1 と文字列の "1"）を持つ map は、片方を比べない
- examples/vera/p5/p5_inputs.py の `hzp_pipeline` は、棒のない出力のパスを 2 回求める。支持格子の層と反射体の格子の棒のない格子の
  ステージも、媒質を書かないのに `rodless/media` を作る（制御棒の分岐を解く格子の `mkdir -p lattice media` も同じ）

関係する ADR: なし

## 完了の条件

各項目を直すかを決め、直すものを直す。

## 経緯

- 2026-11-05: chore/0213-rerun-p5-hzp のレビューの指摘のうち、通常の入力と流れで落ちるテストか誤った数値を示せないものをまとめた
