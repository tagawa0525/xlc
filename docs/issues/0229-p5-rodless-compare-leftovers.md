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
  解いた Slurm 694 が 29 分 28 秒。割り込まれうるパーティション main のジョブの経過時間なので、目安）を足す。炉心はその出力を読まない。1 周目のレビューは、`dvc add` で固めた棒のない出力を比べる形
  （照合の先が古い bin のまま残り、分岐の格子点を変えると軸の違いで `dvc repro` が止まる）を誤りとし、ステージにした。既定の
  `dvc repro -P` から外す形か、照合を 1 度だけ記録して Rust のテスト（`nodes_without_replaced_cells_hold_the_rodless_states`）に任せる形かは
  決めていない
- `compare_rodless.py` の `md5_inputs` と `leaves` は、examples/vera/compare_p5.py の `md5_inputs` と
  examples/depletion/warm-start/compare_warm_start.py の `leaves` と同じ形の複製。入力のファイルを `read` と `md5_inputs` で 2 回読む
- `compare_rodless.py` の入力の端: 2 つのファイルがどちらも `format` を持たないと KeyError で止まる（`unknown format` の SystemExit に
  ならない）。k∞ が 0 の状態で 0 で割る。文字列に直すと同じになるキー（整数の 1 と文字列の "1"）を持つ map は、片方を比べない
- examples/vera/p5/p5_inputs.py の `hzp_pipeline` は、棒のない出力のパスを 2 回求める。支持格子の層と反射体の格子の棒のない格子の
  ステージも、媒質を書かないのに `rodless/media` を作る（制御棒の分岐を解く格子の `mkdir -p lattice media` も同じ）
- examples/verification_tables.py の `TABLES` の JSON の置き場を確かめるテストが無い。`GIT_TRACKED` に無い JSON は、同じディレクトリの
  dvc.lock から md5 を読む（`locked_md5`。dvc.lock が無ければ `<dvc.lock> does not exist`、あってもその JSON の名前の `- path` の行が
  1 つでなければ `<dvc.lock> has N outputs named <JSON> (expected 1)` の `TableError`）ので、同じディレクトリの dvc.yaml がその JSON を
  出力に宣言しているかを確かめるテストを足す。
  3fd28562f で、置き場の誤り（examples/vera/p5/hzp/compare-rodless.json）が、記録用の計算（Slurm 695、696）を流し終えた後の記録を書く
  段で `examples/vera/p5/hzp/dvc.lock does not exist` として分かり、最後のレビューの周より後の挙動の変更になった。テストがあれば、
  `TABLES` に足したコミットの時点で止まる。急がない
- `locked_md5` のエラーの文言の残り: 親の dvc.yaml が下のディレクトリの JSON を出力に宣言した置き場の誤りでは、`dvc.lock does not exist`
  の「run the stage and commit dvc.lock」は当たらない（実行し直しても同じエラー）。`has N outputs named` は、まだ実行していない
  ステージ（N = 0）に実行の指示を示さない。`GIT_TRACKED_NOTE` の「its comment names」は、`GIT_TRACKED` の冒頭のコメント（AWS、燃焼度の軸）
  だけと読め、行ごとのコメントで足した試行の JSON と合わない。`GIT_TRACKED_NOTE` は `GIT_TRACKED` の定義から離れている
- `locked_md5` の正規表現は `deps` の行にも当たる。同じ JSON を出力するステージと deps に持つステージがあると N = 2 で止まり、deps の行だけが
  1 つあるとその md5 で照合を通す

関係する ADR: なし

## 完了の条件

各項目を直すかを決め、直すものを直す。

## 検討して捨てた形

- どの `locked_md5` のエラーにも「ステージの無い JSON は `GIT_TRACKED` に入れる」を示す形: 置き場を誤ったステージの出力を `GIT_TRACKED` に
  入れると md5 の照合（ADR-0176）が外れる（examples/vera/p5/hzp/compare-rodless.json の置き場の誤りが反例。レビューの指摘）

## 経緯

- 2026-11-05: chore/0213-rerun-p5-hzp のレビューの指摘のうち、通常の入力と流れで落ちるテストか誤った数値を示せないものをまとめた
- 2026-11-06: `TABLES` の JSON の置き場を確かめるテストの項を足した（3fd28562f の置き場の直しから。/retro の提案）
- 2026-11-07: `locked_md5` と未追跡のエラーの文言を、ステージの側の修正方法（実行、出力の宣言、`cache: false`）を先に示し `GIT_TRACKED` を例外の置き場として示す形に
  直した（担当 10 本が `dvc.lock does not exist` で止まっていた。/retro の提案）。文言の残りと deps の行の項を足した
