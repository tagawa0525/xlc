---
status: open
labels: [lattice, sensitivity, reproducibility]
opened: 2026-11-12
closed:
milestone: M4
blocked-by: []
resolution:
---

# 一般化随伴の解が同じ実行ファイルの実行の間でビットで揃わない

## 問題

- 2A の一般化随伴の 10 本（examples/depletion/issue-0279/2A-60-generalized-cmfd-gmres.toml、CMFD の反復の GMRES、判定 1e-7）を
  同じ実行ファイル（ブランチ perf/0281-moc-sweep-zeroing の b4bffdbd0 のリリースのビルド、分岐点 0d91e4b2b）で 2 回流すと、感度の
  出力の `orthogonality`（⟨Γ, Fφ⟩ / ⟨|Γ|, Fφ⟩）が 1e-15 の桁で違う（2026-11-13、Slurm 805 と 810）。テーブルの md5、k と k の感度、
  スイープの回数は同じ。main 0d91e4b2b のリリースのビルドと perf 用のビルドの間（Slurm 801 と 802）も同じ桁で違う。このため #0281 では
  一般化随伴の解のビットの一致を実行の比較で示せず、スイープと積の単体のテスト（`sweep_does_not_depend_on_the_contents_left_in_its_work_arrays`、
  `products_do_not_depend_on_the_contents_left_in_the_pool`）で示した
- 原因の候補（確かめていない）: `krylov::inner` は長さ 2^16 以上で rayon の `par_iter().sum()` を使う。rayon の和の分け方は
  仕事の盗み合いで変わりうるので、浮動小数点の和の順が実行ごとに変わりうる。ほかの並列の和は数えていない
- 固有値計算はスレッドの数によらない（`result_does_not_depend_on_the_thread_count` などのテスト）が、一般化随伴では成り立って
  いない。決まった長さの塊ごとに和を取り、塊の順に足せば、スレッドの数と盗み合いによらない

関係する ADR: なし

## 完了の条件

違いを生む和を確かめ、GMRES の内積（と、ほかの並列の和）を分け方によらない形にし、2 回の実行とスレッドの数を変えた実行で
一般化随伴の解がビットで同じことをテストで示す。

## 経緯

- 2026-11-12: #0281 の直しの前後を比べたときに、直しの前どうし（別のビルド）でも `orthogonality` が違うことに気づいた
- 2026-11-13: 同じ実行ファイルの 2 回の実行でも違うことを確かめた
