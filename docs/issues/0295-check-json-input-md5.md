---
status: closed
labels: [dvc, checks]
opened: 2026-11-19
closed: 2026-11-19
milestone:
blocked-by: []
resolution: ADR-0254
---

# 追跡している比較の JSON の入力の md5 が、今の入力と合うかを機械で確かめる

## 問題

比較のステージの JSON（例: examples/depletion/p9-fuel-points/fuel-points-2A.json）は、読んだ入力のパスと md5 を `inputs` に
書き、git で追跡している（`cache: false`）。入力の TOML のコメントだけを直すと、ステージは changed になり、JSON の md5 は古い
まま残るが、それを止める検査が無い。#0291 の語の直し（ブランチ chore/0291-terms-in-example-inputs）では、1 度目の再実行で
p9-2A、p9-product、p9-gd、examples/core の比較 9 本を流し、fuel-points-2A、fuel-points-2O（p9-fuel-points）、hzp-compare（p5）、
compare-6（p6）の 4 本を落とした（レビューで見つかり、後のコミットで再実行した）。次のコメントだけの直しでも同じ落としが起きうる。

再現: 追跡している JSON の `inputs` のパスごとに今のファイルの md5 を計算し、JSON の md5 と比べる。

関係する ADR: ADR-0254（JSON の TOML の入力を正規形の md5 にし、lock の deps と照らす）

## 完了の条件

- 追跡している JSON の `inputs` の md5 が今のファイルと食い違えば止まる検査（hook か `cargo test`、`verification_tables.py --check`
  のどれか）を置くか、置かない理由を決める

## 経緯

- 2026-11-19: #0291 のレビューで、比較の 4 本の再実行の落としが見つかったので起こした
- 2026-11-19: 入力の TOML を正規形で deps にする作業（#0145、ADR-0212 の案 B5）で一緒に片付ける（親の判断）。正規形なら
  コメントだけの直しでステージも JSON の md5 も変わらない
- 2026-11-19: ADR-0254 で決めて実装した。JSON の `inputs` の TOML は値の正規形の md5 にし（コメントだけの修正で古くならない）、
  `scripts/test_normalize_toml.py` の `JsonInputsTest` が、追跡している比較の JSON の `inputs` と、それを書いたステージの dvc.lock の
  deps（TOML は正規形）の md5 を照らす（pre-commit の hook `toml-values`）。今のファイルとの比較（この issue の書いた形）は採らない:
  上流の値を変えて再実行するまでの間（ADR-0130 で changed のまま残してよい）も止まり、JSON と lock が一致すれば、今のファイルとの
  食い違いはステージの changed として `dvc status` に出る。JSON の md5 を変えた 11 本を再実行した（565aec357）。照らすのはステージの
  出力の JSON だけで、手で実行した試行の JSON（examples/depletion/branches/ など）は、照らす lock が無いので見ない（ADR-0254 の帰結）
