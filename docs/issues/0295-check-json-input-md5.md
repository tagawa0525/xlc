---
status: open
labels: [dvc, checks]
opened: 2026-11-19
closed:
milestone:
blocked-by: []
resolution:
---

# 追跡している比較の JSON の入力の md5 が、今の入力と合うかを機械で確かめる

## 問題

比較のステージの JSON（例: examples/depletion/p9-fuel-points/fuel-points-2A.json）は、読んだ入力のパスと md5 を `inputs` に
書き、git で追跡している（`cache: false`）。入力の TOML のコメントだけを直すと、ステージは changed になり、JSON の md5 は古い
まま残るが、それを止める検査が無い。#0291 の語の直し（ブランチ chore/0291-terms-in-example-inputs）では、1 度目の再実行で
p9-2A、p9-product、p9-gd、examples/core の比較 9 本を流し、fuel-points-2A、fuel-points-2O（p9-fuel-points）、hzp-compare（p5）、
compare-6（p6）の 4 本を落とした（レビューで見つかり、後のコミットで再実行した）。次のコメントだけの直しでも同じ落としが起きうる。

再現: 追跡している JSON の `inputs` のパスごとに今のファイルの md5 を計算し、JSON の md5 と比べる。

関係する ADR: なし

## 完了の条件

- 追跡している JSON の `inputs` の md5 が今のファイルと食い違えば止まる検査（hook か `cargo test`、`verification_tables.py --check`
  のどれか）を置くか、置かない理由を決める

## 経緯

- 2026-11-19: #0291 のレビューで、比較の 4 本の再実行の落としが見つかったので起こした
- 2026-11-19: 入力の TOML を正規形で deps にする作業（#0145、ADR-0212 の案 B5）で一緒に片付ける（親の判断）。正規形なら
  コメントだけの直しでステージも JSON の md5 も変わらない
