---
status: open
labels: [build, dvc]
opened: 2026-10-14
closed:
milestone:
blocked-by: []
resolution:
---

# build-binaries.sh の環境変数の確認が、pipefail の下の grep -q で素通りしうる

## 問題

scripts/build-binaries.sh（ADR-0159）は、コンパイルとリンクの設定が環境変数にあれば止める確認を
`if env | grep -qE "$flags"; then` で行う。スクリプトは `set -o pipefail` なので、`grep -q` が最初の一致で抜けたあとに
`env` が書き込むと SIGPIPE（終了コード 141）で終わり、パイプライン全体が失敗とみなされて、条件が偽になりうる。
そのときは、`RUSTFLAGS` などがあってもビルドに進み、DVC の deps に表れない設定で実行ファイルが作られる。起きる条件は
`env` の出力の大きさと一致の位置による（試していない）。

ブランチ perf/0105-code-alignment のレビューで、同じ形を写した examples/build-alignment/build.sh について指摘され、
そちらは `grep -E … >/dev/null` に直した。

関係する ADR: なし

## 完了の条件

scripts/build-binaries.sh の確認を、`-q` を使わない形（`grep -E … >/dev/null`、または `env` の出力を変数に受けてから
`grep`）にする。

## 経緯

- 2026-10-14: perf/0105-code-alignment のレビューで見つけた
