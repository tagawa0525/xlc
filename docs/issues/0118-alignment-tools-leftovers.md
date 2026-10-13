---
status: open
labels: [performance, build]
opened: 2026-10-14
closed:
milestone:
blocked-by: []
resolution:
---

# ADR-0157 の測定の道具のレビューで残った指摘

## 問題

ブランチ perf/0105-code-alignment のレビュー（7 周）で、直さずに残した指摘。ADR-0157 の決定は変えないが、道具を
使い回すとき（ADR-0144 の判定など）に効く。

- examples/build-alignment/build.sh は、scripts/build-binaries.sh の確認（環境変数の rustflags、`NIX_NO_SELF_RPATH`）を
  写している。build.sh は `grep -q` を直したが、build-binaries.sh はまだ（#0113）。build.sh は記号を落とさず
  （`strip="symbols"` なし）、実行ファイルの md5 と大きさは `build` の段階のものと合わない
- `measure-xs.sh`、`measure-lattice-core.sh`、`measure-repeat.sh` は、examples/m0 の段階を DVC を通さずに走らせ、その
  出力の置き場所（`wims69.msgpack`、`pincell-fg.msgpack`）に書く。あれば止め、終わりに消すが、SIGKILL や Slurm の時間切れ
  では残りうる
- `measure-2A.sh` は各回で v0 を先に走らせ、順の影響を分けていない（`measure-repeat.sh` は入れ替える）
- #0112 の表の 284、287、292 は、コミットしていない手元のスクリプトとビルドで出した値で、再現できない
- 置き場所のゆらぎを測り方の側で扱う案（置き場所を変えた複数のビルドで判定する、命令数のカウンタで比べる）は比べて
  いない（ADR-0157 の「決定と理由」）

関係する ADR: なし

## 完了の条件

ADR-0144 の判定の手順を決めるとき（#0112）に、上のそれぞれを直すか、要らないと判断する。

## 経緯

- 2026-10-14: perf/0105-code-alignment のレビューで残った指摘をまとめた
