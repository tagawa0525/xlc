---
status: accepted
date: 2026-10-12
supersedes: []
superseded-by:
issues: ["#0079", "#0094"]
requires: [ADR-0135]
---

# ADR-0147: rust-analyzer の build を target/rust-analyzer に分け、worktree ごとの target は共有しない

## 背景

devShell の rust-analyzer（[ADR-0135](0135-keep-dev-tools-out-of-rust-toolchain.md)）の flycheck（`cargo check`）は、
`target/debug` を run-checks.sh の clippy とテスト、手で走らせる `cargo build` と共有していた（#0079）。cargo は
ビルドディレクトリのロックを取るので、flycheck の間は計算機側のビルドが待たされ、clippy と check は同じ
`target/debug` の fingerprint を作り直し合うことが候補になる（確かめていない）。

あわせて、worktree は 12 個ほどあり、それぞれが `target/`（main で約 20 GB）を持って最初にフルビルドする。`target` を
worktree 間で共有すれば再利用できるかを同じ ADR で決める。判断の軸は、正確さ（互いのビルドを壊さない）、実行速度、
保守の容易さの順。設定は git で追跡できる場所に置き、個人の環境（~/.config、クライアントの設定）に頼らない。

## 検討した案

rust-analyzer の置き場:

- 案 A: リポジトリの `rust-analyzer.toml` に `[cargo] targetDir = true`。rust-analyzer 1.99.0 では効かなかった
  （#0079。`target/rust-analyzer` を作らず `target/debug` を使い続けた）。`cargo.*` は rust-analyzer.toml でなく
  クライアントの initializationOptions で受ける設定のため（候補。確かめていない）
- 案 B: クライアント（lsp-det の `.lsp.json`）の initializationOptions に `cargo.targetDir` を書く。リポジトリの外の
  設定で、git で追跡できず、クライアントごとに書き直しになる
- 案 C: `.cargo/config.toml` の `build.target-dir`、または devShell の `CARGO_TARGET_DIR`。cargo を使うもの全部
  （run-checks.sh、Slurm のジョブ、DVC の段階）の出力先が変わり、分けたい rust-analyzer だけを分けられない。
  ビルドの出力先は DVC の段階の環境にもなる（ADR-0135 の帰結）
- 案 D（採用）: devShell の `rust-analyzer` を、`CARGO_TARGET_DIR=<git の最上位>/target/rust-analyzer` を設定して
  本体を exec する wrapper にする。cargo は子プロセスとして環境を引き継ぐので、flycheck と build script がそこに出る

worktree 間の `target` の共有:

- 案 X: `CARGO_TARGET_DIR` を全 worktree で 1 つにする。ビルドディレクトリのロックがビルド全体に及ぶので、
  並行する worktree のビルド、clippy、テストが 1 本ずつに直列化する。時間を記録する測定（Slurm の排他のジョブ）が
  ほかの worktree のビルドを待つ・待たせることになる。ワークスペースの crate は相対パスで同じ名前になるので、
  ブランチが違うと互いの成果物を上書きし合う候補がある（確かめていない）
- 案 Y（採用）: 今のまま worktree ごとの `target/` にする

## 決定と理由

案 D と案 Y を採る。案 D は、設定が flake.nix にあって git で追跡でき、分けるのが rust-analyzer だけで、cargo を直接
使うものは今のまま `target/debug` を使う。`rust-analyzer.toml` もクライアントの設定も要らない。git の外で起動すると
`git rev-parse` が失敗して止まる（黙ったフォールバックを置かない）。`target/rust-analyzer` は `/target` の下なので
無視される。

案 Y を採るのは、正確さと速度の両方が worktree 間の分離に依るため。共有で得るのは registry の依存の再利用と
ディスク（空き 626 GB に対し 12 worktree で数百 GB 以下）だけで、その代わりに全ビルドが直列化する。registry の
依存の再利用が欲しければ、ロックを取り合わない別の仕組み（コンパイル結果のキャッシュ）を、必要になったときに別の
決定で足す。

## 確認

devShell の wrapper で確かめた（ブランチ feat/0079-target-dir の worktree、2026-10-12。flake.nix の変更と同じ
コミット）。

- LSP（initialize、initialized、`didSave`）を stdio で rust-analyzer に与えると、flycheck が 43 秒で終わり、
  `target/rust-analyzer/{debug,flycheck0}` ができ、`target/debug` はできなかった
- `cargo clippy --workspace --all-targets`（`target/debug`）と、`CARGO_TARGET_DIR=target/rust-analyzer` での
  `cargo check --workspace --all-targets --message-format=json`（flycheck と同じ形）を同時に走らせて、どちらも成功した。
  待ちの表示は `package cache` のロックだけ（~/.cargo の共有で、起動時の短い間。ビルドディレクトリのロックの待ちは
  なかった）
- 共有の案 X の下限: git archive で作った同一内容の 2 つの木で `CARGO_TARGET_DIR` を共有すると、2 つ目のビルドは
  何も Compiling しなかった（1 つ目は 19 件）。内容が同じなら再利用できるが、ブランチが違う場合の上書きと待ちは
  測っていない

## 帰結

- rust-analyzer の索引と flycheck の CPU とメモリの負荷は残る。時間を記録する測定への影響の扱いは #0094 で続ける
- 計算機側の cargo は今まで通り `target/debug`。flycheck が `target/debug` の fingerprint を壊すことはなくなる
- worktree ごとに `target/rust-analyzer` が増える（check だけの成果物なので `target/debug` より小さい）
- worktree-add.sh（nixfiles 側）に手を入れる必要はない。wrapper は git の最上位を実行時に見る
