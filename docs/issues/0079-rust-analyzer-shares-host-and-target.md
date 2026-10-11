---
status: open
labels: [verification]
opened: 2026-10-12
closed:
milestone:
blocked-by: []
resolution:
---

# rust-analyzer が計算機と target/debug を DVC の段階や測定と共有する

## 問題

ADR-0135 で devShell に rust-analyzer を入れたので、devShell の中で起動した Claude Code は、Slurm の計算機（r995）の上で
Slurm の外から rust-analyzer を走らせる。

- ワークスペースの索引と、編集のたびの flycheck（`cargo check --workspace --all-targets`）が CPU とメモリを使う。
  Slurm の排他のジョブはこれを止めないので、時間を記録する測定に重なると、記録に表れずに時間が延びうる
- flycheck は `target/debug` を run-checks.sh のテストや clippy と共有する（2026-10-12 に worktree の
  target/flycheck0/stdout の成果物のパスで確認）。clippy は rustc を clippy-driver で包むので、互いのチェックの結果を
  無効にして作り直し合うことと、`target/debug` のロックを取り合うことが候補（確かめていない）
- リポジトリの `rust-analyzer.toml` に `[cargo] targetDir = true` を書く形は、rust-analyzer 1.99.0 では効かなかった
  （`target/rust-analyzer` を作らず、`target/debug` を使い続けた）
- 親のセッションから worktree に起動したエージェントは親の PATH を引き継ぐので、そのブランチでなく main の devShell の
  rust-analyzer と sysroot を使う

関係する ADR: [ADR-0135](../adr/0135-keep-dev-tools-out-of-rust-toolchain.md)

## 完了の条件

flycheck の置き場を `target/` の外の専用のディレクトリにし（lsp-det の `.lsp.json` の initializationOptions の
`cargo.targetDir` などが候補）、時間を記録する測定の間に rust-analyzer を止めるか、止めずに済む根拠を数値で示す。

## 経緯

- 2026-10-12: ADR-0135 のレビューで挙がった
