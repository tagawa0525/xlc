---
status: open
labels: [verification]
opened: 2026-10-12
closed:
milestone:
blocked-by: []
resolution:
---

# rust-analyzer の負荷が時間の測定に重なることと、worktree に .envrc がないこと

## 問題

issue #0079 は build の置き場（ADR-0147）で閉じた。残る 2 点がある

- rust-analyzer は Slurm の外で動くので、索引と flycheck の CPU とメモリが、Slurm の排他のジョブ（時間を記録する測定）に
  重なると、記録に表れずに時間が延びうる（候補。数値は取っていない）
- 親のセッションから worktree に起動したエージェントは親の PATH を引き継ぐので、そのブランチでなく main の devShell の
  rust-analyzer と sysroot を使う。worktree には `.envrc` がない（git で追跡せず、`.worktreelink` にもない）ので、
  worktree で起動した Claude Code も devShell に入らない。`.envrc` を `.worktreelink` に足すのが候補
  （ADR-0147 の wrapper は git の最上位を実行時に見るので、どの worktree から起動しても、その worktree の
  `target/rust-analyzer` を使う）

関係する ADR: [ADR-0147](../adr/0147-separate-rust-analyzer-target-dir.md)

## 完了の条件

時間を記録する測定の間に rust-analyzer を止めるか、止めずに済む根拠を数値で示す。worktree で起動した Claude Code が
その worktree の devShell に入る。

## 経緯

- 2026-10-12: #0079 を ADR-0147 で閉じるときに分けた
