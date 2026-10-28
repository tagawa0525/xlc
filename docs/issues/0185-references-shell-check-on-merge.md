---
status: open
labels: [dvc, tooling]
opened: 2026-10-28
closed:
milestone:
blocked-by: []
resolution:
---

# references-shell の検査をマージの経路で実行する

## 問題

`scripts/check-references-shell.sh`（ADR-0212）は、pre-commit の hook `references-shell` としてコミットのたびに実行される。
rebase とマージはコミットの hook を通らないので、ブランチと main がそれぞれ hook を通した flake の変更が rebase で混ざると、
混ざった flake の `references` の devShell の store のパスは、どちらの `references-shell.txt` とも違いうる。/dvc-relock は同じ
スクリプトを手で実行するが、`dvc.yaml` の無いブランチは /dvc-relock を通らず、/git-merge の pre-merge のチェック
（`run-checks.sh --merge`）はこのスクリプトを呼ばない。#0183 と同じく、食い違いが `dvc status` に出ないまま main に入りうる。

`examples/verification_tables.py --check` は hook のほかに `cargo test`（crates/docs の検査）からも実行され、マージの前のチェックで
止まる。同じ形にするには crates/docs に検査を足す（#0183 のブランチは crates/ を変えない範囲だった）。

関係する ADR: なし

## 完了の条件

rebase とマージで混ざった flake でも、マージの前のチェックが `references-shell.txt` の食い違いで止まる。

## 経緯

- 2026-10-28: #0183 のブランチ（fix/0183-references-shell）のレビューで見つかった
