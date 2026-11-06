#!/usr/bin/env bash
# pre-commit の hook `docs-consistency` が実行する、文書の整合性の検査（crates/docs/tests/consistency.rs の
# `docs_are_consistent`。ADR-0233）。prek がステージしていない変更を退避した作業ツリーと、git が渡す
# GIT_INDEX_FILE の index（`git ls-files`）を読むので、ステージした中身を見る。
#
# 排他の Slurm のジョブの有無によらず、いつも 1 コアに限る。hook の中の cargo は claude-hooks の
# guard-exclusive-slurm に見えず、ジョブの時間の記録を乱しうる。xlc-docs は依存のクレートを持たず、1 コアでも
# 初回のビルドは数秒、キャッシュがあれば 0.1 秒ほどなので、squeue で判定するより常に限る方が単純で
# 失敗の道も少ない。CPU は、固定の番号でなく、このプロセスに許された CPU の最初の 1 つ（Slurm の
# `-c` のジョブの中のように許された CPU が限られていると、外の番号では taskset が失敗する）
set -euo pipefail

cd -- "$(git rev-parse --show-toplevel)"

allowed=$(sed -n 's/^Cpus_allowed_list:[[:space:]]*//p' /proc/self/status)
if [[ -z "$allowed" ]]; then
  echo "check-docs-consistency.sh: /proc/self/status から許された CPU を読めません" >&2
  exit 1
fi

exec taskset -c "${allowed%%[-,]*}" \
  cargo test -q -p xlc-docs --test consistency -- --exact docs_are_consistent
