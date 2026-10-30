---
status: open
labels: [tooling, dvc]
opened: 2026-10-31
closed:
milestone:
blocked-by: []
resolution:
---

# 記録が引く追跡しない生成物に .dvc があるかを検査する

## 問題

ADR-0222 は、記録（examples/*/verification*.md）と ADR が md5 か名前で引く追跡しない生成物を `dvc add` で追うと決めたが、
新しい生成物を追うのは、記録を書く担当が `dvc add` を覚えていることに頼る（CLAUDE.md の「計算」の 1 行）。忘れると、
/git-merge の後片付けが #0207 と同じく黙って消す。記録の md5 と `.dvc` ファイルの md5 を突き合わせる検査は無い。

あわせて、ADR-0222 のブランチのレビュー（2026-10-31）で、次の小さい点が残った:

- `scripts/dvc-checkout.py` は、作り直して `dvc add` をやり直していない生成物で dvc が出す「Can't remove the following
  unsaved files」を読まず、「戻らなかった出力の一覧が無い」の一般の誤りで止まる（dvc の出力はその前に表示される）
- `scripts/count-changed-stages.py` の `total:` の `stages` は、コマンドを持たない `*.dvc` の行も足す（生成物が無いか
  変わった worktree で、走るステージの数の見積もりより大きく出る）

関係する ADR: ADR-0222

## 完了の条件

記録か ADR に新しく書いた md5 か名前が、追跡しない生成物を指すのに `.dvc` が無いか md5 が合わないとき、チェック
（pre-commit の hook か `run-checks.sh` の中の検査）が止める。#0207 の /retro の試し（md5 の拾い方、DVC の出力の除き方）の
指摘を検査の側で扱う。上の 2 点は、検査を足すときに合わせて直す。
