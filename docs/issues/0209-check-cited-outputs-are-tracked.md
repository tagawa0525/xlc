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

あわせて、ADR-0222 のブランチのレビュー（2026-10-31、3 周）で、次の点が残った:

- `scripts/dvc-checkout.py` は、作り直して `dvc add` をやり直していない生成物で dvc が出す「Can't remove the following
  unsaved files」を読まず、「戻らなかった出力の一覧が無い」の一般の誤りで止まる（dvc の出力はその前に表示される）
- `scripts/count-changed-stages.py` の `total:` の `stages` は、コマンドを持たない `*.dvc` の行も足す（生成物が無いか
  変わった worktree で、走るステージの数の見積もりより大きく出る）

- 作り直して `dvc add` をやり直していない生成物のある進行中のブランチは、ADR-0222 のマージに rebase した後の
  `worktree-restore-outputs.sh` とスキル dvc-relock の前半 1 で checkout 全体が止まる（対処はスキル dvc-relock の前半 1。上の
  一般の誤りの文で原因が見えにくい）
- examples/depletion/verification.md の「分岐のテーブルの交差の項」は、`table-2B.msgpack` と
  `table-2A-fine.msgpack` を main の作業ツリーからコピーすると書く（今は `.dvc` で戻る。記録の本文はこのブランチで変えない）
- スキル topic-start の「生成物を戻す」は、終了コード 255 と、戻らない出力を `dvc repro` で作ると書く（スクリプトは
  dvc-checkout.py の 0 か 1 を返し、`.dvc` の対象は `dvc repro` で作れない）
- CLAUDE.md の `dvc repro -P` と `dvc repro <yaml>` は、変わった `.dvc` の対象の md5 を書き直しうる（ADR-0222 の帰結）が、
  その注意は CLAUDE.md に無い
- スキル dvc-relock の後半の `.dvc` の checkout は `HEAD^1 HEAD` で、マージのコミットが HEAD である前提（後にほかのマージや
  relock のコミットが載ると、別の差分を拾う）
- スキル verification-record に、記録に書く前の `dvc add` の手順が無い（CLAUDE.md の「計算」にだけある）

関係する ADR: ADR-0222

## 完了の条件

記録か ADR に新しく書いた md5 か名前が、追跡しない生成物を指すのに `.dvc` が無いか md5 が合わないとき、チェック
（pre-commit の hook か `run-checks.sh` の中の検査）が止める。#0207 の /retro の試し（md5 の拾い方、DVC の出力の除き方）の
指摘を検査の側で扱う。上の点は、検査を足すときに合わせて直す。
