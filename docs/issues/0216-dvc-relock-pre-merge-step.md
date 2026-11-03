---
status: open
labels: [tooling, dvc]
opened: 2026-11-03
closed:
milestone:
blocked-by: []
resolution:
---

# dvc-relock のマージの直前の段を、手順の入口と許可に通す

## 問題

docs/dvc-relock-removed-stage-outs は、CLAUDE.md の作業の流れに dvc-relock の「マージの前に」（main の作業ツリーに残る古い
出力を消す段）を足した。その 1 周目の /code-review（2026-11-03、HEAD b83289912）が挙げた次の点は、手順のとおりに実行して
失敗する誤りではないので直していない。

- /topic-start（委任を受けた担当も読む）は dvc-relock の後半だけを挙げ、マージの直前の段に触れない。同じスキルの「作業中は
  main の作業ツリーに触らない」とも、例外を書かずに食い違う
- dvc-relock の入口の文（前半と後半を worktree とマージ済みかで選ぶ）とフロントマターの description が、マージの直前の段を
  含まない。「マージの前に」の節は `## 前半` の下の `###` で、「この後は /git-merge の「チェック」に進む。」の後に置かれている
- dvc-relock の allowed-tools に、main の worktree での `git -C … status`、`ls`、`rm` が無い

## 完了の条件

- /topic-start とスキルの入口からマージの直前の段にたどれ、その段の操作が許可の確認なしに実行できる

関係する ADR: なし
