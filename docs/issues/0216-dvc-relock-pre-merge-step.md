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

b83289912 は、CLAUDE.md の作業の流れに dvc-relock の「マージが止まったら」（main の作業ツリーに残る古い
出力を消す段）を足した。その /code-review（2026-11-03 の 1〜3 周目）が挙げた次の点は、手順のとおりに実行して
失敗する誤りではないので直していない。

- /topic-start（委任を受けた担当も読む）は dvc-relock の後半だけを挙げ、マージの直前の段に触れない。同じスキルの「作業中は
  main の作業ツリーに触らない」とも、例外を書かずに食い違う
- dvc-relock の入口の文（前半と後半を worktree とマージ済みかで選ぶ）とフロントマターの description が、マージの直前の段を
  含まない。「マージが止まったら」の節は `## 前半` の下の `###` で、「この後は /git-merge の「チェック」に進む。」の後に置かれている
- dvc-relock の allowed-tools に、main の worktree での `git -C … status`、`ls`、`rm` が無い
- dvc.yaml を手で直して `.gitignore` の行が残ったとき、持ち主の無い出力は無視されたまま main に残る（マージは止まらない）
- `cache: false` に変えた出力で、main の無視された版がブランチの追跡する版と違っても、git は黙って上書きする（前の手順の
  「差があるファイルは消さず原因を調べる」は、上書きが止まる前提で書かれていた）

関係する ADR: なし

## 完了の条件

- /topic-start とスキルの入口からマージの直前の段にたどれ、その段の操作が許可の確認なしに実行できる
- 上の残りの 2 点を、直すか、直さない理由を書いて閉じる

## 経緯

- 2026-11-03: #0093 のマージ（dcc9036a2）が、examples/core の持ち主の無くなった古い出力 6 本で pre-merge の hook に止まった。
  /retro の docs/dvc-relock-removed-stage-outs で dvc-relock の節を hook の一覧を使う形に直し、残りをここに書いた
