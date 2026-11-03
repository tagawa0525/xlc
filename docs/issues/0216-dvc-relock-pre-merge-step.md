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

b83289912 は、CLAUDE.md の作業の流れに dvc-relock の「マージの前に」（main の作業ツリーに残る古い
出力を消す段。b33deca17 で「マージが止まったら」に改めた）を足した。その /code-review（2026-11-03 の 1〜3 周目）が挙げた次の点は、手順のとおりに実行して
失敗する誤りではないので直していない。

- /topic-start（委任を受けた担当も読む）は dvc-relock の後半だけを挙げ、マージの直前の段に触れない。同じスキルの「作業中は
  main の作業ツリーに触らない」とも、例外を書かずに食い違う
- dvc-relock の入口の文（前半と後半を worktree とマージ済みかで選ぶ）とフロントマターの description が、マージの直前の段を
  含まない。「マージが止まったら」の節は `## 前半` の下の `###` で、「この後は /git-merge の「チェック」に進む。」の後に置かれている
- dvc-relock の allowed-tools に、main の worktree での `git -C … status`、`ls`、`rm` が無い
- dvc.yaml を手で直して `.gitignore` の行が残ったとき、持ち主の無い出力は無視されたまま main に残る（マージは止まらない）
- `cache: false` に変えた出力で、main の無視された版がブランチの追跡する版と違っても、git は黙って上書きする（前の手順の
  「差があるファイルは消さず原因を調べる」は、上書きが止まる前提で書かれていた）

3 周の後の記録の周（2026-11-03、HEAD 9a911d605）が、直さずに残った blocking を挙げた（歯止めの外なので直していない。マージは
ユーザーの判断）:

- 手順 2 は「キャッシュから戻せる」を `dvc status` で見るが、`modified`（main の古い出力の中身が dvc.lock の md5 と違う。後半の節の
  とおり DVC を入れる前の出力が main に残りうる）の出力も消す対象に読める。消すと中身は git にもキャッシュにも無く、手順 3 は
  lock の版を置き、消したものは戻らない
- 手順 1 は hook のパスを `.gitignore` の消えた行と手で突き合わせるので、ディレクトリの outs（`-/outdir` に対して
  `outdir/a` などが並ぶ）で迷う。`git merge --abort` の後の `git status --short --ignored <パス>` で `!!`（main で無視されていた）
  かを見れば機械的に分かる
- 節は hook の一覧を「マージの後に未追跡になるパスそのもの」と書くが、hook は `git diff --name-only` と無視されない未追跡の
  全部を並べるので、前からある未追跡のファイルや変更も混じる
- 手順 2 の `dvc status` は、/topic-start の「main の作業ツリーでは dvc を実行しない」と食い違う
- /git-merge は hook の拒否に「abort、ブランチで直し、rebase とレビューとチェックからやり直す」を求め、この節に向かわない

捨てた方式: 前もって `.gitignore` の差を手で読む方式（2 周目までの版）は、相対パス、`.dvc` ファイル、直下のパターン、消して
足し直した行で漏れ、`cache: false` の出力は消す要が無いのに消して戻せないことがある（2 周目と 3 周目の指摘）ので捨てた。main と
ブランチの dvc.lock の outs を比べる案は、`.dvc` ファイルとディレクトリの outs を別に扱う手間と、hook の一覧で足りるかの見込みから
試していない（確かめていない）。

関係する ADR: なし

## 完了の条件

- /topic-start とスキルの入口からマージの直前の段にたどれ、その段の操作が許可の確認なしに実行できる
- 上の残りの 2 点を、直すか、直さない理由を書いて閉じる

## 経緯

- 2026-11-03: #0093 のマージ（dcc9036a2）が、examples/core の持ち主の無くなった古い出力 6 本で pre-merge の hook に止まった。
  /retro の docs/dvc-relock-removed-stage-outs で dvc-relock の節を hook の一覧を使う形に直し、残りをここに書いた
