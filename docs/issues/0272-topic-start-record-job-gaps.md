---
status: open
labels: [docs]
opened: 2026-11-08
closed:
milestone:
blocked-by: []
resolution:
---

# /topic-start の記録の計算の手順が /git-merge と /delegate の項と合わないところ

## 問題

ブランチ docs/align-review-and-wait-rules（CLAUDE.md と /topic-start を、nixfiles の /git-merge のレビューの周の規則と、/delegate の
5 分より長い待ちの規則を指す形にした）の 2 周目のレビューに出た指摘。どれもこのブランチの前からある /topic-start の
「記録に書く計算の前に main に rebase する」の文で、手順の実行が失敗する誤りではないので non-blocking として直していない:

- 手順 5 は計算の前のレビューを /git-merge の「直す指摘の基準」の歯止めの項だけで指し、「3 周目で blocking を直したら」の項
  （直した head に記録の周を回す）を指さない。3 周目の直しを読まれないまま 1 時間を超える計算に流しうる
- 手順 5 の「流した後のレビューで出た … blocking は直して再実行し」は、流した後の周が記録の周のとき、/git-merge の
  「3 周目で blocking を直したら」（記録の周の指摘は直さず、(a) か (b) は issue に書いてユーザーに判断を求める）と食い違う
- 「計算を投げる直前」の確認は `git log <分岐点>..main` で見るが、/delegate の「作業中の規則」の計算を投げる直前の項は
  `base-moved.sh` で比較元を取ってから `rev-list --count` で見る。同じ確認に 2 つの手順がある（xlc にリモートは無いので、
  今は結果が同じ）
- 「計算を投げる直前」に rebase したとき、3〜5（`scripts/dvc-checkout.py`、チェック、レビュー）をやり直すかを書いていない

## 修正方法

言い換えをやめて /git-merge と /delegate の項を指すか、xlc に固有の部分（rebase をいつも行う、dvc-checkout の順）だけを書く。
どれを残すかは、/delegate の計算を投げる直前の項と照らして決める。

関係する ADR: なし
