---
status: open
labels: [docs]
opened: 2026-11-08
closed:
milestone:
blocked-by: []
resolution:
---

# /topic-start と CLAUDE.md のレビューと記録の計算の文が /git-merge と /delegate の項と合わないところ

## 問題

手順の実行が失敗する誤りではないので、/git-merge の「直す指摘の基準」で non-blocking として残した点:

- /topic-start の「記録に書く計算の前に main に rebase する」の手順 5 は、計算の前のレビューを /git-merge の「直す指摘の基準」の
  歯止めの項だけで指し、「3 周目で blocking を直したら」の項（直した head に記録の周を回す）を指さない。3 周目の直しを読まれないまま
  1 時間を超える計算に流しうる
- 同じ手順 5 の「流した後のレビューで出た … blocking は直して再実行し」は、流した後の周が記録の周のとき、/git-merge の
  「3 周目で blocking を直したら」（記録の周の指摘は直さず、(a) か (b) は issue に書いてユーザーに判断を求める）と食い違う。
  「それ以外は issue に残すか、まとめて 1 度だけ再実行する」は、non-blocking を 1 件の issue にまとめて別のトピックで扱う
  /git-merge の non-blocking の項より広い
- 同じ節の「計算を投げる直前」の確認は `git log <分岐点>..main` で見るが、/delegate の「作業中の規則」の計算を投げる直前の項は
  `base-moved.sh` で比較元を取ってから `rev-list --count` で見る。同じ確認に 2 つの手順がある
- 同じ節は、「計算を投げる直前」に rebase したとき、3〜5（`scripts/dvc-checkout.py`、チェック、レビュー）をやり直すかを書いていない
- 同じ節の「2 の rebase は /delegate の同じ節の計算を投げる直前の項より広い」は、/delegate のその項が rebase を「この場合だけ行う」と
  限ることを言い換えで弱めている
- /topic-start の「マージ」の節の括弧は、文書だけのブランチを「文書だけのブランチ」の項に、ほかを「3 周目で blocking を直したら」の項に
  分けるが、/git-merge の「文書だけのブランチ」の項は記録の周と (a) の扱いで「3 周目で blocking を直したら」の項を使う
- CLAUDE.md の「作業の流れ」の「計画と issue だけのブランチ」は、/git-merge の「文書だけのブランチ」の項（説明の文書も入る）より
  狭い名前で対象を挙げる
- CLAUDE.md の「計画と issue だけのブランチは、事実の誤りと矛盾を直したらマージする」を /git-merge の項を指す形に替えたので、
  引いた先の一次情報の無い文書の中の矛盾（計画の 2 つの順序の食い違いなど）は (c) に入らず non-blocking になった。xlc でそれを
  blocking に残すかを決める
- #0222 の `resolution` は最後のコミットだけを指し、#0222 の点の多くを直した前のコミットを指さない

関係する ADR: なし

## 完了の条件

上の点ごとに、/git-merge と /delegate の項を指す形に直すか、xlc に固有の部分（rebase をいつも行う、dvc-checkout の順）として
残す理由を書いて閉じる。CLAUDE.md の矛盾の扱いはユーザーが決める。

## 経緯

- 2026-11-08: docs/align-review-and-wait-rules（CLAUDE.md と /topic-start を /git-merge と /delegate の項を指す形にした）のローカルの
  レビューの 2、3 周目の直さない指摘から起こした
