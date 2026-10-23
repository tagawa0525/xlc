# ADR（設計判断の記録）

何をなぜ決めたかを、1件1つの決定として記録する。書式は MADR を簡略にしたもの。

- 計画（docs/plans/）には、ステージごとの目的、範囲、作業の順序、検証の対象だけを書き、決定はここへのリンクにする
- 未解決の問題とやることは docs/issues/ に書く。issue で決まったことは ADR にし、issue から参照する
- ADR は `ADR-NNNN`、issue は `#NNNN` と書いて参照する（計画、issue、README、ほかの ADR、コードのコメントのどこでも）

## ファイル

`NNNN-short-name.md`（4桁の連番と、決めたことを動詞で始めた英語の短い名前の kebab-case。例:
`0068-solve-pinwise-slowing-down.md`）。

```markdown
---
status: accepted      # proposed | accepted | rejected | superseded | deprecated
date: 2026-10-08      # 最後に status が変わった日（accepted なら決めた日）
supersedes: []        # この ADR が置き換えた ADR（例: [ADR-0002]）
superseded-by:        # この ADR を置き換えた ADR（例: ADR-0007。複数なら [ADR-0068, ADR-0074]）
amends: []            # この ADR が一部を置き換えた ADR（例: [ADR-0013]）。なければ書かなくてよい
amended-by: []        # この ADR の一部を置き換えた ADR（例: [ADR-0129]）。なければ書かなくてよい
issues: []            # きっかけになった、またはこの ADR で解決した issue（例: ["#0001"]）
requires: []          # 前提にしている ADR（例: [ADR-0068]）。なければ空。design-records check が求める
---

# ADR-NNNN: 題名（決めたことを一文で）

## 背景

何を決める必要があったか。前提と制約。判断の軸（精度、速さ、保守性、評価済みデータだけで決まるか、など）を
2〜4個の箇条書きで添えてもよい。

## 検討した案

- 案 A: 内容。長所と短所
- 案 B: 内容。長所と短所

## 決定と理由

どの案にしたか、なぜか。測った数値や出典があれば書く。

## 確認

決定が成り立つことを何で示しているか。テストの名前、例、比較（examples/vera/verification.md の OpenMC との比較など）と、
あれば判定の閾値。確かめるものがない約束事なら、そう1行で書く。

## 帰結

この決定で何が起きるか。残る近似、あとで見直す条件、関係する ADR と issue。
```

## 書き方

- 測った数値は、その時点の値として、出典のファイルと、数値を出したコミット（checkout すれば再現できるもの。
  ADR-0122）と日付を添える。例: 「VERA 1B で −3080 pcm（examples/vera/verification.md、a764ab0、2026-10-07）」。コミットが簡単に分からなければ、
  ファイルと日付を書く
- 測った数値の詳細は examples/*/verification.md に置き、ADR からは参照する
- 同じ核データ（ENDF/B-VIII.0 を NJOY で処理したもの）の OpenMC との比較は、コード同士の検証（verification）であり、
  実測との比較による妥当性確認（validation）ではない。区別が要るところではそう書く
- 一つの決定の中で比べて採らなかった方式は、その ADR の「検討した案」に書く。決定として本格的に検討して採らなかった
  ものは、`status: rejected` の ADR にしてもよい

## 移行の基準

2026-10-08 に計画（docs/plans/）から ADR へ移したときの基準。この移行にだけ当てはめる。以後の新しい決定は、実装の前に
`status: proposed` の ADR として書く（「運用」）。

- 決定は、main にマージしたコードに実装されたときだけ、独立した ADR にした
- 計画にだけあった決定と、main にマージしていないブランチでだけ実装した決定は、独立した ADR にしない。置き換えた、
  またはいまの決定の ADR の「検討した案」に、「計画で一時採った（YYYY-MM-DD、<コミット>）」または「ブランチ <名前> で
  試した（…）」と書く
- この基準を決める前に番号を振り、基準により独立した ADR にしないもの（ADR-0091、ADR-0092）は、番号を再利用しない
  ために `status: rejected` のファイルとして残し、記録した先を示す

## 運用

- 新しい決定は、実装の前に `status: proposed` の ADR として書く。proposed の ADR は下書きで、main にマージしたあとも、
  accepted になるまで自由に直してよい
- 実装をマージするブランチで、proposed の ADR の本文（「確認」に実際のテストと結果、「帰結」）を実装に合わせて直し、
  同じ変更で `status: accepted` にして `date` をその日にする。採らないと決めたら `status: rejected` にし、理由を「帰結」に
  短く書き、`date` をその日にする
- accepted、rejected、superseded になった ADR は変えられなくなる。直してよいのは、誤字、リンクの直し、`status`、
  `superseded-by`、`amends`、`amended-by` の変更だけである。決定の全体を変えるときは新しい ADR を作り、新しい方の `supersedes` と古い方の
  `superseded-by` で両方向に結び、古い方を `status: superseded` にする
- main にマージする前は、ブランチの上の ADR を status によらず自由に直してよい
- 置き換える ADR は、古い決定のうちまだ有効な部分を書き直して含む（採用中の ADR だけで、いまの状態が分かるようにする）
- 古い決定の一部だけを変える ADR は、古い ADR を `requires` に入れ、「決定と理由」か「帰結」に「これは ADR-NNNN の
  「…」の部分を置き換える。ほかの部分は変わらない」と、置き換える部分を古い ADR の語で書く。ADR 全体の置き換えではない
  ので、`supersedes` / `superseded-by` は使わず、新しい方の `amends` と古い方の `amended-by` で両方向に結び
  （`design-records check` が一致を確かめる）、古い方の status は変えない（例: ADR-0129、ADR-0132、ADR-0183）
- 採用中の ADR は、決定を正しく書いていればよい。実装が決定に追いついていないときは、その差を「帰結」に書き、
  issue に結ぶ（例: ADR-0013 と #0012、ADR-0076 と #0014）
- フロントマターの `issues` には、その ADR のきっかけになった issue と、その ADR で解決する issue だけを書く。issue の側は
  「関係する ADR:」の行に、`issues` にその issue を書いた ADR を挙げ、両方向をそろえる。本文の中で触れるだけの
  issue や ADR は、どちらにも入れなくてよい
- 番号は `design-records new <slug> <title>` が振り、注釈つきのタグ `adr-NNNN`（メッセージは slug）で予約する。
  番号は、全ブランチの履歴、作業ツリー、タグにある最大の番号 + 1 で、再利用しない（ファイルを消した記録の番号も、
  履歴にある限り飛ばす）。このリポジトリにはリモートがないので、タグはローカルにだけ打つ。main に残る飛び番は、
  未マージや放棄のブランチの番号である。予約の仕組みより前に作った未マージの記録は、手で同じ形のタグを打って
  予約する（`git tag -a adr-NNNN -m <slug> main`）
- 番号がぶつかったら（rebase の後など）、`design-records check` と crates/docs の検査（`cargo test`。#0013）が見つける。
  /git-merge が rebase の直後に `design-records -C <worktree> renumber <旧番号>` で振り直す。マージの後は
  `post-merge-cleanup.sh` が `design-records prune` で、main の履歴に入った記録のタグを消す（slug がタグの
  メッセージと一致するものだけ。予約の後に slug を変えた記録のタグは手で消す。マージしなかったブランチのタグは
  残す）。仕組みの決定は nixfiles の ADR-0018
- 一覧のファイルは作らない。採用中の ADR は `grep -l '^status: accepted$' docs/adr/[0-9]*.md` で探す
