# ADR（設計判断の記録）

何をなぜ決めたかを、1件1つの決定として記録する。書式は MADR を簡略にしたもの。

- 計画（docs/plans/）には、段階ごとの目的、範囲、作業の順序、検証の対象だけを書き、決定はここへのリンクにする
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
issues: []            # きっかけになった、またはこの ADR で解決した issue（例: ["#0001"]）
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

決定が成り立つことを何で示しているか。テストの名前、例、比較（examples/vera/README.md の OpenMC との比較など）と、
あれば判定の閾値。確かめるものがない約束事なら、そう1行で書く。

## 結果

この決定で何が起きるか。残る近似、あとで見直す条件、関係する ADR と issue。
```

## 書き方

- 測った数値は、その時点の値として、出典のファイルと、数値を記録したコミットと日付を添える。例:
  「VERA 1B で −3080 pcm（examples/vera/README.md、fcd48f1、2026-10-07）」。コミットが簡単に分からなければ、
  ファイルと日付を書く
- 測った数値の詳細は examples/*/README.md に置き、ADR からは参照する
- 同じ核データ（ENDF/B-VIII.0 を NJOY で処理したもの）の OpenMC との比較は、コード同士の検証（verification）であり、
  実測との比較による妥当性確認（validation）ではない。区別が要るところではそう書く
- 一つの決定の中で比べて採らなかった方式は、その ADR の「検討した案」に書く。決定として本格的に検討して採らなかった
  ものは、`status: rejected` の ADR にしてもよい

## 運用

- 採用した ADR は書き換えない。直してよいのは、誤字、リンクの直し、置き換えられたときの `status` と
  `superseded-by` だけである。決定を変えるときは新しい ADR を作り、新しい方の `supersedes` と古い方の
  `superseded-by` で両方向に結ぶ
- 番号は、ブランチの上で既存の最大の番号 + 1 を振り、使い回さない。マージの前に main に rebase し、番号が
  ぶつかっていたら振り直す（重複を見つける検査は #0013 で作る。それまでは手で確かめる）
- 一覧のファイルは作らない。採用中の ADR は `grep -l '^status: accepted$' docs/adr/[0-9]*.md` で探す
