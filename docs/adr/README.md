# ADR（設計判断の記録）

何をなぜ決めたかを、1件1つの決定として記録する。書式は MADR を簡略にしたもの。

- 計画（docs/plans/）には、段階ごとの目的、範囲、作業の順序、検証の対象だけを書き、決定はここへのリンクにする
- 未解決の問題とやることは docs/issues/ に書く。issue で決まったことは ADR にし、issue から参照する

## ファイル

`NNNN-short-name.md`（4桁の連番と、英語の短い名前の kebab-case）。番号は作った順に振り、再利用しない。

```markdown
---
status: accepted      # proposed | accepted | superseded | deprecated
date: 2026-10-08      # 決めた日
superseded-by:        # 置き換えた ADR の番号（superseded のとき）
---

# 題名（決めたことを一文で）

## 背景

何を決める必要があったか。前提と制約。

## 検討した案

- 案 A: 内容。長所と短所
- 案 B: 内容。長所と短所

## 決定と理由

どの案にしたか、なぜか。測った数値や出典があれば書く。

## 結果

この決定で何が起きるか。残る近似、あとで見直す条件、関係する ADR と issue。
```

## 運用

- 採用した ADR は書き換えない（誤字や、リンク先の移動に合わせた直しは除く）。決定を変えるときは新しい ADR を
  作り、古い方を `status: superseded` と `superseded-by` にする
- 測った数値の詳細は examples/*/README.md に置き、ADR からは参照する
- 一覧のファイルは作らない。採用中の ADR は `grep -l '^status: accepted$' docs/adr/[0-9]*.md` で探す
