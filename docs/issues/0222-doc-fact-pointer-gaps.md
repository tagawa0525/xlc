---
status: closed
labels: [docs]
opened: 2026-11-04
closed: 2026-11-08
milestone: M4
blocked-by: []
resolution: コミット 4bf526b84
---

# CLAUDE.md と /topic-start の (c) の指し方がスキルの原文より狭い・広い

## 問題

ブランチ docs/claude-md-doc-fact-blocking（CLAUDE.md「作業の流れ」の「計画と issue だけのブランチ」の項と、/topic-start の
レビューの項を、nixfiles の /git-merge の blocking の基準 (c)（文書の本文が引いた先の一次情報と食い違う誤りと、規則が求める出どころが
無いこと。nixfiles の ADR-0066 の「## 変更」2026-11-04）と「3 周目で blocking を直したら」の例外の項を指す形にした）の記録の周
（3 周の後。head 085375ccb）に出た指摘。どれもスキルの原文を言い換えた文が原文より狭いか広いもので、手順の実行が失敗する誤りでは
ないので non-blocking として直していない:

- /topic-start の「出典で再現する」は、(c) の再現の定義（食い違いなら出典、出どころが無いなら規則が求める項目と無い箇所と見つけた
  引いた先）より狭い
- /topic-start と CLAUDE.md の「記録の周の後に残ったとき」は、原文の「記録の周で blocking が出たら」（3 周目の直しが新しく作った
  誤りを含む）より狭く読める
- /topic-start の担当が報告に書くものが「blocking」だけで、原文の「blocking と引いた先」の引いた先が無い
- /topic-start は担当の側（直さず止まる）だけで、親が自分で書いたブランチの側（親が引いた先を開いて確かめてから直し、例外の周を回す）が
  無い
- CLAUDE.md の (c) の参照が「計画と issue だけのブランチ」の項の中にあり、(c) がそのブランチだけの規則に読める（原文はブランチの種類を
  限らない）
- CLAUDE.md の「例外の項に従う（ユーザーに判断を求める条件を含む）」は、(a) か (b) が 1 件でもあるときの扱いが例外の項でなく本項に
  あるのと食い違う

## 修正方法

言い換えをやめて、/git-merge の (c) と例外の項を指すだけにする（決定は nixfiles の 0066 の ADR。nixfiles の #576 と同じ型の問題で、原文を 2 か所にコピーすると細部がずれる）。
CLAUDE.md の参照は「計画と issue だけのブランチ」の項の外（レビューの一般の規則）に置く。

関係する ADR: なし

## 経緯

- 2026-11-04: docs/claude-md-doc-fact-blocking の記録の周の直さない指摘から起こした
- 2026-11-08: docs/align-review-and-wait-rules で CLAUDE.md と /topic-start を /git-merge の項を指す形にし、(c) の参照を
  「計画と issue だけのブランチ」の項の外に移して閉じた。残った言い換えは #0272
