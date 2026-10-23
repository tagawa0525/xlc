---
status: open
labels: [verification]
opened: 2026-10-13
closed:
milestone:
blocked-by: []
resolution:
---

# 手で流した OpenMC の statepoint のモデルと、extract の入力のモデルの一致を確かめない

## 問題

5A-2D の OpenMC の実行はステージでなく、ステージ `openmc-5A-2D-extract`（examples/vera/dvc.yaml）が手で流した
`references/openmc/vera-5A-2D/` の statepoint から結果を抜く。extract は statepoint の粒子数、世代数、非アクティブバッチ数、乱数の
種を入力（`p5/openmc-5A-2D.toml`）と突き合わせて、違えば止まる（feat/0086-periphery-constants）。しかし `problem`、
`cross-sections` と、モデルを作るコード（`openmc_core_2d.py`、`openmc_model.py`。どちらも extract の deps）の変更は確かめない。
それらを変えて再実行しないと、extract は古いモデルの statepoint から結果を抜き、ステージは最新として扱われる。

作業ディレクトリの `model.xml` と、入力から作り直したモデルを比べる案は、`model.xml` が実行した worktree の絶対パス
（ライブラリ）を含むので（examples/openmc_run.py のコメント）、そのままでは比べられない。確かめていない。

関係する ADR: なし

## 完了の条件

extract が、statepoint を出したモデルと入力のモデルの一致を確かめる（パスを除いた model.xml の比較など）か、OpenMC の
実行をステージにして（外部の計算機で流す手順。ADR-0151）この問題をなくす。

## 経緯

- 2026-10-13: feat/0086-periphery-constants のレビューで挙がった
