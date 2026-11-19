---
type: Issue
status: closed
labels: [verification]
opened: 2026-10-13
closed: 2026-11-19
milestone:
blocked-by: []
resolution: コミット 63358d37d
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
- 2026-11-19: openmc_core_2d.py と openmc_model_3d.py の extract が、作業ディレクトリの model.xml と入力から作ったモデルを、
  ライブラリのパスと statepoint を書く世代を除いて比べ、違えば止まるようにした（openmc_resume.py の `check_model`。テストは
  examples/test_openmc_resume.py と examples/vera/test_openmc_extract.py）。手元の statepoint の 5 本（5A-2D、AWS の m3 の 3A、3B、5B-2D、
  5C-2D）で確かめ、5A-2D だけが、標準のタリーの組（ADR-0152）を足す前のモデルで流したので止まった。入力に `standard-tallies = false` を
  書いて流したモデルに合わせ（verification.md の標準のタリーの組の段落も直した）、`openmc-5A-2D-extract` を再実行して、出力の md5 が
  前と同じことを確かめた
