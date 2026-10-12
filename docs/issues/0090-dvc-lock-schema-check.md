---
status: open
labels: [dvc, tooling]
opened: 2026-10-12
closed:
milestone:
blocked-by: []
resolution:
---

# dvc.lock の schema 行を検査する

## 問題

eed4106a（chore: relock the staged xs stages on top of main）で作り直した `examples/vera/dvc.lock` から、先頭の
`schema: '2.0'` の行が落ちた。dvc はこの行の無い lock を読めず、`dvc status`、`dvc repro`、`dvc commit`、`dvc stage list`
が `'./dvc.lock' validation failed. required key not provided, in schema` で止まった。ca94f387 で 1 行足して直したが、
追跡している lock を検査する仕組みがないので、同じ欠落は次の relock でも通る。

落ちた操作について確かめたこと:

- dvc は schema の無い lock を読めず、書きもしない。dvc で作った lock から `schema` の行を `sed` で消す実験（scratch の
  2 段階の dvc.yaml）では、`dvc status`、`dvc commit -f`、`dvc repro -s` がすべて同じ検証エラーで止まり、lock は消えた
  ままだった。dvc が lock を作るときは、この行を必ず出す。
- したがって欠落は、dvc の外での書き込みで起きた。候補は、main の lock と走らせた lock の段階を手で合わせた書き換え。
  eed4106a の作業の記録は残っていないので、どの操作だったかは確かめていない。`dvc-relock` の手順（dvc の出力を
  コミットする操作）には、lock を手で書く箇所はない。

検査に使う道具:

- `dvc status` は出力のハッシュを計算し、vera では約 27 秒かかる（examples/vera の worktree。m0 と reflector は 1.5 秒以内）。
- `dvc stage list <dvc.yaml>` は lock を読み込んで検証し、ハッシュを計算しない。vera で約 2 秒。壊れた lock では終了コード 255。

関係する ADR: なし

## 完了の条件

- 追跡している `examples/*/dvc.lock` を含むコミットで、dvc が読めない lock は pre-commit の hook（`dvc-lock`、
  `scripts/check-dvc-lock.sh`）で止まる。schema の行を消した vera の lock で止まることを確かめる。
- `dvc-relock` の手順の最後に、作り直した lock の先頭が `schema: '2.0'` で、`dvc stage list` が終了コード 0 で終わることの確認がある。
- run-checks.sh では走らない（それは ~/.claude の道具）。rebase は pre-commit を走らせないので、rebase で lock が壊れたときは
  `dvc-relock` の確認で見つける。

## 経緯

- 2026-10-12: eed4106a で欠落し、ca94f387 で 1 行を足して直した（マージ 8dfc1af3）。ADR-0123（段階の駆動）、ADR-0138（xs の
  段階の分割）、#0075（核種ごとの段階）の作業で起きた。
