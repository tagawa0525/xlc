---
status: closed
labels: [dvc, tooling]
opened: 2026-10-12
closed: 2026-10-12
milestone:
blocked-by: []
resolution: コミット 329efb44
---

# dvc.lock の schema 行を検査する

## 問題

3fe1a45a（chore: relock the staged xs stages on top of main）で作り直した `examples/vera/dvc.lock` から、先頭の
`schema: '2.0'` の行が落ちた。dvc はこの行の無い lock を読めず、`dvc status`、`dvc repro`、`dvc commit`、`dvc stage list`
が `'./dvc.lock' validation failed. required key not provided, in schema` で止まった。c51c6c74 で 1 行足して直したが、
追跡している lock を検査する仕組みがないので、同じ欠落は次の relock でも通る。

落ちた操作について確かめたこと:

- dvc は schema の無い lock を読めず、書きもしない。dvc で作った lock から `schema` の行を `sed` で消す実験（scratch の
  2 段階の dvc.yaml）では、`dvc status`、`dvc commit -f`、`dvc repro -s` がすべて同じ検証エラーで止まり、lock は消えた
  ままだった。dvc が lock を作るときは、この行を必ず出す。
- したがって欠落は、dvc の外での書き込みで起きた。候補は、main の lock と実行した lock の段階を手で合わせた書き換え。
  3fe1a45a の作業の記録は残っていない。ただし 3fe1a45a の diff は、先頭の行だけが消え、段階の書式は dvc のものである（手で段階を
  組み合わせた形と合う）。どの操作で組み合わせたかは確かめていない。`dvc-relock` の手順には、lock を手で書く箇所はない。

検査に使うツール:

- `dvc status` は出力のハッシュを計算し、vera では約 27 秒かかる（examples/vera の worktree。m0 と reflector は 1.5 秒以内）。
- `dvc stage list <dvc.yaml>` は lock を読み込んで検証し、ハッシュを計算しない。vera で約 2 秒。壊れた lock では終了コード 255。

関係する ADR: なし

## 完了の条件

- 追跡している `examples/*/dvc.lock` を含むコミットで、dvc が読めない lock は pre-commit の hook（`dvc-lock`、
  `scripts/check-dvc-lock.sh`）で止まる。schema の行を消した vera の lock で止まることを確かめる。
- `dvc-relock` の手順の最後に、作り直した lock の先頭が `schema: '2.0'` で、`dvc stage list` が終了コード 0 で終わることの確認がある。
- run-checks.sh では走らない（それは ~/.claude のツール）。rebase は pre-commit を実行しないので、rebase で lock が壊れたときは
  `dvc-relock` の確認で見つける。

## 検査の方式

`dvc stage list <dvc.yaml>` を使う。lock を読み込んで検証するが、出力のハッシュを計算しないので、hook に置ける時間で終わる
（vera で約 2 秒。`dvc status` は約 27 秒）。

## 確かめた結果

完了の条件（schema のない lock で hook が止まる）は、コミット 329efb44 の上で確かめた: `prek run dvc-lock --files examples/vera/dvc.lock` は、
先頭の schema を消した lock で終了コード 1、戻すと 0。

- `scripts/check-dvc-lock.sh` と pre-commit の hook（`prek run dvc-lock --files ...`）は、現在の m0、reflector、vera の lock を
  通す（終了コード 0）。
- 先頭の `schema: '2.0'` を消した vera の lock は、script も hook も終了コード 1 で止まる（hook の出力は
  `dvc stage list` の検証エラーと FAILED の行）。元に戻すと通る。
- 手順は、`sed -i '/^schema/d' examples/vera/dvc.lock` のあと `prek run dvc-lock --files examples/vera/dvc.lock`、
  それから `git checkout -- examples/vera/dvc.lock`。

テストは書いていない。リポジトリに shell の検査のテストの仕組みがなく（crates/docs のテストは Markdown の決まりを見る）、
作るには dvc の実行を含むテストが要るため。上の手順を、lock を変えるたびに手で行う。

## 残る限界

- hook は commit のときだけ効く。rebase と merge では走らないので、rebase の解決で壊れた lock は止まらない。
  dvc-relock の手順 6 が、relock のあとに同じ確認をする。
- hook が呼ぶ `dvc` は PATH のものである。direnv（devShell）の中では devShell の dvc（wrapper）になるが、devShell の外では
  `python3.14-dvc-3.67.1` の素の dvc が使われる。m0、reflector、vera の 3 つで、`nix develop -c dvc stage list` と素の
  `dvc stage list` の出力は一致した（出力の md5 を比べた）。ほかの段階の出力は確かめていない。
- prek が PATH に無いと、グローバルの pre-commit は `.pre-commit-config.yaml` の検査を飛ばして通す（警告だけ）。この環境では
  prek が入っているので効くが、PATH の違う shell や GUI からの commit では hook が効かない。
- 欠落の原因（dvc の外の書き込み）は確かめていない（上の「問題」）。hook は原因を直さず、壊れた lock を commit させないだけ。

## 経緯

- 2026-10-12: 3fe1a45a で欠落し、c51c6c74 で 1 行を足して直した（マージ f21bcc14）。ADR-0123（段階の駆動）、ADR-0138（xs の
  段階の分割）、#0075（核種ごとの段階）の作業で起きた。
