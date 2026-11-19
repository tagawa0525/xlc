---
type: Issue
status: open
labels: [dvc, tooling]
opened: 2026-11-19
closed:
milestone:
blocked-by: []
resolution:
---

# ADR-0258 の移行の段（正規形の hook、-s の経路、付け直しのコマンド、形式の検査）と ADR-0260 のパスのツールを実装する

## 問題

ADR-0258 と ADR-0260 は決定だけを書いた。ADR-0260 の「移行」のツールと ADR-0258 の段 1〜5 は入った（下の各項）。残りは、パスの
ツールの項にある試行の移し方である。段の中身、段ごとに直す ADR とスキルは ADR-0258 の「移行の順」と ADR-0260 の「移行」に、通すことはそれぞれの
「確認」にある（ここにコピーしない）。着手の順に並べる。段の番号は ADR-0258 のものである。

- パスのツール（ADR-0260）: ファイルを動かしたときに、中身のハッシュで同じファイルかを判定して lock のパスを書き換えるツール。
  ブランチ feat/0298-rewrite-lock-paths で `scripts/rewrite-lock-paths.py` として入れ、ADR-0260 を accepted にした。残りは、既に
  ある試行（`trials/`）の出力を examples に移すときに再計算にしない形で、ADR-0259 の「形を求める時点」（移すときは再実行する）と
  ADR-0260 の「試行を examples に移すとき」（再計算にしない仕組みが要る）が食い違う。ADR-0259 と合わせて決める
- 段 1: 正規形を書き直す pre-commit の hook。ブランチ feat/0298-normalized-form-hook で hook `normalized-forms`
  （`scripts/normalized_forms.py`）として入れた
- 段 2: `scripts/dvc-repro.py` の `-s` の経路（段 3 より先）。ブランチ feat/0298-skip-upstream-route で入れた。手元の生の `-s` は hook で
  止めないと決めた（ADR-0258 の案 H2）
- 段 3: rebase の後の付け直しのコマンド。ブランチ feat/0298-relock-command で `scripts/dvc-relock.py` として入れた
- 段 4: changed の扱いをスキルと計画に書く。ブランチ docs/0298-changed-as-backlog で入れた
- 段 5: 形式を変える修正の検査。ブランチ feat/0298-readable-outputs-check で `scripts/check-readable-outputs.py`（読み手は
  crates/outputs の `read-outputs`。run-checks.sh が `scripts/run-rust-docs-tests.sh` から呼ぶ）として入れ、ADR-0258 を accepted に
  した。Python のツールが読む出力は検査に入れないと決めた（ADR-0258 の検討した案）。それらが古い出力を読めるかは、差分を
  読んで判断する残りである

関係する ADR: ADR-0258、ADR-0260

## 完了の条件

- [x] パスのツールと段 1〜5 が、それぞれのブランチで main に入っている
- [ ] 既にある試行を examples に移すときの扱いが、ADR-0259 と ADR-0260 で決まっている
- [x] ADR-0258 と ADR-0260 の帰結の「実装は決定に追いついていない」の項が今の姿に直っている

## 経緯

- 2026-11-19: ADR-0258 の決定を書き直したときに、実装の残りとして起こした
- 2026-11-20: パスのツールを ADR-0260 に分け、着手の順をパスのツールが最初、形式の検査が最後に直した
- 2026-11-20: パスのツールを入れた（ADR-0260 を accepted）。試行を examples に移すときの扱いを残りに分けた
- 2026-11-20: 段 2（`-s` の経路）を入れた
- 2026-11-20: 段 3（付け直しのコマンド）を入れた
- 2026-11-20: 段 1（正規形の hook）を入れた
- 2026-11-20: 段 4（changed の扱いをスキル、CLAUDE.md、計画、open の issue に書いた）を入れた
- 2026-11-20: 段 5（形式を変える修正の検査）を入れ、ADR-0258 を accepted にした
