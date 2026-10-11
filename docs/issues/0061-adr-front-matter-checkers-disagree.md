---
status: closed
labels: [docs]
opened: 2026-10-11
closed: 2026-10-11
milestone:
blocked-by: []
resolution: コミット 7adeadec
---

# ADR のフロントマターの検査が design-records と crates/docs で食い違う（`requires`）

## 問題

`design-records new`（nixfiles の CLI。/adr スキルが使う）はフロントマターに `requires: []` を付け、`design-records check` は
それがないと error にする。一方、このリポジトリの検査（crates/docs の `ADR_KEYS`。`cargo test` が実行する）は
`status`、`date`、`supersedes`、`superseded-by`、`issues` の 5 つだけを受け付け、main の ADR のどれにも `requires` は
ない。ADR を書くたびに、`new` が付けたキーを手で消すことになる（ADR-0122、0124、0126 で起きた）。
`design-records check` は、このリポジトリのチェック（run-checks.sh）では走らない（`.pre-commit-config.yaml` もない）。

関係する ADR: なし

## 完了の条件

- `requires` を採るなら、crates/docs の `ADR_KEYS` に足し、既存の ADR のフロントマターに `requires: []` を足して、
  チェック（run-checks.sh）で `design-records check` を実行する
- 採らないなら、nixfiles 側で `requires` を任意のキーにする（nixfiles の issue）
- どちらにするかは、`requires` が何のためのキーか（nixfiles の design-records の文書）を確かめて決める

## 結末

`requires` を採った。crates/docs の `ADR_KEYS` に足し、項目が `ADR-NNNN` の形で存在する ADR であることを検査する。
既存の ADR すべてに `requires: []` を足した。あわせて、design-records が推奨する節名に合わせ、ADR の `## 結果` を
`## 帰結` に改名した（本文は変えていない）。`design-records check` の error は status が rejected の ADR 0091 と 0092
の節の不足だけで、nixfiles 側で必須の節を免除する変更を待つ。`design-records check` を run-checks.sh で実行することは、
このリポジトリの外（`~/.claude`）の変更なので、ここでは入れていない。

## 経緯

- 2026-10-11: ADR-0126 を `design-records new` で作ったときに気づいた
- 2026-10-11: `requires` を採ると決めた（ユーザーの判断）。ADR のフロントマターと検査、節名を直して閉じた
