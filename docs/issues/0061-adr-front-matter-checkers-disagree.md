---
status: open
labels: [docs]
opened: 2026-10-11
closed:
milestone:
blocked-by: []
resolution:
---

# ADR のフロントマターの検査が design-records と crates/docs で食い違う（`requires`）

## 問題

`design-records new`（nixfiles の CLI。/adr スキルが使う）はフロントマターに `requires: []` を付け、`design-records check` は
それがないと error にする。一方、このリポジトリの検査（crates/docs の `ADR_KEYS`。`cargo test` が実行する）は
`status`、`date`、`supersedes`、`superseded-by`、`issues` の 5 つだけを受け付け、124 件の ADR のどれにも `requires` は
ない。ADR を書くたびに、`new` が付けたキーを手で消すことになる（ADR-0122、0124、0126 で起きた）。
`design-records check` は、このリポジトリの `.pre-commit-config.yaml` には入っていない。

関係する ADR: なし

## 完了の条件

- `requires` を採るなら、crates/docs の `ADR_KEYS` に足し、既存の ADR のフロントマターに `requires: []` を足して、
  `.pre-commit-config.yaml` で `design-records check --changed` を実行する
- 採らないなら、nixfiles 側で `requires` を任意のキーにする（nixfiles の issue）
- どちらにするかは、`requires` が何のためのキーか（nixfiles の design-records の文書）を確かめて決める

## 経緯

- 2026-10-11: ADR-0126 を `design-records new` で作ったときに気づいた
