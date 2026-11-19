---
type: Issue
status: open
labels: [docs, tooling]
opened: 2026-11-19
closed:
milestone:
blocked-by: []
resolution:
---

# issue の status の値が OKF の status と意味の違う同じキーにある

## 問題

docs/ は OKF v0.3 のバンドルで（ADR-0256）、OKF の `status` は文書の鮮度の合図（`draft`、`stable`、`deprecated`）である。
OKF は知らない値を `stable` でないものとし、今の知識に絞る consumer はその文書を外す（OKF v0.3 §5.4）。issue の `status` は
進み具合（`open`、`closed`）で、OKF の consumer は開いている issue も閉じた issue も外す。

- ADR の `status` の同じ食い違いは、claude-harness の #261 がキーの改名で扱う（2026-11-19 に open。claude-harness の ADR-0092）
- issue の `status` を読むのは、crates/docs の検査（#0013）、scripts/check-exclusions.py、docs/issues/README.md の
  `grep -l '^status: open$'`

関係する ADR: ADR-0256

## 完了の条件

- [ ] claude-harness の #261 の改名の形が決まった後に、issue の `status` も同じ形で改名するかを決め、決めたとおりに直す

## 経緯

- 2026-11-19: docs/ を OKF のバンドルにしたとき、`status` を変えずに残した
