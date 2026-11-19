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

# OKF の索引の ADR と issue の行に説明が無い

## 問題

docs/ は OKF v0.3 のバンドルで（ADR-0256）、索引（`index.md`）の行の説明は、文書のフロントマターの `description` から
`design-records index` が取る。OKF v0.3 の §8 は、索引の行に `description` を含めることを推奨（SHOULD）にしている。

- ADR と issue はフロントマターに `description` を持たないので、docs/adr/index.md と docs/issues/index.md の行は題だけになる。
  `design-records check` は説明の無い行を違反にも注意にもしない（2026-11-19 に xlc の docs/ で確かめた）
- ADR に `description` を足すと、`design-records check` が「unexpected front matter key」の違反にする（claude-harness の
  crates/design-records/src/lib.rs の `OPTIONAL_KEYS` は `type`、`amends`、`amended-by` だけ。lock 3d1a4f85）。
  crates/docs の検査も ADR と issue の知らないキーを違反にする
- 件数は ADR 249 件、issue 292 件（2026-11-19 の main）

関係する ADR: ADR-0256

## 完了の条件

- [ ] 説明を持たせるか、題だけの索引でよいかを決める
- [ ] 持たせるなら、design-records と crates/docs が ADR と issue の `description` を受け付け、既存の文書に説明を付ける
  （題や「背景」「問題」の最初の文から機械的に取るか、1 件ずつ書くか）

## 経緯

- 2026-11-19: docs/ を OKF のバンドルにしたとき、説明を足さずに題だけの索引にした（数百件を書く前に決めるため）
