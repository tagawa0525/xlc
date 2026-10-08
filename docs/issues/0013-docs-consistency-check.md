---
status: closed
labels: [docs]
opened: 2026-10-08
closed: 2026-10-08
milestone: M3
blocked-by: []
resolution: crates/docs（cargo test）
---

# 文書の整合性の検査（番号の重複、状態、置き換えのリンク、参照先）

## 問題

docs/adr/README.md と docs/issues/README.md は、番号の重複を検査で見つけることを前提にしているが、その検査はまだない。
ADR と issue の前付け（`status` の値、`supersedes` と `superseded-by` が両方向にそろっていること）や、`ADR-NNNN` と
`#NNNN` の参照先が存在すること、issue と ADR の参照が両方向にあることも、いまは手で確かめている。

関係する ADR: なし（ADR と issue の運用は docs/adr/README.md と docs/issues/README.md）

## 完了の条件

次を確かめる検査を作り、pre-commit かチェックのスクリプトで走らせる。

- ADR と issue の番号の重複
- `status` の値が決められたものか
- `supersedes` と `superseded-by` が両方向にそろっているか
- `ADR-NNNN` と `#NNNN` の参照先が存在するか

## 経緯

- 2026-10-08: ADR を入れたときの READMEs が、まだない検査を前提にしていた
- 2026-10-08: crates/docs に検査を作った。番号の重複、前付けのキーと `status`、置き換えの両方向のリンク、閉じた issue の日付と
  解決、ADR の `issues` と「関係する ADR:」の対応、`ADR-NNNN` と `#NNNN` の参照先、相対リンクを確かめる。ワークスペースの
  `cargo test` で走るので、チェックのスクリプト（run-checks.sh）が文書だけの変更でも走らせる
