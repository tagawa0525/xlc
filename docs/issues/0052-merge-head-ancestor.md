---
status: closed
labels: [docs]
opened: 2026-10-10
closed: 2026-10-10
milestone:
blocked-by: []
resolution: コミット 2aff807
---

# マージのときの検査が、閉じた issue のブランチの上のコミットを祖先と認めない

## 問題

閉じた issue の `resolution: コミット <hash>` は、HEAD の祖先のコミットを指すことを crates/docs の検査で確かめる。
docs/issues/README.md はブランチの上のコミットを書いてよいとしているが、マージのコミットを作る前の検査（hook の
`run-checks.sh --merge`）では HEAD がまだ main なので、ブランチの上のコミットが祖先にならず、マージが止まる
（fix/0051-grazing-segments のマージで、2026-10-10）。

関係する ADR: なし

## 完了の条件

- マージの途中（MERGE_HEAD がある）では、マージするブランチの祖先のコミットも認める

## 経緯

- 2026-10-10: fix/0051-grazing-segments のマージで見つけた
- 2026-10-10: マージの途中では MERGE_HEAD の祖先も認めるようにした（2aff807）
