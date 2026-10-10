---
status: closed
labels: [docs]
opened: 2026-10-10
closed: 2026-10-10
milestone:
blocked-by: []
resolution: コミット 15ac482
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
- 2026-10-10: マージの途中では MERGE_HEAD の祖先も認めるようにした（15ac482）
- 2026-10-10: 一時的なリポジトリを使うテストは、hook の中では git が渡す GIT_DIR と GIT_INDEX_FILE を受け継ぎ、本物の
  リポジトリに作用する（レビューで指摘）。一時的なリポジトリへの git だけ、GIT_ で始まる環境変数を外した。本物の
  リポジトリへの検査では外さない（hook の中では GIT_INDEX_FILE がマージの結果の index を指す）。GIT_DIR と
  GIT_INDEX_FILE を存在しないパスにして実行し、テストが通ることを手で確かめた（環境変数を変えるテストは、同じ
  プロセスのほかのテストに効くので書いていない）
- 2026-10-10: 試して採らなかった方式: 利用者の hook を `-c core.hooksPath=/dev/null` で止める。コミットの件名の hook が
  止まらなかったので、コミットとマージに `--no-verify` を付け、件名も Conventional Commits にした
