---
status: open
labels: []
opened: 2026-10-12
closed:
milestone:
blocked-by: []
resolution:
---

# 履歴の書き換えで残った古いオブジェクトを消す

## 問題

2026-10-11 に履歴を書き換えた（openmc-*.msgpack を履歴から消し、語を直し、記録の中の hash を付け直した）。書き換える前の
オブジェクトは reflog から参照されていて、リポジトリに残っている。元の履歴は GitHub の private リポジトリ
tagawa0525/xlc-history-backup にある。そのときのツール（規則、filter-repo の callback、hash の付け直しを繰り返す iterate3.sh、
全コミットの検証 verify2.py。git で追跡しない references/history-rewrite/）は、#0141 で nixfiles の history-rewrite に
置き換えて消した。2026-10-21 に 2 度目の書き換え（#0141。nixfiles の history-rewrite で、ja-terms の
辞書で語を直した）をした。その作業ディレクトリ ~/.local/share/history-rewrite/xlc-2026-10-21/ には、書き換える前のミラー
（original.git、889 MB）、元の番号から新しい番号への対応表（commit-map）、worktree の未コミットの変更の控えのパッチがあったが、
2026-11-07 に消した（経緯）。書き換える前のオブジェクトの一部は、このリポジトリにまだ残る（旧ブランチ feat/core-cmfd の書き換える前の
先端 13fddc36e は、2026-11-07 に `git cat-file -t` で commit と出た。ほかにどれだけ残るかは確かめていない）。

関係する ADR: なし

## 完了の条件

ユーザーが、書き換えた後の履歴で困ることがないと判断したら
`git reflog expire --expire=now --all && git gc --prune=now` を行う（取り消せないので、行う前にユーザーに確かめる）。
2 度目の書き換えの作業ディレクトリは、2026-11-07 に消した（経緯）。

## 経緯

- 2026-10-11: 書き換えの後、しばらく様子を見ることにした（メモリに置いていたものを 2026-10-12 に issue にした）
- 2026-10-21: 2 度目の書き換え（#0141）をした。掃除は 2 回分をまとめて行う
- 2026-11-06: feat/core-cmfd の未コミットの変更（控えのパッチ `worktree-2.patch` と `worktree-2.after.patch` の lib.rs の 1 行 `pub mod cmfd;` と、パッチでなく別に控えた `cmfd.rs.before`
  の `crates/numeric/src/cmfd.rs`）は、そのままはコミットしていない。
  ソルバの切り出しの試しは ADR-0232 の案 S1 に書き、新しいクレート `xlc-cmfd` に取り込んだ（6bf685f7c）。ブランチと worktree は
  ユーザーの判断で片付けた（#0076）
- 2026-11-07: 控えのコピーを残すか消すか（#0249）の問いに「2消していいのでは？」（ユーザー 2026-11-07）を受け、親が作業ディレクトリ
  ~/.local/share/history-rewrite/xlc-2026-10-21/ を消した。消す前に確かめることにしていた控えのパッチの変更は、上の 2026-11-06 の項のとおり、
  そのままはコミットせず、ソルバの切り出しの試しを ADR-0232 の案 S1 に書いて `xlc-cmfd` に取り込んだ。完了の条件から作業ディレクトリの項を外した
