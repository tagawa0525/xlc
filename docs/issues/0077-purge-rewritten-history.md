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
tagawa0525/xlc-history-backup にある。2026-10-21 に 2 度目の書き換え（#0141。nixfiles の history-rewrite で、ja-terms の
辞書で語を直した）をした。その作業ディレクトリ ~/.local/share/history-rewrite/xlc-2026-10-21/ に、書き換える前のミラー
（original.git、889 MB）、元の番号から新しい番号への対応表（commit-map）、worktree の控えのパッチがある。

関係する ADR: なし

## 完了の条件

ユーザーが、書き換えた後の履歴で困ることがないと判断したら
`git reflog expire --expire=now --all && git gc --prune=now` を行い、上の作業ディレクトリを消す（どちらも取り消せないので、
行う前にユーザーに確かめる）。

## 経緯

- 2026-10-11: 書き換えの後、しばらく様子を見ることにした（メモリに置いていたものを 2026-10-12 に issue にした）
- 2026-10-21: 2 度目の書き換え（#0141）をした。掃除は 2 回分をまとめて行う
