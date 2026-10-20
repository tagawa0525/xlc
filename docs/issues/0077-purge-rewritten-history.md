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
tagawa0525/xlc-history-backup にある。そのときの書き換えのツール（規則、filter-repo の callback、hash の付け直しを繰り返す
iterate3.sh、全コミットの検証 verify2.py）は references/history-rewrite/ にある（git で追跡しない）。次の書き換え（#0141）は
nixfiles の history-rewrite（nixfiles の ADR-0045）で行い、references/history-rewrite/ はそのときに消す。

関係する ADR: なし

## 完了の条件

ユーザーが、書き換えた後の履歴で困ることがないと判断したら
`git reflog expire --expire=now --all && git gc --prune=now` を行う（取り消せないので、行う前にユーザーに確かめる）。
references/history-rewrite/ は #0141 の書き換えで消す（#0018 で references/ を片付けるときはそのまま）。

## 経緯

- 2026-10-11: 書き換えの後、しばらく様子を見ることにした（メモリに置いていたものを 2026-10-12 に issue にした）
- 2026-10-21: 2 度目の書き換え（#0141。辞書で語を直す）を予定。掃除はその後にまとめて行う
