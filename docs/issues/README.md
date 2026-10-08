# issues

未解決の問題とやることを、1件1ファイルで管理する。何をなぜ決めたかは ADR（docs/adr/）に書き、ここには
まだ決まっていないこと、分かっていないこと、これからすることを書く。決まったことは ADR にし、issue からは ADR を
`ADR-NNNN` と書いて参照する。

## ファイル

`NNNN-short-name.md`（4桁の連番と、英語の短い名前の kebab-case）。番号は作った順に振り、使い回さない。

```markdown
---
status: open          # open | closed
labels: [lattice, shielding]
opened: 2026-10-08
closed:
milestone:            # 任意。M3、M4 など
blocked-by: []        # 先に閉じる必要のある issue（例: ["#0001"]）
resolution:           # 閉じるとき: ADR-NNNN | マージのコミット | wontfix: 理由
---

# 題名

## 問題
何が起きているか、数値、再現の手順。

関係する ADR: ADR-NNNN（その ADR の前付けの `issues` にこの issue があるもの）

## 完了の条件
何が分かれば、何ができれば閉じるか。

## 経緯
- 2026-10-08: 起こした経緯や途中で分かったこと
```

## 運用

- 閉じるときは `status: closed` と `closed` の日付と `resolution` を書き、ファイルは消さずに残す。閉じる変更は、対応したブランチの
  中でコミットする
- ブランチ名とコミットの本文で番号を示す（例: ブランチ `fix/0004-frame-gap`、本文に `Refs: #0004`。決定に結びつくなら `Refs: ADR-0002` も）
- 番号は、ブランチの上で既存の最大の番号 + 1 を振る。マージの前に main に rebase し、番号がぶつかっていたら
  振り直す（重複を見つける検査は #0013 で作る。それまでは手で確かめる）
- 一覧のファイルは作らない（手で書く一覧は各ファイルとずれる）。開いている issue は
  `grep -l '^status: open$' docs/issues/[0-9]*.md` で探す
