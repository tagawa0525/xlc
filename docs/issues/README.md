# issues

未解決の問題とやることを、1件1ファイルで管理する。何をなぜ決めたかは ADR（docs/adr/）に書き、ここには
まだ決まっていないこと、分かっていないこと、これからすることを書く。決まったことは ADR にし、issue からは ADR を
`ADR-NNNN` と書いて参照する。

## ファイル

`NNNN-short-name.md`（4桁の連番と、英語の短い名前の kebab-case）。番号は `design-records new` が作った順に振り、再利用しない。

```markdown
---
status: open          # open | closed
labels: [lattice, shielding]
opened: 2026-10-08
closed:
milestone:            # 任意。M3、M4 など
blocked-by: []        # 先に閉じる必要のある issue（例: ["#0001"]）
resolution:           # 閉じるとき: ADR-NNNN | コミット <hash> | wontfix: 理由
---

# 題名

## 問題
何が起きているか、数値、再現の手順。

関係する ADR: ADR-NNNN（その ADR のフロントマターの `issues` にこの issue があるもの）

## 完了の条件
何が分かれば、何ができれば閉じるか。

## 経緯
- 2026-10-08: 起こした経緯や途中で分かったこと
```

## 運用

- 検証の記録（examples/*/verification.md）の数値を引くときは、記録のファイルと、数値を出したコミット（ADR-0122）と日付を書く。
  記録は今の姿に書き直されるので（ADR-0117）、節の題だけでは引き先が消えうる
- 閉じるときは `status: closed` と `closed` の日付と `resolution` を書き、ファイルは消さずに残す。閉じる変更は、対応したブランチの
  中でコミットする
- `resolution` は、`ADR-NNNN`（決めた ADR。その ADR のフロントマターの `issues` にこの issue を挙げる）、`コミット <hash>`
  （解決したコミット。7〜40 桁の英小文字の16進）、`wontfix: 理由` のどれかにする。マージのコミットは、閉じる変更を
  コミットする時点では分からないので、ブランチの上のコミットを書いてよい。マージの前にブランチを rebase したら
  hash を書き直す（検査が、HEAD の祖先でない hash を見つける）
- ブランチ名とコミットの本文で番号を示す（例: ブランチ `fix/0004-frame-gap`、本文に `Refs: #0004`。決定に結びつくなら `Refs: ADR-0002` も）
- issue は `design-records new --kind issue <slug> <title>` で作る。番号を振ってタグ `issue-NNNN` で予約し、
  表題の行だけを書くので、フロントマターと節は手で埋める。番号の振り方、衝突の修正方法（`renumber --kind issue`）、タグの
  解放（`prune`）は ADR と同じ（docs/adr/README.md）。フロントマターの検査は crates/docs が行う（#0013）
- 一覧のファイルは作らない（手で書く一覧は各ファイルとずれる）。開いている issue は
  `grep -l '^status: open$' docs/issues/[0-9]*.md` で探す
