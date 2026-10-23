---
status: accepted
date: 2026-10-14
requires: [ADR-0155]
supersedes: []
superseded-by:
issues: ["#0166"]
---

# ADR-0168: Sync the newest statepoint once it has been written, and let a relaunch shorten the statepoint interval

## 背景

ADR-0155 の同期は、走っている間は最新の statepoint（書いている途中でありうる）を S3 にコピーしないので、中断で残るのは 1 つ前の
statepoint までである。2026-10-13〜14 の本番 m3（c7a.4xlarge のスポット 4 台）では中断が 30〜60 分に 1 回起き、statepoint の間隔が
100 世代の 5B と 5C-2D（1 つに約 45 分）では、中断の間隔より 2 つ分の時間が長く、4 時間で進みが残らなかった（5C-2D は 0300 のまま）。

- 中断のたびに失う時間を、statepoint の間隔の約 2 倍から約 1 倍に縮める
- 参照解の物理の結果（k、tally）は statepoint の間隔に依らない

## 検討した案

- **A. 最新の statepoint も、書き終えていればコピーする**: ファイルの mtime が今から 120 秒より古ければ書き終えているとみなす。
  OpenMC が statepoint を書くのは数秒で、書いている間は mtime が更新される。短所は、120 秒という値が書く時間の見積もりに依ること
  （書く時間がそれを超える大きさの statepoint では、書いている途中のものを写しうる。写っても、再開は壊れたものを消して
  1 つ前から続ける）
- **B. 起動し直しのときに `statepoint-interval` を上書きする**: 中断が続く間は間隔を縮めて、失う分を減らす。S3 の `src/` の入力は
  書き換えない。起動し直しの記録（`launches/`）に使った間隔を残す。起動し直しごとに明示し、前の起動の値は引き継がない（間隔は
  結果に影響しないが、黙って変わる設定を作らない）
- 採らなかった案（ADR-0155 の「検討した案」で試して捨てた、書かれて 60 秒たたない最新の statepoint を同期から外す方式）は、
  完了の直後の同期でも最後の世代が外れた。A は走っている間だけに適用し、ジョブが成功した（`done` が 0）ときは従来どおりすべてコピーする
- 書き終わりを OpenMC の出力のログ（`Creating state point`）で判断する案: ログの形式への依存が増えるので採らない

## 決定と理由

A と B を採る。A は同期の関数 `sync_out` の変更だけで、`--delete` の扱いは変えない（走っている間はしない）。B は
`aws_run.sh --statepoint-interval <世代>` から `run_remote.sh`、`openmc_resume.py` の `run --statepoint-interval` に渡す。
中断の多い時期は A だけでも失う分は約 1 つ分（間隔 100 で約 45 分）に縮み、B で 50 にすれば約 22 分になる。

## 確認

- flake の check `aws-user-data`（examples/vera/test_aws_user_data.sh）: `sync_out` が、mtime が 120 秒より古い最新の statepoint を
  含め、120 秒以内のものは除くこと。`run_remote.sh` への間隔の引数が user-data に入ること
- `openmc_resume.py` の `--statepoint-interval` は、参照解の Python で `set_statepoints` が入力の値より指定の値で statepoint の世代を決めることを
  手元で確かめた（Python のテストはこのリポジトリにない。実機の実行は確かめていない）
- `aws_run.sh --dry-run --job openmc-3A --statepoint-interval 50 m3`（実機の起動はしていない）

## 帰結

- 中断で失うのは、最後にコピーした statepoint からの分（最大で間隔の約 1 つ分と同期の間隔）になる。ADR-0155 の「間隔の約 2 倍」を置き換える
- 参照解の statepoint の間隔は物理の結果に影響しないが、記録（verification.md）には起動し直しで使った間隔を書く
- 起動し直しの VM は、手元の `aws_run.sh` が作る新しい user-data で起動する（`launches/` の `TOOLS_COMMIT`）。走っている VM は
  古い user-data のまま
