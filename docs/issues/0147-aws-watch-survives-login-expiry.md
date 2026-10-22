---
status: closed
labels: [aws]
opened: 2026-10-22
closed: 2026-10-22
milestone:
blocked-by: []
resolution: コミット a464c1633
---

# aws_watch.sh が認証切れで止まり、起動し直しが止まった

## 問題

examples/aws/aws_watch.sh は、`aws login` のセッションが切れると終了コード 4 で止まる設計だった（`check_expired`）。
2026-10-21 11:06 に m4-assembly の監視が `reading done failed (rc=1)` で止まった。ログ（~/aws-watch-m4-assembly.log の
722 行目付近）の出力は次のとおりで、`expired` の正規表現（`ExpiredToken|RequestExpired|session has expired`）に合わず、
想定外の失敗（終了コード 3）で終わった。

```text
download failed: s3://.../openmc-2O/done to - An error occurred (400) when calling the HeadObject operation: Bad Request
```

原因は認証切れと見られる（`aws login` は約 12 時間で切れ、S3 の HeadObject は切れた認証に 400 Bad Request を返す）が、
400 がそれだけとは限らず、確かめていない。

監視が止まると、欠けたジョブの起動し直しが止まる。VM と S3 の同期は続くので結果は失われないが、2P、2M、2A-density の
3 本が VM の 24 時間の上限で止まったまま約 20 時間起動し直されなかった（2026-10-22 07:15 に手で監視を再開した）。
止まったことはユーザーが気づかない限り分からない。

関係する ADR: なし

## 完了の条件

- 読み取り（aws_status.sh と done）の失敗で監視が終了せず、`ATTENTION:` の行を出して次の間隔で読み直す
- 読めない間は起動し直さない（欠けたジョブの判定が信用できない）
- 認証切れの疑い（既存の正規表現と、HeadObject の 400 Bad Request）は `aws login` のやり直しを求める文言にする
- 読み直せたら通常に戻る

## 経緯

- 2026-10-22: 起こした
- 2026-10-22: 直した。aws_status.sh と done の読み取りの失敗で終了せず、`ATTENTION:` の行（認証切れの疑いは aws login の
  やり直しを求める文言。HeadObject の 400 Bad Request を含む）を出して次の間隔で読み直す。読めない間は起動し直さない。
  aws_run.sh が認証切れで失敗したときも同じ。起動時の実行の読み込みの認証切れだけは、利用者がその場で見るので 4 で終わる。
  実機での認証切れの再現はしていない（偽の aws_common.sh などによるテスト）
