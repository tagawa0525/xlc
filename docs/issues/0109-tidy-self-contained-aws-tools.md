---
status: open
labels: [aws, tooling]
opened: 2026-10-13
closed:
milestone:
blocked-by: []
resolution:
---

# Tidy the self-contained AWS tools after the ADR-0155 reviews

## 問題

ADR-0155 のブランチのレビュー（13 回）で出た指摘のうち、課金が止まらない、結果が失われる、中断から再開できない、のどれにも当たらない
ものを、まとめて残す。どれも今の動作（examples/vera/aws_*.sh）を誤らせるものではない。

関係する ADR: ADR-0155

- テストの書き方: examples/vera/test_aws_user_data.sh の多くは、生成した user-data の行の有無と順を見ている（How）。振る舞い（What）を
  見ているのは sync_out だけで、finish（done、INTERRUPTED、READY）と待ちのループは行の照合だけ。sync_out と同じく、スタブで実行する形にする
- TDD の履歴: user-data の変更は、テストと実装を同じコミットに入れた（このブランチの依頼で TDD は対象外としたが、flake の check の
  対象になったので、以後は RED と GREEN を分ける）
- 重複: ジョブの名前からスクリプトへの対応（run_remote.sh の case と aws_common.sh の job_script）、VM に送るファイルの一覧
  （remote_send.sh と aws_run.sh）、statepoint を世代の順に並べる処理（openmc_resume.py、user-data、aws_common.sh）。1 か所にまとめる
- ジョブの生死を `pgrep -f` のコマンド行の照合で見ている。run_remote.sh が切り離したジョブの PID を書くようにして、`kill -0` で見る
- aws_status.sh と起動し直しは、ジョブごとに S3 を 3 回呼ぶ。実行ごとに `list-objects-v2 --prefix runs/<名前>/` の 1 回にまとめられる
- 同じ実行の起動し直しを 2 つの端末から同時に実行すると、同じジョブの VM が 2 台立ちうる（排他がない）
- aws_launch.sh も live_instance を使うようになり、終了中（shutting-down）の同名のインスタンスがあると起動を断る（以前は通した）

## 完了の条件

上の項目を、直すか、直さない理由を書いて閉じる。

## 経緯

- 2026-10-13: ADR-0155 のブランチ（feat/aws-autonomous）のレビューから
