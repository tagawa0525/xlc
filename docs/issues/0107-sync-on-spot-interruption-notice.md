---
status: open
labels: [aws, reference]
opened: 2026-10-13
closed:
milestone:
blocked-by: []
resolution:
---

# Sync to S3 on the spot interruption notice

## 問題

ADR-0155 の VM は、SIGTERM を受けると finish で最後の同期とログの置き場へのコピーを試みるが、確認（2026-10-13、手元からの終了 4 回）では
1 度も `boot-*.log` が S3 に置かれなかった。原因は確かめていない。候補は、cloud-final.service が KillMode=process で、user-data が
シャットダウンの最後の kill でしか SIGTERM を受けず、そのときにはネットワークが止まっていること（レビューの指摘。確かめていない）。

影響: 中断で失うのは S3 にコピーした最後の statepoint からの分（最大で `statepoint-interval` の約 2 倍と同期の間隔）で、中断のログが残らない。
本番の間隔 100 世代、1 世代約 9 秒（r995。c7a では確かめていない）なら、最大で約 30 分。

関係する ADR: ADR-0155

## 完了の条件

次のどちらかを入れ、中断（AWS FIS のスポットの中断の模擬、または手元からの終了）で、最後の書き終えた statepoint とログが S3 に置かれることを
確かめる。

- 待ちのループで IMDS の `spot/instance-action`（2 分前の通知）を見て、通知が来たら同期とログのコピーを行う
- ジョブを systemd の unit で実行し、ExecStop の同期をネットワークの停止より前に並べる

## 経緯

- 2026-10-13: ADR-0155 の確認とレビューで見つけた。ADR-0155 のブランチでは入れていない（確かめる手段の用意が要るため）
