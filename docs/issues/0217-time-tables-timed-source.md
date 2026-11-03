---
status: open
labels: [tooling, verification]
opened: 2026-11-03
closed:
milestone:
blocked-by: []
resolution:
---

# 時間の表の出どころが timed のジョブかを機械で確かめない

## 問題

時間の列を持つ記録の表（examples/depletion/verification.md の depletion-1c-times、depletion-assembly-times。数え方はスキル
verification-record の「計算は Slurm で、時間は `--timed` の排他のジョブのものだけ」）は、比較の JSON の `seconds` から
`examples/verification_tables.py` が描く。`--check` は表と JSON の食い違いを止めるが、`seconds` を書いた上流の格子計算か OpenMC が
割り込まれない timed のジョブで走ったかは見ない。`--timed` を付けずに上流を流すと、CLAUDE.md の「計算」が記録に書かないとする
時間が、表を描き直すだけで記録に入り、チェックは通る。今はスキルの規則（区分 timed か `--timed` で流す）だけで防いでいる。

案（作っていない。決定ではない）: 格子計算の履歴と OpenMC の抽出が、走ったジョブのパーティション（`SLURM_JOB_PARTITION`）と
排他かを出力に書き、比較のステージが `seconds` と一緒に JSON にコピーし、時間を描く関数が timed でない出どころで止まる。履歴の形式を
変えるので、読めなくなる出力と作り直す時間を数えてから決める（CLAUDE.md の「計算」）。

## 完了の条件

時間を描く表の出どころが timed のジョブでなければ、`--check` か描く段で止まる。または、止めない理由を決めて閉じる。

## 経緯

- 2026-11-03: 1C の再実行（#0152）で、担当が依頼の外で 1C の格子計算と比較を timed のジョブ（Slurm 658）で再実行した。規則を
  スキル verification-record に足したとき（ブランチ docs/timed-for-time-tables）、機械の確認は作らずにここに残した
