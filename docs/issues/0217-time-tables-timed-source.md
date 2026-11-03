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

時間の列を持つ記録の表（どの表かは、スキル verification-record の「計算は Slurm で、時間は `--timed` の排他のジョブのものだけ」の
コマンドで数える）は、比較の JSON の `seconds` から
`examples/verification_tables.py` が描く。`--check` は表と JSON の食い違いを止めるが、`seconds` を書いた上流の格子計算か OpenMC が
割り込まれない timed のジョブで走ったかは見ない。`--timed` を付けずに上流を流すと、CLAUDE.md の「計算」が記録に書かないとする
時間が、表を描き直すだけで記録に入り、チェックは通る。今はスキルの規則（区分 timed か `--timed` で流す）だけで防いでいる。

案（作っていない。決定ではない）: 格子計算の履歴と OpenMC の抽出が、走ったジョブのパーティション（`SLURM_JOB_PARTITION`）と
排他かを出力に書き、比較のステージが `seconds` と一緒に JSON にコピーし、時間を描く関数が timed でない出どころで止まる。履歴の形式を
変えるので、読めなくなる出力と作り直す時間を数えてから決める（CLAUDE.md の「計算」）。

決める前に要ること:

- OpenMC の行: 1C の時間の表（examples/depletion/verification.md）の OpenMC の時間は Slurm 409（`sacct` で Partition は main、
  2026-10-15。`sacct -S 2026-10-01 -r timed` の最も古いジョブは 514 で 2026-10-24）の値で、案のとおりに止めると今の表が落ちる。CLAUDE.md の
  「計算」は数時間の OpenMC の参照解を AWS のスポットに回し、集合体の時間の表（`depletion_assembly_times`）は「OpenMC の時間は
  AWS の機種によるので載せない」（関数の docstring）。1C の表も OpenMC の時間を外すか、timed のパーティションの前の排他のジョブを
  認めるかを決める
- run-cache とチェックポイント: パーティションを出力に書いても、run-cache から戻した出力と、チェックポイントから再開した履歴は、
  前のジョブの時間を持つ（スキルの同じ節）。書く値は、`seconds` を測ったすべての実行のものにする

関係する ADR: なし

## 完了の条件

時間を描く表の出どころが timed のジョブでなければ、`--check` か描く段で止まる。または、止めない理由を決めて閉じる。

## 経緯

- 2026-11-03: 1C の再実行（#0152）で、担当が依頼の外で 1C の格子計算と比較を timed のジョブ（Slurm 658）で再実行した。規則を
  スキル verification-record に足したとき（ブランチ docs/timed-for-time-tables）、機械の確認は作らずにここに残した
