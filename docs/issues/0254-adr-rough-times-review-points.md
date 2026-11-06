---
status: open
labels: [docs]
opened: 2026-11-07
closed:
milestone:
blocked-by: []
resolution:
---

# ADR の目安の時間の規則のレビューで直さなかった指摘

## 問題

CLAUDE.md の「計算」の節に、ADR の「検討した案」には `--timed` でないジョブの時間と作業中の版で出した値を目安として書いてよい
とする 1 文を足したブランチ（docs/claude-md-adr-rough-times）の 1 周目のレビューで、/git-merge の「直す指摘の基準」の blocking に
当たらず直さなかった指摘:

1. docs/adr/README.md の「書き方」は、測った数値に「数値を出したコミット（checkout すれば再現できるもの。ADR-0122）と日付」を
   添えるとする。CLAUDE.md の足した文は、作業中の版の値と `commit=<hash>-dirty` を許し、checkout しても再現できない。
   2 つの規則のどちらを優先するかが書かれていない
2. 足した文は 2026-11-07 のユーザーの判断によるが、CLAUDE.md のほかの決定のような「（ユーザー <日付>）」の印と ADR が無い
3. ADR-0234 のほかにも、accepted の ADR の「検討した案」が、パーティション main のジョブの時間を割り込まれうる印なしに引いている
   （指摘の挙げた例: ADR-0204 の Slurm 473、ADR-0228 の Slurm 672、ADR-0201 の Slurm 531 と 532、ADR-0215 の Slurm 625。
   確かめていない）
4. ADR-0234 の「決定と理由」も Slurm 709 と 710 の時間を理由に使うが、足した文は「検討した案」だけに掛かる。#0253 も「決定と理由」に
   触れていない
5. `-c` のジョブは排他でなく CPU を共有する（`TIME:` 行の `shared (-c N)`）が、足した文は割り込まれうることだけを添えさせる
6. 足した文はコミットだけを求め、日付を求めない（docs/adr/README.md の「書き方」は日付も求める）
7. 足した文は「ジョブのコメントの `commit=`」と書くが、終わったジョブの `sacct` の `Comment` の列は空で、`commit=` は SubmitLine の
   `--comment` に残る（`~/.claude/scripts/slurm-run.sh` のコメント。元は nixfiles のリポジトリ）。読む場所を書いていない

関係する ADR: なし

## 完了の条件

各項について、CLAUDE.md か docs/adr/README.md を直すか、直さない理由を書いて閉じる。1 はユーザーの判断（作業中の版の値を
目安として書いてよい）と README の再現の要件のどちらを規則にするかの決定が要る。

## 経緯

- 2026-11-07: docs/claude-md-adr-rough-times の 1 周目のレビューの指摘から起こした
