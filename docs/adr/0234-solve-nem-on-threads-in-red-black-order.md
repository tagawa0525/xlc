---
status: proposed
date: 2026-11-06
requires: [ADR-0067, ADR-0137, ADR-0146, ADR-0214, ADR-0232]
supersedes: []
superseded-by:
issues: ["#0156"]
---

# ADR-0234: 炉心の固有値計算を、ノードを市松の 2 色に分けた順で解き、同じ色のノードを複数のスレッドで解く

## 背景

全炉心の燃焼（Problem 9 の形）は、CMFD（[ADR-0232](0232-accelerate-nem-with-node-cmfd.md)）の後も全体の 95.9 % が固有値計算で、1 回の解の固有値計算は
1 つのスレッドで解く（examples/depletion/verification.md の「固有値計算の CMFD（小さな炉心と全炉心）」。Slurm 704、77a720cd6、2026-11-06。全体 374.2 s、
固有値計算 45 回 358.9 s）。全炉心の燃焼の速度の対策は、Xe とホウ素の探索の同時の収束、CMFD、複数のスレッドの順に進める（ユーザー 2026-11-06 の判断の要約。
ADR-0232 の背景）。これはその 3 番目。

NEM のスイープ（`nem::Core::sweep_group`）と粗メッシュの固有値問題（`xlc_cmfd::Coarse::solve` のブロック Gauss-Seidel）は、ノード（セル）を直方体の並びの順に
1 つずつ解き、先に解いた隣の値をすぐ使う。この順のままでは、ノードを並列に解けない。

判断の軸: 正確さ > 実行速度 > 保守（CLAUDE.md の「方針」）。

- 収束の判定（`nem::Options` の許容）は緩めない。解く順を変えた解は、同じ判定を満たす解として許容の 2 倍の中で一致すること
- 結果はスレッドの数によらないこと（[ADR-0137](0137-compute-options-as-input.md) は、スレッドの数を結果を変えない前提で環境変数から読むことを許す）
- 前の解き方（直方体の並びの順に 1 つずつ）を入力で選べる形で残す（ADR-0231 の `initial-guess = "flat"`、ADR-0219 の `xenon-boron = "nested"`、
  ADR-0232 の `acceleration = "none"` と同じ扱い）

## 検討した案

（測定の後に書く）

## 決定と理由

（測定の後に書く）

### 帰結

### 確認
