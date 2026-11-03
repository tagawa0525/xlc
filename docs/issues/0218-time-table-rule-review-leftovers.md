---
status: open
labels: [tooling, verification]
opened: 2026-11-03
closed:
milestone:
blocked-by: []
resolution:
---

# 時間の表の規則のレビューで残した指摘

## 問題

スキル verification-record の「計算は Slurm で、時間は `--timed` の排他のジョブのものだけ」に時間の表の規則を足したブランチ
（docs/timed-for-time-tables、#0217）の /code-review の 3 周で、手順のとおりに実行して失敗する誤りでないものとして直さなかった指摘:

- 「上流（xs など）は先に揃えておく」が揃え方を書いていない。p2 の格子計算の上流の xs-nuclide は、CLAUDE.md の「計算」のとおり
  `examples/stage_runner.py` で流すと書くとよい
- #0217 が引く 1C の時間の表の OpenMC の値（Slurm 409、2026-10-15）に、数値を出したコミットを添えていない（README の「運用」）
- 同じ節で、チェックポイントの `seconds` が前後の実行の和である説明が、新しい項とチェックポイントの項の 2 か所にある
- 時間を描く表の数え方は、初めに `grep -n -B6 '<関数名>,$' examples/verification_tables.py` で `TABLES` の項を引く形にしたが、
  `depletion_assembly_times` では前の内包表記の `depletion_assembly_free_gas_k` と 2A、2O の JSON まで出し、JSON の一覧が 6 行を
  超える項では表の名前が切れるので、`TABLES` の関数の本文を `inspect.getsource` で読む 1 行の Python に替えた（1 周目の指摘）。
  その 1 行も、補助の関数を通して `seconds` を読む表は拾わない（スキルに書いた）

3 周目の直し（59ac456ea）に回した記録の周が挙げた blocking（/git-merge の「直す指摘の基準」の歯止めで直さず、ユーザーの判断に回す）:

- スキルの「格子計算のステージの出力を消してから `dvc repro --no-run-cache <yaml>:<ステージ>`」は、格子計算を実行しない。DVC 3.67.1
  の `run_stage`（dvc/stage/run.py）は force でなければ `stage_cache.restore` を呼び、`restore`（dvc/stage/cache.py 198〜203 行）は
  `changed_stage()` でなく deps が揃い dvc.lock に outs の hash があれば、run-cache の有無（`--no-run-cache`）より前に dvc.lock から
  ステージを作り、キャッシュの前の出力を戻す（確かめた: インストールされた DVC のソース）。`seconds` は前のジョブのままになる。
  Slurm 658 が走ったのは lock を戻したため（#0152）。代わりの候補は `dvc repro -f --downstream <yaml>:<格子計算のステージ>`
  （force が目標とその下流だけに掛かるかは確かめていない）か、#0152 のように lock を戻す形
- 続けて流す比較のステージは、上流の OpenMC のステージが changed なら、数時間の参照解を timed のジョブの中で実行しうる。
  この 2 段にも CLAUDE.md の「計算」の数え方（`scripts/count-changed-stages.py`）を求めていない
- 「`--no-run-cache` は上流の changed のステージも実行する」の下で上流に changed の xs-nuclide があると、CLAUDE.md の「計算」が
  `examples/stage_runner.py` で流すとする xs-nuclide を `dvc repro` で実行する（上の 1 項目目の揃え方の欠けと同じ根）

関係する ADR: なし

## 完了の条件

上の各項を直すか、直さない理由を書いて閉じる。

## 経緯

- 2026-11-03: ブランチ docs/timed-for-time-tables のレビューの 3 周で残した
