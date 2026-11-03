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

関係する ADR: なし

## 完了の条件

上の各項を直すか、直さない理由を書いて閉じる。

## 経緯

- 2026-11-03: ブランチ docs/timed-for-time-tables のレビューの 3 周で残した
