---
status: open
labels: [vera, dvc]
opened: 2026-11-06
closed:
milestone:
blocked-by: []
resolution:
---

# p5 の格子計算のステージの dvc.lock が今の vera.msgpack に追いついていない

## 問題

`examples/vera/p5/dvc.lock` は、deps の `../vera.msgpack` を 109 行で記録し、そのうち 92 行が md5 0dfa3aca、17 行が 55be71c9（2026-11-06 の main の
93e677a0）。今の main のライブラリ（`examples/vera/vera.msgpack`）は db5c15e8 なので、p5 の格子計算のステージとその下流は changed で、それを上流に持つ
ステージを `dvc repro <yaml>:<ステージ>` で流すと、p5 の格子計算も一緒に走る。

関係する ADR: なし

- 例: `dvc repro --dry examples/depletion/p9-timing/dvc.yaml:timing-p9`（ブランチ perf/0156-warm-start-nem、2026-11-06）は、p5 の 45 ステージと `core-p9`、
  `timing-p9` を出す（`core-p9` は p5 の軸方向の平板と反射体の格子計算の出力を deps に持つ）。p5 の格子計算は Slurm 695（2:14:20）と 696（6:22:23）の規模
- ADR-0230 の帰結は、ライブラリの作り直しで p5 の 109 ステージが changed になり、後回しの規則で再実行しないとする（md5 は 0dfa3aca から db5c15e8）。
  17 行の 55be71c9 がいつのライブラリかは確かめていない
- 再実行は後回しにする（CLAUDE.md の「計算」の「読めなくなった出力の再実行は、それを読む作業が次に来るまで後回し」と「高速化の途中は、高速化の
  測定以外のジョブを流さない」）。#0156 の全炉心の時間の測定（ADR-0231）は、このため p9-timing のステージでなく、p5 の lock どおりのテーブルを使う
  試行で流した（examples/depletion/verification.md の「全炉心の燃焼の時間の内訳（Problem 9 の形）」）

## 完了の条件

p5 の格子計算のステージを今のライブラリで流し直し（それを読む作業が来たとき）、p5 の dvc.lock と、それを上流に持つステージ（p9-timing を含む）を
今の出力で記録する。

## 経緯

- 2026-11-06: #0156 の全炉心の測定を p9-timing のステージで流そうとして、p5 の格子計算が巻き込まれると分かった（ブランチ perf/0156-warm-start-nem）
