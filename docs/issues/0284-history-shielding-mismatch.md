---
status: open
labels: [lattice, depletion]
opened: 2026-11-13
closed:
milestone:
blocked-by: []
resolution:
---

# 格子計算の燃焼の履歴と、それから作るテーブルの自己遮蔽の設定の食い違いを止める仕組みが無い

## 問題

- テーブルの `[history]`（crates/lattice/src/history.rs）は、履歴のライブラリの温度と輸送する核種が入力と合うかは確かめるが、履歴を
  解いた自己遮蔽の設定（自由ガスの上限など）は履歴の msgpack に無く、確かめられない。既定が変わった後に古い履歴から作ったテーブルは、
  テーブルの点の k∞ が履歴から黙ってずれる（大きさは #0283 の「切り分け」）
- examples/depletion/p2/lattice-2A.msgpack（examples/depletion/p2/dvc.yaml の `lattice@2A`）は ADR-0243 と ADR-0246 の前の既定で解かれたまま
  （#0283）。このパスを書いた追跡している入力は、#0283 の切り分けの 6 つのほかに examples/depletion/ の下に 44（branches/ 22、history/ 12、
  issue-0279/core-micro/ 4、micro/ 2 など。`git ls-files examples | grep '\.toml$' | xargs grep -l p2/lattice-2A.msgpack`、注釈だけのものを含む）あり、
  今の既定でこの履歴から作るテーブルはどれも同じ違いを持つ。作り直すと、下流の DVC のステージは同じ yaml の `compare@2A`（OpenMC との比較）
- 候補（どれも確かめていない）: 履歴に自己遮蔽の設定（`ufg-free-gas` など、k を変える自己遮蔽の入力）を書き、テーブルの入力と違えば
  止める。比べる項目の範囲（自己遮蔽のほかの設定、MOC の離散化）と、意図して違う設定で解く入力（#0283 の caps-400.toml の
  ような切り分け）の書き方は決めていない
- 上の候補の形式の変更で読めなくなるもの（実装していない）: 設定を履歴の msgpack（`xlc-depletion-history`。crates/lattice/src/depletion.rs の
  `HistoryFile`）の必須の項目にすると、今の履歴はどれも読めなくなる（後方互換を作らない。CLAUDE.md の「方針」）。DVC のステージの履歴は
  13（examples/depletion/p2 の `lattice` 8、examples/depletion/dvc.yaml の `lattice-2A`、`lattice-1C@rings1`、`lattice-1C@rings4`、
  `lattice-mox1` の 4、examples/depletion/p9-product の `lattice-2A` 1。`dvc.lock` の `bin/lattice` のステージのうち、deps の入力に
  `[depletion]` のあるもの）で、ほかに `.dvc` で追う試行の履歴がある。p2 の 8 本は前の作り直し（Slurm 664、666〜672、400 kT の既定、
  2026-11-03〜04）の経過時間の和が 10 時間 41 分で、自己遮蔽が 400 kT の約 1.2 倍（ADR-0246 の帰結）として約 13 時間、ほかの 2A の 2 本を
  足して約 14〜15 時間の見込み（1C と mox1 のピンセルは短い。確かめていない）。履歴を読む入力（`[history]` を持つ追跡している TOML の
  101）は、作り直した履歴を読むまで動かない。作り直しは読む作業が来たときに行う（CLAUDE.md の「計算」）

関係する ADR: なし

## 完了の条件

- 履歴とテーブルの設定の食い違いが、テーブルを作るときに分かる（止めるか、入力で認めたときだけ通す）
- `lattice@2A` を今の既定で作り直すか、作り直さずに読む入力に違いを書くかが決まっている（作り直しは、それを読む作業が来たときに行う。
  CLAUDE.md の「計算」）

## 経緯

- 2026-11-13: #0283 で、#0279 の段 5 の参照とのずれが、履歴とテーブルの自由ガスの上限の既定の違いと分かり、同じ違いを止める仕組みが
  無いので起こした
