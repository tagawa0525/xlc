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
  （#0283）。このパスを書いた追跡している入力は、#0283 の切り分けの 6 つと examples/depletion/issue-0284/lattice-2A.toml（注釈だけ）のほかに
  examples/depletion/ の下に 40（branches/ 22、history/ 12、micro/ 2 など。`git ls-files examples | grep '\.toml$' | xargs grep -l p2/lattice-2A.msgpack`、
  注釈だけのものを含む）あり、今の既定でこの履歴から作るテーブルはどれも同じ違いを持つ。#0279 の段 5 のテーブル（examples/depletion/issue-0279/core-micro/）は、
  `lattice@2A` の入力のコピーを今の既定で解いた examples/depletion/issue-0284/lattice-2A.msgpack を読む（下の「経緯」の 2026-11-13）。作り直すと、下流の DVC のステージは同じ yaml の `compare@2A`（OpenMC との比較）。
  上流の `xs-nuclide` と `xs` も、`bin/xs` が deps で changed（`xlc-format` のテーブルの形式の変更で `bin/xs` が変わる。#0117）
- 候補（どれも確かめていない）: 履歴に自己遮蔽の設定（`ufg-free-gas` など、k を変える自己遮蔽の入力）を書き、テーブルの入力と違えば
  止める。比べる項目の範囲（自己遮蔽のほかの設定、MOC の離散化）と、意図して違う設定で解く入力（#0283 の caps-400.toml の
  ような切り分け）の書き方は決めていない
- 上の候補の形式の変更で読めなくなるもの（実装していない）: 設定を履歴の msgpack（`xlc-depletion-history`。crates/lattice/src/depletion.rs の
  `HistoryFile`）の必須の項目にすると、今の履歴はどれも読めなくなる（後方互換を作らない。CLAUDE.md の「方針」）。履歴を書く DVC のステージは
  13（examples/depletion/p2 の `lattice` 8、examples/depletion/dvc.yaml の `lattice-2A`、`lattice-1C@rings1`、`lattice-1C@rings4`、
  `lattice-mox1` の 4、examples/depletion/p9-product の `lattice-2A` 1。`dvc.lock` の `bin/lattice` のステージのうち、deps の入力に
  `[depletion]` のあるもの）で、ほかに `.dvc` で追う試行の履歴がある。履歴を読む DVC のステージは examples/depletion/p9-product の
  `branch-2A`（ログの `total_seconds` 14,958 秒、Slurm 752）と、examples/depletion/p9-2A の `branch-2A`（../p9-product/lattice-2A.msgpack を読む。
  `dvc.lock` が無く、実行していない）。Python の比べるツール（compare_depletion.py、compare_depletion_assembly.py、compare_core_micro.py など）は
  履歴の `format` と使う項目だけを見るので、項目を足すだけなら読める（`format` の名前を変えれば止まる）。
  履歴を書く 13 本の時間の目安: p2 の 8 本は前の作り直し（Slurm 664、666〜672、400 kT の既定、2026-11-03〜04）の経過時間の和が 10 時間 41 分
  （比較のステージを含む。668〜672 は割り込まれうるジョブで、記録に使える時間でない）。examples/depletion/dvc.yaml の 4 本はログの `total_seconds`
  の和 2,530 秒、p9-product の `lattice-2A` はログの `total_seconds` 4,818 秒（Slurm 747。割り込まれうるジョブで目安）。和は約 12.7 時間で、
  今の既定の自己遮蔽の 400 kT からの比（ADR-0246 の帰結の、2A の燃焼の最初の 1 段で 1.23 倍。燃えた燃料の比は測っていない）を
  自己遮蔽の分に掛けると、もう少し長い（全部に掛けた上限で約 15.6 時間。見込みで、確かめていない）。この和は、上流のライブラリの
  作り直しを含まない。`bin/xs` が changed なので、p2 の 545 核種と examples/depletion の 65 核種の `xs-nuclide` と、それぞれの `xs` も走る
  （`xs-nuclide` の CPU 時間の記録の和は p2 29,209 秒、examples/depletion 11,616 秒。各 yaml のディレクトリの stage-cpu-seconds.json。runner が
  同時に流すので経過時間はこれより短く、`xs` の多群の作成の時間は含まない）。履歴を読む側（上の 2 つの `branch-2A` と、`[history]` を持つ
  追跡している TOML の 101）は、作り直した履歴を読むまで動かない。作り直しは読む作業が来たときに行う（CLAUDE.md の「計算」）

関係する ADR: なし

## 完了の条件

- 履歴とテーブルの設定の食い違いが、テーブルを作るときに分かる（止めるか、入力で認めたときだけ通す）
- `lattice@2A` を今の既定で作り直すか、作り直さずに読む入力に違いを書くかが決まっている（作り直しは、それを読む作業が来たときに行う。
  CLAUDE.md の「計算」）

## 経緯

- 2026-11-13: #0283 で、#0279 の段 5 の参照とのずれが、履歴とテーブルの自由ガスの上限の既定の違いと分かり、同じ違いを止める仕組みが
  無いので起こした
- 2026-11-13: #0279 の段 5 のやり直しで、`lattice@2A` を作り直さず、その入力のコピー（examples/depletion/issue-0284/lattice-2A.toml。違いは
  ライブラリのパスの書き方と出力の置き場だけ）を今の既定で DVC のステージの外で解き、段 5 のテーブルの履歴にした（xlc 62e0786d3 の
  `cargo build --release`、crates は main aef2e9dca と同じ。Slurm 823、main のパーティションで途中で止められたので時間は記録しない。md5 0f5c5d1f、
  `dvc add`）。`lattice@2A` を作り直すと、上流の p2 のライブラリ（545 核種の `xs-nuclide` と `xs`）も `bin/xs` の変更で changed で走り、
  前の担当の数えで合計約 4.3 時間（確かめていない。`scripts/count-changed-stages.py examples/depletion/p2/dvc.yaml` は 5 分で終わらず
  止めた）。ライブラリは中身が同じバイトになる見込み（確かめていない）で、作り直しはこの比較の結論を変えないので、`lattice@2A` と
  下流の `compare@2A`（OpenMC との比較の記録）の作り直しは、それを読む作業が来るまで後回しにした（CLAUDE.md の「計算」）
