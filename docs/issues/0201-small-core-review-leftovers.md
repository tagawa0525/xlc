---
status: open
labels: [core, depletion]
opened: 2026-10-29
closed:
milestone: M4
blocked-by: []
resolution:
---

# 小さな炉心の燃焼のブランチのレビューの残り

## 問題

ブランチ feat/small-core-depletion のレビューの 1〜3 周目と記録の周（/code-review high）で、誤った数値や落ちるテストにならない次の指摘を残した。

- `examples/depletion/compare_core_steps.py` は、中点の固有値計算の回数を core のログ（tracing の `depletion midpoint` の行の `solves=`）から読む。
  ログの形（イベントの名前、項目、ANSI の色）が変わると読めず、ログを残さない再実行と比べられない。ADR-0191 の案 O2 は中点の値を履歴に
  書かないと決めたので、回数を履歴に持たせるなら ADR-0191 の出力の決定を直す
- `crates/core/src/depletion.rs` の履歴の読み込み（`FileIn`）は、`midpoint_boron` を文字列で読み、`"search"`、`"start"` を手で列挙子に対応させる
  （rmpv の `from_value` が文字列を単位の列挙子として読まないため）。書く側の serde の名前と 2 か所になる。`IntoDeserializer` で列挙子の
  `Deserialize` を通せば 1 か所になる
- `examples/depletion/branches/core-small-*.toml` の 8 つのケースは、使い捨てのスクリプトで基準のサイクルを分けて書き、生成のツールと `--check` を
  残していない。基準のサイクルを変えると 8 つを書き直す（p2_inputs.py の形が先例）
- 記録の core の 8 つの燃焼（分割 16 は約 25 分、ほぼ 1 コア）は、依頼の「数分なら手元、数十分なら Slurm」と、リポジトリの CLAUDE.md の
  「xlc の集合体以上は Slurm」に反して手元で流した。結果は決定的で、再実行しても数値は変わらないので再実行していない。#0156 の時間は
  8 つを同時に流した目安
- `compare_core_steps.py` の `key` は、分割の間の点を日数を 1e-6 日に丸めて組にする。分けた時間の和の丸め（1e-13 ほど）が丸めの境を
  またぐと、正しい組が「no points at」で止まる（誤った数値にはならない）。点の数を確かめるようにしたので、点の番号で組にして日数の差を
  許容で確かめる形にできる（2 周目）
- `verification_tables.py` の `core_steps_ratio` は、丸めて 0 になる負の比を「−0.00」と描く（`signed` は符号を付けない）（3 周目）
- `compare_core_steps.py` の `midpoint` は、中点のホウ素が `search` と `start` のちょうど 2 つのときだけ作り、ほかの名前は黙って `None` にする。
  `core_steps_midpoint` は `steps` に `search` が無い要約で `TableError` でなく `KeyError` になる（3 周目）
- `compare_core_cycle.py` の、中点を初めの点のホウ素で解いた履歴を止める検査（`default_midpoints`）は、参照では `single_state` と 2 度呼ばれ、
  サイクルの履歴の側は test_compare_core_cycle.py にテストがない（3 周目）
- crates/core/src/depletion/tests.rs の `midpoints_at_the_starting_boron_fall_to_first_order` の入力で選んだ案 S1 の確かめは、最後の点の集合体 A1 の
  燃焼度だけを比べ、点のホウ素や状態を変えた後の中点のホウ素を比べない（3 周目）
- `core_history.py` の `default_midpoints` は `"start"` だけを止め、`midpoint_boron` の無い履歴と知らない値を既定の探索として通す（Rust の読み込みと
  compare_core_steps.py は止める）（記録の周）
- 文の誤り: crates/core/src/input.rs のモジュールの説明の `start` は「ADR-0191 の案 S1 との差を測るため」でなく案 S1 そのもの（採用の案 S2 との差を
  測るため）。`CoreHistory::midpoint_boron` と `MidpointBoron` の説明の「中点を解いたホウ素濃度」は、濃度でなく決め方（記録の周）
- `compare_core_steps.py` の `differences` は start − search にも使うが、片方だけ出力のあるノードのエラーが「one split only」と分割を指す（記録の周）
- crates/core/src/depletion.rs の中点のホウ素（`given`）を、ステップの中点のクロージャの中で `depletion.midpoint_boron` から決める（ケースの設定で、
  `solve` の隣で 1 度決められる）（記録の周）

関係する ADR: ADR-0191

## 完了の条件

- 各項目を直すか、直さない理由を書いて閉じる

## 経緯

- 2026-10-29: feat/small-core-depletion のレビューの 1〜3 周目と記録の周で挙がった
