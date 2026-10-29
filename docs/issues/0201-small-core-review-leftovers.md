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

ブランチ feat/small-core-depletion のレビューの 1 周目（/code-review high）で、誤った数値や落ちるテストにならない次の指摘を残した。

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

関係する ADR: ADR-0191

## 完了の条件

- 各項目を直すか、直さない理由を書いて閉じる

## 経緯

- 2026-10-29: feat/small-core-depletion のレビューの 1 周目で挙がった
