---
status: open
labels: [core, depletion]
opened: 2026-10-29
closed:
milestone: M4
blocked-by: []
resolution:
---

# サイクルのステップ（ADR-0218）のレビューの残り

## 問題

ADR-0218 のブランチ（feat/core-cycle-steps）のレビューの 1 周目で、誤った数値や落ちるテストにならない次の指摘を残した。

- `crates/core/src/depletion.rs` の `deplete`: 使わない状態の並びを `None` で持つ `Vec<Option<Layout>>` と、「state … is not laid out」の
  エラーは、入力（`DepletionInput::new`）がどのステップも使わない状態を止めるので main からは届かない。`MicroNodes::new` の `files_of` の
  作り直しと `table {file}` の代わりの名前、並びの `files` の比べ合わせも、`with_banks` が同じファイルの並びを作るので届かない。
  ステップの使う状態の並びだけを `Vec<Layout>` で持てば、分岐が減る
- 同じ: バンクの位置が同じでホウ素や温度だけが違う状態も、状態ごとに `with_banks` で並び全体とミクロ燃焼のノードの使い方を作る。
  Problem 9 の 30 状態でバンク D の位置は数種で、バンクの位置で並びを共有すれば作る回数とメモリが減る（大きさは測っていない。#0156）
- `examples/depletion/compare_core_cycle.py`: 参照の点を EFPD だけで組にし、状態を見ない。状態を変える点（ADR-0218 の案 P2 で次の
  ステップの状態で解く）の参照との Δρ は棒の価値を含み、`--report` の Δρ の範囲の行に混ざる（記録の本文は点ごとに読んでいる）
- 同じ: `rho_pcm` が compare_core_micro.py、compare_core_partial_rods.py に続く 3 つ目の複製（#0197 の比較のツールの複製と同じ）。
  `core_history.py` にまとめられる
- 入力: 崩壊だけのステップに書いた `state` は、後の `state` を省いたステップに引き継がれる（ADR-0218 の「省けば前のステップと同じ」のとおり）。
  停止の間だけの状態（HZP など）を書くと、次の出力のステップもその状態になる。崩壊だけのステップの `state` を止めるか、引き継がないかは、
  Problem 9 の入力を作るときに決める

関係する ADR: ADR-0218

## 完了の条件

- 各項目を直すか、直さない理由を書いて閉じる

## 経緯

- 2026-10-29: feat/core-cycle-steps のレビューの 1 周目（/code-review high）で挙がった
