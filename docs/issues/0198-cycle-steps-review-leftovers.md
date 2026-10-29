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
- 同じ: 崩壊だけのステップに書いた `state` は、崩壊そのものに効かず、案 P2 で前の出力のステップの終わりの点を解く状態になる（停止の前の
  点が HZP の k とホウ素で記録される）。上の項目と一緒に決める
- 標準出力の `point` の行（ADR-0191 の形のまま）は日数と状態を持たないので、崩壊だけのステップの点は前の点と同じ EFPD と燃焼度の行になり、
  順でしか見分けられない（履歴のファイルは `days` と `state` を持つ）
- `deplete` の崩壊だけのステップは、チェーンだけで決まる崩壊の行列をステップごとに作り直す。`point_state` の `map_or(0, …)` と `used` の
  `|| point_state(0) == k` は、ステップが空でないことを確かめた後では効かない
- `compare_core_cycle.py` の `densities` は、核分裂のないノード（数密度が空）を持つ履歴で `IndexError` になる（集合体 1 体の確認では起きない）

関係する ADR: ADR-0218

## 完了の条件

- 各項目を直すか、直さない理由を書いて閉じる

## 経緯

- 2026-10-29: feat/core-cycle-steps のレビューの 1、2 周目（/code-review high）で挙がった
