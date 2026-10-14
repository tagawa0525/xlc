---
status: open
labels: [xs, format, chain]
opened: 2026-10-14
closed:
milestone: M4
blocked-by: []
resolution:
---

# チェーンの反応の群の断面積の処理の残り（温度ごとの読み直し、重複した補助関数）

## 問題

ADR-0172（feat/library-reactions）のレビューで、正しさに関わらない次の点が残った。

- xs の `chain_reactions`（crates/xs/src/process.rs）は温度ごとに呼ばれ、MF8、MF9、MF10 をテープから読み直し、温度によらない
  MF10 の表を作り直す。代表の MT の表（U-238 の MT102 などの大きな表）を温度ごとに複製して平均し直す。MF9 の表の外の警告も
  温度ごとに出る。MF8〜10 は `read_evaluation`（`Evaluation`）で 1 度だけ読める。時間への効きは測っていない
- `nil_or` と同じ補助関数が crates/format/src/library.rs、table.rs、crates/depletion/src/chain.rs にあり、`nil_or_resolved` も残る
- 整数欄の `count` が crates/xs/src/endf/mf8.rs と crates/chain/src/endf.rs にある
- チェーンの `REACTIONS`（crates/chain/src/build.rs）は運動学を `CHAIN_REACTIONS` の位置で対応させる（テスト
  `reaction_kinematics_belong_to_their_mt` が並びの変更を見つける）
- lattice の `shield_cell` はリングの捕獲を UFG の値にするが、`shielded_nuclide` は捕獲を捨てる（燃焼の反応率で使うときに運ぶ）

関係する ADR: なし

## 完了の条件

上の点を直すか、直さない理由をここに書く。xs の読み直しは、燃焼チェーンの核種のライブラリの時間（ADR-0169）で効きを測ってから決める。

## 経緯

- 2026-10-14: ADR-0172 のレビューの 5 周目で、正しさの指摘が尽きたあとに残った点をまとめた
