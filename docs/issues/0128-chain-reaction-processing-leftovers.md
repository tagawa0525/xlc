---
type: Issue
status: closed
labels: [xs, format, chain]
opened: 2026-10-14
closed: 2026-11-19
milestone: M4
blocked-by: []
resolution: "wontfix: レビューで直さなかった指摘を集めた issue は作らない決定（nixfiles の ADR-0066）に合わせ、点ごとに確かめた。温度ごとの読み直しの時間は xs の U-238 で 0.01%、Am-241 で 1.3% 以下で、ほかは 2 か所の重複と形の好み"
---

# チェーンの反応の群の断面積の処理の残り（温度ごとの読み直し、重複した補助関数）

## 問題

ADR-0172（feat/library-reactions）のレビューで、正しさに関わらない次の点が残った。

- xs の `chain_reactions`（crates/xs/src/process.rs）は温度ごとに呼ばれ、MF8、MF9、MF10 をテープから読み直し、温度によらない
  MF10 の表を作り直す。代表の MT の表（U-238 の MT102 などの大きな表）を温度ごとに複製して平均し直す。MF9 の表の外の警告も
  温度ごとに出る。MF8〜10 は `read_evaluation`（`Evaluation`）で 1 度だけ読める。時間への効きは測っていない
- `nil_or_resolved`（crates/format/src/library.rs）が、共通の `nil_or`（crates/format/src/msgpack.rs）と別に残る
- 整数欄の `count` が crates/xs/src/endf/mf8.rs と crates/chain/src/endf.rs にある
- チェーンの `REACTIONS`（crates/chain/src/build.rs）は運動学を `CHAIN_REACTIONS` の位置で対応させる（テスト
  `reaction_kinematics_belong_to_their_mt` が並びの変更を見つける）

関係する ADR: なし

## 完了の条件

上の点を直すか、直さない理由をここに書く。xs の読み直しは、燃焼チェーンの核種のライブラリ（M4 の計画の「進め方」3）の時間で効きを測ってから決める。

## 経緯

- 2026-10-14: ADR-0172 のレビューの 5 周目で、正しさの指摘が尽きたあとに残った点をまとめた
- 2026-11-19: issue の棚卸しで、済んだ指摘を消した（`nil_or` の共通化と、燃焼の反応率が自己遮蔽した捕獲を使うこと。`nil_or_resolved` は残る）
- 2026-11-19: 点ごとに確かめて閉じた（扱いは閉じるコミットの本文）。温度ごとの読み直しは、U-238 と Am-241 の `xs nuclide` の平坦な profile（perf）で、チェーンの反応の平均の積分が 0.01%、表の値の引き（`Tab1::value_within`、全部の使い道の和）が MF9 の分岐のある Am-241 で 1.3% だった（U-238 は MF8〜10 が無い）
