---
status: open
labels: [lattice, depletion]
opened: 2026-10-23
closed:
milestone: M4
blocked-by: []
resolution:
---

# MOC の前の解から始めるとき、チェックポイントに前の解を持たせる

## 問題

ADR-0196 のチェックポイントは、燃焼では解き終えた点と次の点の組成を、分岐では解き終えた状態の結果を持ち、再開して完成した
結果が一気に流した結果とビットで一致する。ブランチ perf/moc-warm-start の ADR-0195（2026-10-23 の時点で proposed、main に未マージ）は、
MOC の反復を前の輸送の収束した解（`moc::Start`）から始める。燃焼は直前の輸送の解を、分岐は軸ごとの解（`branches::Starts`）を
持ち越す。チェックポイントがその解を持たないと、再開した最初の輸送は平らな初期値から始まり、結果は収束の判定の程度で
一気に流したものと違う。

ブランチ perf/0137-reuse-shielding-preprocessing（ADR-0194）の Dancoff 係数の持ち越しは、入力が同じ間だけ前の値を使うので、
再開の最初の輸送で同じ入力から解き直す。ビットで同じになるかは確かめていない。

関係する ADR: ADR-0196

## 完了の条件

後からマージする側のブランチで、rebase のあとに `depletion::tests::a_run_resumed_from_a_checkpoint_matches_the_straight_run_bit_for_bit`
と `checkpoint::tests::states_resume_after_the_saved_ones_and_match_the_straight_run`（分岐は main の状態のループを通らないので、
持ち越す解の扱いを足すときにテストも足す）が通ること。持ち越す解をチェックポイントの形式に足す（燃焼は最後の輸送の `Start`、
分岐は軸ごとの `Start`）。

## 経緯

- 2026-10-23: ADR-0196 の実装（feat/lattice-checkpoint-resume）で、並行する perf/moc-warm-start の ADR-0195 を読んで起こした
