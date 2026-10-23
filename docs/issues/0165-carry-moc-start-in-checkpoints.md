---
status: closed
labels: [lattice, depletion]
opened: 2026-10-23
closed: 2026-10-23
milestone: M4
blocked-by: []
resolution: ADR-0196
---

# MOC の前の解から始めるとき、チェックポイントに前の解を持たせる

## 問題

ADR-0196 のチェックポイントは、燃焼では解き終えた点と次の点の組成を、分岐では解き終えた状態の結果を持ち、再開して完成した
結果が一気に流した結果とビットで一致する。ADR-0195 は、
MOC の反復を前の輸送の収束した解（`moc::Start`）から始める。燃焼は直前の輸送の解を、分岐は軸ごとの解（`branches::Starts`）を
持ち越す。チェックポイントがその解を持たないと、再開した最初の輸送は平らな初期値から始まり、結果は収束の判定の程度で
一気に流したものと違う。

ADR-0194 の Dancoff 係数の持ち越しは、入力が同じ間だけ前の値を使うので、再開の最初の輸送で同じ入力から解き直す。

関係する ADR: ADR-0194、ADR-0195、ADR-0196

## 完了の条件

後からマージする側のブランチで、rebase のあとに `depletion::tests::a_run_resumed_from_a_checkpoint_matches_the_straight_run_bit_for_bit`
と `checkpoint::tests::states_resume_after_the_saved_ones_and_match_the_straight_run`（分岐は main の状態のループを通らないので、
持ち越す解の扱いを足すときにテストも足す）が通ること。持ち越す解をチェックポイントの形式に足す（燃焼は最後の輸送の `Start`、
分岐は軸ごとの `Start`）。

## 経緯

- 2026-10-23: ADR-0196 の実装（feat/lattice-checkpoint-resume）で、並行する perf/moc-warm-start の ADR を読んで起こした
- 2026-10-23: ADR-0195 のマージ（main f3a143d7）の後に rebase し、燃焼は直前の（中点の）輸送の解を、分岐は軸ごとの解と次に解く
  状態をチェックポイントに書くようにした。完了の条件の 2 つのテストと、`checkpoint::tests::saved_states_and_starts_read_back_bit_for_bit`
  が通る（ADR-0194 の Dancoff 係数の持ち越しを含む main の上で、燃焼の再開は時間の内訳を除いてバイトで一致する）
