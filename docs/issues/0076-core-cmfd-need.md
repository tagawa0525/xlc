---
status: open
labels: [core]
opened: 2026-10-12
closed:
milestone: M3
blocked-by: []
resolution:
---

# 炉心の CMFD の要否を Problem 3〜5 の反復の回数と時間で決める

## 問題

炉心の NEM をノードの CMFD で加速する案は、ブランチ feat/core-cmfd の ADR の下書き（proposed、未マージ。2c162bc9）にある。
ユーザーの方針で、この作業は止めている。CMFD なしで VERA Problem 3〜5 を解き、外側の反復の回数と時間を測ってから
要否を決める。Problem 3〜5 の入力（docs/plans/m3-core-hzp.md の項目 8）が先に要る。

feat/core-cmfd の worktree には、コミットしていない実装の途中（crates/numeric/src/cmfd.rs と crates/numeric/src/lib.rs の
1 行）がある。追跡しない OpenMC の出力（examples/mox/openmc-mox{1,2}.msgpack、examples/vera/p2/openmc-2{B,C,D}.msgpack）も
あるが、main の checkout の同じファイルとバイト単位で一致する（2026-10-12 に cmp で確認。追跡しない参照の出力の置き場は #0018）。粗メッシュのソルバを共通のクレートに切り出す計画（docs/plans/m3-core-hzp.md の項目 6）と関係する。
lattice の `cmfd::CoarseMesh::solve` の解き方を変える #0069 と同じソルバに触るので、切り出すときは #0069 の変更に合わせる。

関係する ADR: なし

## 完了の条件

Problem 3〜5 を CMFD なしで解いた反復の回数と時間を記録し、CMFD を入れるかを決める。入れるならその下書きを
書き直してマージし、入れないなら rejected にしてマージし、cmfd.rs で試した切り出しを ADR の「検討した案」かこの issue に書いてから、ブランチと worktree を片付ける。

## 経緯

- 2026-10-12: ユーザーが CMFD を止めると決めた方針と、ブランチに残る作業をメモリから issue に移した
- 2026-11-06: CMFD なしの全炉心の反復の回数と時間（#0156 のブランチ chore/0156-p9-timing。CMFD を入れるかは決めていない）。どれも 1 回の解を 1 スレッドで、
  平らな初期値から解く:
  - Problem 5 の 3D の全炉心の HZP（examples/vera/p5/hzp/core-hzp.log、Slurm 696、b080059c、2026-11-05）: 34 状態の外側の反復は 4,226〜9,215 回、
    span `nem` 42 個（探索の解を含む）の time.busy は 156〜532 秒、全体は `core seconds=1082.67`。パーティション main のジョブで、34 状態を 32 スレッドで
    並列に解いたので、時間は目安
  - Problem 9 の形の全炉心の燃焼（examples/depletion/verification.md の「全炉心の燃焼の時間の内訳（Problem 9 の形）」、35e6eead、Slurm 697 の timed、
    2026-11-06）: 固有値計算 25 回の外側の反復は 4,213〜6,068 回、1 回 132〜203 秒（全体の 99.8 %）。BOC の 1 回は 118 秒、反復 4,165〜4,183 回
    （区分なしの手元の試行で目安）。30 ステップの燃焼は約 2,300 回の固有値計算で、約 113 時間の見積もり（#0156）
