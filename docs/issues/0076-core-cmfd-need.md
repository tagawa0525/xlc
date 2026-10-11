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
要否を決める。

feat/core-cmfd の worktree には、コミットしていない実装の途中（crates/numeric/src/cmfd.rs と crates/numeric/src/lib.rs の
1 行）がある。追跡しない OpenMC の出力（examples/mox/openmc-mox{1,2}.msgpack、examples/vera/p2/openmc-2{B,C,D}.msgpack）も
あるので、worktree を片付ける前に、ほかに置き場があるかを確かめる。粗メッシュのソルバを共通のクレートに切り出す計画（docs/plans/m3-core-hzp.md の項目 6）と関係する。
lattice の `cmfd::CoarseMesh::solve` の解き方を変える #0069 と同じソルバに触るので、切り出すときは #0069 の変更に合わせる。

関係する ADR: なし

## 完了の条件

Problem 3〜5 を CMFD なしで解いた反復の回数と時間を記録し、CMFD を入れるかを決める。入れるならその下書きを
書き直してマージし、入れないなら rejected にしてブランチと worktree を片付ける。

## 経緯

- 2026-10-12: ユーザーが CMFD を止めると決めた方針と、ブランチに残る作業をメモリから issue に移した
