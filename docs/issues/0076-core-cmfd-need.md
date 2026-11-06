---
status: closed
labels: [core]
opened: 2026-10-12
closed: 2026-11-06
milestone: M3
blocked-by: []
resolution: ADR-0232
---

# 炉心の CMFD の要否を Problem 3〜5 の反復の回数と時間で決める

## 問題

炉心の NEM をノードの CMFD で加速する案は、ブランチ feat/core-cmfd（2026-11-06 に消し、タグ archive/feat-core-cmfd で残す）の ADR の下書き（proposed、未マージ。de5b93436）にあった。
ユーザーの方針で、この作業は止めている。CMFD なしで VERA Problem 3〜5 を解き、外側の反復の回数と時間を測ってから
要否を決める。Problem 3〜5 の入力（docs/plans/m3-core-hzp.md の項目 8）が先に要る。

feat/core-cmfd の worktree（2026-11-06 に消した）には、コミットしていない実装の途中（crates/numeric/src/cmfd.rs と crates/numeric/src/lib.rs の
1 行）があった。追跡しない OpenMC の出力（examples/mox/openmc-mox{1,2}.msgpack、examples/vera/p2/openmc-2{B,C,D}.msgpack）も
あったが、main の checkout の同じファイルとバイト単位で一致する（2026-10-12 に cmp で確認。追跡しない参照の出力の置き場は #0018）。粗メッシュのソルバを共通のクレートに切り出す計画（docs/plans/m3-core-hzp.md の項目 6）と関係する。
lattice の `cmfd::CoarseMesh::solve` の解き方を変える #0069 と同じソルバに触るので、切り出すときは #0069 の変更に合わせる。

関係する ADR: ADR-0232

## 完了の条件

Problem 3〜5 を CMFD なしで解いた反復の回数と時間を記録し、CMFD を入れるかを決める。入れるならその下書きを
書き直してマージし、入れないなら rejected にしてマージし、cmfd.rs で試した切り出しを ADR の「検討した案」かこの issue に書いてから、ブランチと worktree を片付ける。

## 経緯

- 2026-10-12: ユーザーが CMFD を止めると決めた方針と、ブランチに残る作業をメモリから issue に移した
- 2026-11-06: CMFD なしの全炉心の反復の回数と時間（#0156 のブランチ chore/0156-p9-timing。CMFD を入れるかは決めていない）。どれも 1 回の解を 1 スレッドで、
  平らな初期値から解く:
  - Problem 5 の 3D の全炉心の HZP（examples/vera/p5/hzp/core-hzp.log、Slurm 696、3575aff8、2026-11-05）: 34 状態の外側の反復は 4,226〜9,215 回、
    span `nem` 42 個（探索の解を含む）の time.busy は 156〜532 秒、全体は `core seconds=1082.67`。パーティション main のジョブで、34 状態を 32 スレッドで
    並列に解いたので、時間は目安
  - Problem 9 の形の全炉心の燃焼（examples/depletion/verification.md の「全炉心の燃焼の時間の内訳（Problem 9 の形）」、ca59c64b、Slurm 697 の timed、
    2026-11-06）: 固有値計算 25 回の外側の反復は 4,213〜6,068 回、1 回 132〜203 秒（全体の 99.8 %）。BOC の 1 回は 118 秒、反復 4,165〜4,183 回
    （区分なしの手元の試行で目安）。30 ステップの燃焼は約 2,300 回の固有値計算で、約 113 時間の見積もり（#0156）
- 2026-11-06: CMFD を入れると決めた（ADR-0232。全炉心の燃焼の速度の対策を Xe とホウ素の同時の収束、CMFD、複数のスレッドの順に進める判断。ユーザー
  2026-11-06 の判断の要約）。下書き（feat/core-cmfd の 0119 番の ADR）を書き直し、旧 worktree の `crates/numeric/src/cmfd.rs` の切り出しの試しを
  ADR-0232 の検討した案（案 S1）に書いた。ソルバは共通のクレート `xlc-cmfd` に置いた（lattice は #0241）。反復の回数と時間は
  examples/depletion/verification.md の「固有値計算の CMFD（小さな炉心と全炉心）」: Problem 9 の形の全炉心の燃焼（Slurm 704 の timed、77a720cd6）で、
  全体 4,664.7 s から 374.2 s、外側の反復の和 151,099 回から 4,638 回、固有値計算 1 回の外側の反復 27〜384 回。残りはブランチ feat/core-cmfd と
  worktree の片付けで、マージの後に行う
- 2026-11-06: 「古いCMFDは消しても良い。」（ユーザー 2026-11-06）を受け、親が worktree xlc-feat-core-cmfd とブランチ feat/core-cmfd を消し、
  先端 4d82b4ce8 にタグ archive/feat-core-cmfd を付けた（ADR-0232 の検討した案が引く de5b93436 を履歴に残すため。docs/branches.md）。
  未コミットの `crates/numeric/src/cmfd.rs` と `lib.rs` の変更はタグに入らない。そこから取り込んだのは、lattice の粗メッシュの固有値問題の
  ソルバ（ブロック Gauss-Seidel のべき乗法）を切り出す試し（ADR-0232 の案 S1）で、置き場を新しいクレート `xlc-cmfd` に変えて core が使う
  （案 S3。6bf685f7c のマージの本文の「旧 worktree の未コミットの cmfd.rs から取り込んだ形」。lattice は #0241）。完了の条件のうち CMFD なしの反復の回数と時間は Problem 5 の HZP と Problem 9 の形の燃焼で記録し（上の 2026-11-06 の項。
  Problem 3 と 4 の値はこの issue に無い）、CMFD を入れると決め（ADR-0232）、切り出しの試しを案 S1 に書き、ブランチと worktree を片付けたので閉じる
