---
status: open
labels: [core, depletion, verification]
opened: 2026-11-06
closed:
milestone:
blocked-by: []
resolution:
---

# 同時の反復のレビューで残した記録と入力の指摘

## 問題

ブランチ perf/0156-xenon-boron-joint（ADR-0219 の Xe とホウ素の同時の反復）に親が回したローカルのレビューの 3 周目で、再現を書けず直さなかった指摘
（/git-merge の「直す指摘の基準」の non-blocking）。

関係する ADR: なし

1. 「固有値計算の初期値（小さな炉心）」の `compare-core-small-guess.json` が引く `core-small-equilibrium-1.log`（md5 02a10843）は、ADR-0231 の
   ブランチの worktree の後片付けで消えていて（`dvc add` していなかった）、同時の反復のブランチの `core-small-equilibrium-1.log.dvc` が追うのは
   同時の反復のログ（md5 96d9f7c1）。履歴（md5 0c3828ae）は `core-small-equilibrium-nested-1.msgpack.dvc` で戻るが、ログは戻らない。
   比較の JSON が引く生成物の棚卸しは、ADR-0222 の対象を広げるブランチ docs/0222-compare-json-generated の issue の対象
2. 同時の反復では `EQUILIBRIUM_SOLVES`（50）が Xe とホウ素を合わせた反復の上限になり、組成が釣り合ってホウ素だけが動いているときも
   「the equilibrium xenon did not converge in 50 solves」で止まる。文言がホウ素の残差を示さない（解き直すたびに探す形は Xe の解き直し 50 回 × 探索の
   `max-solves`）。そうなる入力は見つかっていない
3. `compare_core_xenon_boron.py` の `not all(str(p) for p in paths)` は `str(Path(""))` が `.` なので空のパスを止めず、`IsADirectoryError` になる
   （compare_core_guess.py のコピー。#0236 の重複の対象）
4. `xenon-boron = "nested"` は、`search` の状態があれば受け付け、efpd のステップがその状態を解かないケースでは黙って効かない（`midpoint-boron = "start"` は
   同じ条件をエラーにする）
5. ADR-0219 の「確認」と記録の小さな炉心の固有値計算の時間（454.2 s → 137.2 s）は、Slurm の外で 2 本を同時に流した手元の目安で、記録は「目安」と
   書くが ADR の行の注は短い（#0234 の 4 と同じ扱い）
6. 記録の再現のコマンドは 2 本を順に流す `for` だが、「計算」は `taskset` で別のコアに置いて同時に流したと書く。時間の再現の条件がコマンドから分からない

## 完了の条件

1 を docs/0222-compare-json-generated の棚卸しに引き継ぎ、2〜6 を直すか、直さない理由を書く。

## 経緯

- 2026-11-06: 親の 3 周目（HEAD 33749b7b7）の指摘から起こした。関係する ADR を「なし」にしたのは、ADR-0219 の `issues` に入れると ADR の変更になり、
  マージの前のレビューをやり直すため（#0236 と同じ扱い）。1 周目の指摘は親が直し（#0236 の経緯）、2 周目の指摘も直した
