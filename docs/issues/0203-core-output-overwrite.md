---
status: open
labels: [core, lattice, xs, tooling]
opened: 2026-10-30
closed:
milestone:
blocked-by: []
resolution:
---

# core、lattice、xs は出力のファイルが別のプロセスに書かれていても黙って上書きする

## 問題

core、lattice、xs は、出力のファイルを別のプロセスが書いているかを見ずに開き、黙って上書きする。同じ出力のパスを
指す 2 本のジョブが重なると、片方の結果がもう片方に混ざるか消える。

2026-10-29 に #0200 の担当が起こした例（Slurm 600、601、602 の時刻は `sacct -j 600,601,602` で確かめた）:

- Slurm 600（`xlc-small-core-equilibrium-splits`、2026-10-29 18:07 〜 20:35）は、平衡 Xe（ADR-0219）で目標の k を
  1.075 にした `examples/depletion/branches/core-small-equilibrium-<分割>.toml` を流していた。
- その間に、目標の k だけを 1.071 に変えた同じ名前のケースを Slurm 601（`xlc-small-core-k1071-splits`、19:13 〜 21:42）で、
  同じ `.out`、`.log`、`.msgpack` に流した。600 を取り消そうとしたが、`scancel` が auto mode の分類器に拒まれた。
- 2 本が同じ出力に書き、601 の平衡のログに NUL が入った。担当は行数を確かめて NUL を除き、履歴は 601 が上書きしたことを
  `states` の k で確かめた。その後、Slurm 602（`xlc-small-core-equilibrium-final`、22:14 〜 翌 00:41）で全部を作り直したので、
  記録の値は無事だった。

出典は #0200 の担当の最終の報告（2026-10-30 00:54）。報告の内容のうち、NUL の除去と `states` の確認は、
ここでは確かめていない。

2 本が重なったのは、排他のジョブでなかったため: `sacct` では 600 が 5 CPU（割り当て 6）、601 が 16 CPU で、どちらも
パーティション main、ノード r995。`-c` の付いた排他なしのジョブ（CLAUDE.md の「計算」）で、同じ出力のパスを
2 本が使うことを、投げる側が見ていない。なぜ排他なしで投げたかは確かめていない。

コードの現状（2026-10-30 の main）: crates の src に `flock` はない。書き方は 3 通りある。

- core: `crates/core/src/main.rs` の `std::fs::write` 2 箇所（133、272 行）と、xs の `crates/xs/src/main.rs` の
  `std::fs::write`（358 行）。開くときに切り詰める
- lattice の `write_synced`（`crates/lattice/src/checkpoint.rs` 406 行）: `File::create` で切り詰める
- lattice のチェックポイントの `write`（同 198 行）: 固定名の `<path>.tmp` に書いて rename で置き換える

どれも、書く先が別のプロセスに開かれているかを見ない。

候補の対策（どちらを採るかは決めていない）:

- (a) 出力の msgpack に排他の `flock` を取り、取れなければエラーで止める。シェルのリダイレクトで書く `.out` と `.log` は
  守れない（今回 NUL が入ったのはログ）。実装の制約: 切り詰める書き方（`File::create`、`std::fs::write`）は、ロックを
  取る前に中身を消すので、切り詰めずに開いて（`OpenOptions`）ロックを取ってから切り詰める。rename で置き換える
  チェックポイントは、ロックが古い inode に付くので、対象の出力でなく別のロック用のファイルに取り、`<path>.tmp` が
  固定名で 2 本が衝突することも扱う。共有のファイルシステムでの `flock` の効き方は確かめていない。
- (b) 条件を変えたらケースの名前を変える規則を、/topic-start か CLAUDE.md の「計算」に書く。コードは変えないので、
  規則を破れば同じことが起きる。

関係する ADR: なし

関係する issue: #0200、#0202

## 完了の条件

対策を ADR で決めるか、規則を書く。(a) を採るなら、守る出力（core、lattice、xs の msgpack）の範囲と、ログを
守れないことの扱いを決め、テストで、2 本目がエラーで止まることを示す。(b) を採るなら、規則を書いて閉じる。

## 経緯

- 2026-10-29: #0200 の担当が、Slurm 600 と 601 が同じ出力に書いて 601 のログに NUL が入ったことを報告した。
- 2026-10-30: issue に起こした。#0200、#0202 の側には、この issue への参照を足していない。
