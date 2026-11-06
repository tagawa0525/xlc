---
status: open
labels: [tooling, docs]
opened: 2026-11-06
closed:
milestone:
blocked-by: []
resolution:
---

# 文書の整合性の hook のレビューで残した指摘

## 問題

ブランチ chore/0233-docs-consistency-hook（ADR-0233 の pre-commit の hook `docs-consistency`）に担当が回したローカルのレビューの
1〜3 周目で、再現を書かず直さなかった指摘（/git-merge の「直す指摘の基準」の non-blocking）。

関係する ADR: ADR-0233

1. rebase の後の片付けのコミットが hook で止まりうる。main とブランチが合わさってはじめて崩れる参照（main が ADR の `issues` から
   ブランチの issue を外した、ブランチがリンクするファイルを main が消したなど）が rebase の後の木にあると、post-rewrite の hook が
   残した hash の付け替え（`chore: remap the commit refs after the rebase`）のコミットはその違反で止まる。先に違反を直すコミットは、
   prek が付け替えを退避するので、ブランチの issue の `resolution` が rebase の前のコミットを指したままになり、止まる。両方を 1 つの
   コミットに入れるしか道が無い。番号の重複（design-records の hook が止める）と付け替えが重なったときも、`design-records renumber`
   が汚れた作業ツリーで止まるので、同じ形で詰まる（この hook の前からある）。再現は確かめていない。扱いは /git-merge の
   「番号の重複（rebase の後）」の手順の側（nixfiles）で決める
2. `scripts/check-docs-consistency.sh` は許された CPU の最初の 1 つ（r995 では 0）に限り、`run-checks.sh` はオンラインの CPU の
   最後の 1 つ（`getconf _NPROCESSORS_ONLN` − 1）に `-j 1` で限る。排他のジョブの間に検査を置く規則が 2 つある。CPU 0 は割り込みを
   受けやすいという見方は確かめていない
3. マージを締めるコミット（競合の解決の後の `git commit`）では、グローバルの pre-commit が pre-merge-commit（`run-checks.sh --merge`
   の `cargo test`。`docs_are_consistent` を含む）を実行した後に prek がこの hook を実行し、同じ検査が 2 度走る（キャッシュありで
   0.11 秒）。同じ worktree で別の cargo のビルドが走っている間は、`git commit` が target のロックを待つ
4. `--exact docs_are_consistent` が 0 本に当たっても `cargo test` は 0 で終わる（レビューが `--exact no_such_test` で確かめた）。
   テストを改名するか `#[ignore]` を付けると、hook は何も検査せずに通る
5. テストが読むリポジトリの根はコンパイルの時の `env!("CARGO_MANIFEST_DIR")` で決まる。別の worktree の target を `CARGO_TARGET_DIR` で
   共有すると、その worktree で作ったテストのバイナリが再利用され、別の worktree の作業ツリーを読む（1 周目のレビューが
   この worktree の target を共有して起き、テストが消えた worktree を指して落ちた）。`run-checks.sh` の `cargo test` も同じ
6. hook の `cargo test` に `--locked` が無い。Cargo.toml の依存の変更だけをステージして Cargo.lock を残すと、cargo が Cargo.lock を
   書き換え、prek がファイルの変更として hook を失敗にする
7. tests/consistency.rs に、まだ無い API を呼ぶ RED のテストをコミットすると、コンパイルが通らずに hook が止まる（実行して落ちる RED は、
   `docs_are_consistent` だけを実行するので止まらない。crates/docs/src の単体テストはこのターゲットでコンパイルしない）
8. ADR-0233 の時間は手元の `/usr/bin/env time` の測定で、Slurm の timed のジョブではない（依頼が手元で測ると決めた）
9. `pass_filenames: false` と `always_run: true` で hook は 1 回しか呼ばれないので、`require_serial: true` は効かない
10. 5 は、共有した先の worktree が残っていれば落ちずに通り、違反を見逃す（3 周目のレビューが小さな crate を 2 つのコピーで
    target を共有して再現した）。ADR-0233 の帰結は、4 と 5 の、検査が黙って効かなくなる道を挙げていない
11. ADR-0233 の案 A の「どれもコミットの中身に依らない」は、直前の「残る 2 つ」を指すが、段落の全体を指すとも読める
    （`production_sources_hold_no_test_code` と `record_tables_match_the_compare_outputs` は中身を読む）

## 完了の条件

1 の手順を /git-merge で決めるか、起きないことを確かめる。2 を 1 つの規則にそろえるか、そろえない理由を書く。3〜11 を直すか、直さない理由を書く。

## 経緯

- 2026-11-06: 担当の 1 周目（HEAD c1a7b3840。後で ADR を実装の前に並べ直したので、今のブランチの祖先ではない）の指摘から起こした
- 2026-11-06: 担当の 2 周目（HEAD 023d127fc）の指摘から 4〜9 を足した。競合を解いた後の `git rebase --continue` でも hook が走るという
  指摘は、一時のリポジトリで確かめると pre-commit が呼ばれなかったので外した
- 2026-11-06: 担当の 3 周目（HEAD 1bad0bc92）の指摘から 10 と 11 を足した
