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
1 周目で、再現を書かず直さなかった指摘（/git-merge の「直す指摘の基準」の non-blocking）。

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

## 完了の条件

1 の手順を /git-merge で決めるか、起きないことを確かめる。2 を 1 つの規則にそろえるか、そろえない理由を書く。3 を許すか直すかを決める。

## 経緯

- 2026-11-06: 担当の 1 周目（HEAD c1a7b3840。後で ADR を実装の前に並べ直したので、今のブランチの祖先ではない）の指摘から起こした
