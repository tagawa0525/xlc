---
status: accepted
date: 2026-10-11
supersedes: []
superseded-by:
issues: ["#0018"]
---

# ADR-0126: worktree の references は .worktreelink で main への symlink にして共有する

## 背景

`references/`（核データ、道具、DVC のキャッシュ `references/dvc-cache`）は git で追跡せず、main の worktree にだけ
実体がある（[ADR-0123](0123-drive-stages-with-dvc.md)、[ADR-0124](0124-pin-reference-assets-with-nix.md)）。
linked worktree に `references` がないと、DVC は worktree の中に別のキャッシュを黙って作り、ENDF を deps にする段階は
deps がなくて止まる（ADR-0123 の「結果」）。これまでは docs/environment.md の手順で、人が worktree ごとに symlink を
張っていた。ADR-0123 は、worktree を作る道具で symlink を作る仕組みと、置き場が違えば止める仕組みがまだないと
書いている。

判断の軸:

- worktree を作るたびに人が覚えていなくても、キャッシュと核データが共有されること
- 共有されていない worktree が黙って動かず、止まること
- worktree を作る道具（`~/.claude/scripts/worktree-add.sh`。nixfiles の `.claude/scripts/`）は xlc 以外でも使うので、
  xlc のパスを決め打ちしないこと

## 検討した案

- 案 A: 手で張る（これまで）。忘れても何も止まらず、別のキャッシュができてから気づく
- 案 B: worktree ごとに `dvc config --local cache.dir` で main のキャッシュの絶対パスを書く。キャッシュは共有できるが、
  ENDF の deps（`../../references/endf/...`）は解決せず、`.dvc/config.local` も worktree ごとに手で書くので、
  案 A と同じ忘れ方をする
- 案 C: `references/` を worktree にコピーする。数 GB で、main と同期しない
- 案 D: worktree-add.sh が `references` を決め打ちで張る。汎用のスクリプトに xlc のパスが入る
- 案 E（採用）: リポジトリのルートの `.worktreelink` に共有するパスを書き、worktree-add.sh がそれを読んで main の
  worktree への絶対パスの symlink を作る。加えて、linked worktree の `references` が main を指す symlink でなければ
  `cargo test` を落とす

## 決定と理由

案 E にする。

- `.worktreelink` は 1 行 1 パス（`#` はコメント）で、`references` だけを書く。worktree-add.sh は、そのパスが git で
  無視されていること、main で追跡されていないことなどを確かめてから symlink を作り、`SHARED: references` を出す。
  満たさなければエラーで止まる。仕組みと検討した案は nixfiles の `docs/adr/0016-share-ignored-paths-with-worktreelink.md`
  （名前は PostHog Desktop の `.worktreelink` から借りた）にある
- symlink の先は main の worktree の絶対パス。DVC は symlink をたどって中身の md5 を取るので、deps の md5 は main と
  同じになる（ADR-0124）
- `.worktreelink` のないブランチ（main に rebase する前のもの）や、ほかの方法で作った worktree には symlink が
  できない。そこは crates/docs のテスト `linked_worktrees_share_references` が、`references` が symlink でない、または
  main の `references/` に解決しないときに落ちて止める。テストは、その worktree の `cargo test` と、/git-merge が
  マージの前に worktree で走らせる `run-checks.sh --base main` が走らせる。マージ前の hook は main の worktree で
  走るので、このテストでは止まらない。手で張るコマンドは docs/environment.md に残す
- main の worktree では `references/` は実体なので、テストは確かめない

案 B と D を退けたのは判断の軸のとおりで、案 A と C は共有されない状態を止められない。

## 確認

- `linked_worktrees_share_references`（crates/docs/tests/consistency.rs）。symlink のない linked worktree で落ち、張ると
  通ることを確かめた（2026-10-11、このブランチの worktree）
- worktree-add.sh の `.worktreelink` の挙動は nixfiles の `.claude/tests/scripts/git-workflow.sh`
- このブランチから使い捨ての worktree を作ると `SHARED: references` が出て、symlink が main の `references/` を指し、
  git は無視した（2026-10-11、レビューで確認）

## 結果

- `.worktreelink`、CLAUDE.md の「計算」の節、docs/environment.md の手順
- ADR-0123 の「結果」に書いた「まだない」仕組みは、この ADR で満たす
- `.worktreelink` を足す前に作った worktree（feat/core-cmfd など）は、手で張るまで `cargo test` が落ちる
- `post-merge-cleanup.sh` は変えない。無視された symlink は `git worktree remove` の妨げにならず、消えるのは symlink
  だけで main の実体は残る（nixfiles のテストで確認）
