---
name: topic-start
description: xlc で 1 つのトピックの作業を始め、マージするまでの手順。記録の番号、worktree、references の symlink、生成物のコピー、ビルド、TDD のコミットの分け方、記録に書く計算の前の main への rebase、チェック、/git-merge の呼び方、マージ前後の片付け。
argument-hint: <branch>
---

# トピックを始めてマージする

1 つのブランチは 1 つのトピックで閉じる（/topic-triage）。xlc には GitHub のリモートがないので、PR の代わりに
/git-merge でマージする。CLAUDE.md の「作業の流れ」にある規則は繰り返さず、その手順の細部を書く。

## 始める

```sh
# 1. 記録の番号。ADR か issue が要るときだけ（スキル adr。docs/issues/README.md）。
#    手で最大 + 1 を振らない。タグで予約するので、ほかの worktree と番号が重ならない
design-records new <slug> "<決めたことの 1 文>"            # ADR（proposed で始まる）
design-records new --kind issue <slug> "<題名>"           # issue（題名の行だけ。フロントマターと節は手で埋める）

# 2. worktree。main から上流なしで作る。ブランチ名に issue の番号を入れる（例 fix/0004-frame-gap）
~/.claude/scripts/worktree-add.sh <branch>
```

`design-records new` は、書きたい worktree の中で実行する（実行した場所にファイルができる）。先に worktree を作って、そこで
実行する。

worktree の出力に `SHARED: references` が出ていることを確かめる。出ていなければ、`ls -l references` で main への
symlink かを見て、なければ worktree の直下で手で張る（ADR-0126）:

```sh
ln -s "$(git rev-parse --path-format=absolute --git-common-dir)/../references" references
```

symlink がないと、DVC は worktree の中に別のキャッシュを黙って作る。`cargo test` の
`linked_worktrees_share_references` も落ちる。

## 生成物を main から持ってくる

git で追跡しない生成物は worktree に無い。使うものを main からコピーする（btrfs なので reflink で容量を食わない）:

```sh
cp --reflink=auto ~/github/xlc/examples/vera/vera.msgpack examples/vera/
cp --reflink=auto ~/github/xlc/examples/vera/p2/openmc-2K.msgpack examples/vera/p2/    # 比べるツールが読む参照解のピン出力（ADR-0121）
dvc checkout dvc.yaml                                                                 # 段階が呼ぶ bin/ の実行ファイル（ADR-0159）
dvc checkout examples/vera/dvc.yaml                                                   # DVC のステージの出力が要るとき（openmc-library/ など）
```

そのあと、worktree でビルドする。`target/` は worktree ごとにあり、最初は全体を作るので数分かかる（手元でよい。ADR-0026）。
DVC のステージが呼ぶ実行ファイルは、これとは別に `dvc repro` の `build` が `target/dvc` に作る（ADR-0159）:

```sh
cargo build --release
```

ワークスペースの全体のビルドは、`build` の 3 つのパッケージだけのビルドと共有の依存の feature の統合が違いうるので、
`target/release` の実行ファイルは `bin/` と同じバイトとは限らない。段階と同じ機械語が要るとき（perf、段階のパニックの
再実行）は `cargo build --release -p xlc-xs -p xlc-lattice -p xlc-core` で作り直す（ADR-0163）。

## TDD のコミット

テストを先に書いて、失敗することを確かめてからコミットする。実装は別のコミットにする。

- RED と GREEN の件名は別の観点から書く。テストは期待すること（`test(core): expect f* sides facing reflector sites to be errors`）、
  実装はどうするか（`fix(core): reject f* sides that face reflector sites`）。「test: X」と「feat: X」を並べない
- 本文に理由を書き、決定は ADR に書いて `Refs: ADR-NNNN`（issue は `Refs: #NNNN`）で引く。件名と本文は英語
- 実装を始めたあとでテストだけをコミットするとき（テストが実装と同じファイルの `mod tests` にある場合など）は、テストの
  部分だけをステージする。`git stash` は使わない（ほかのセッションと共有で、取り違える）。`git add -p` は対話なので、
  hunk の順に `y`/`n` を標準入力で渡す: `printf 'n\ny\nn\n' | git add -p <file>`。そのあと `git diff --cached` で、
  テストの hunk だけが入っていることを確かめてからコミットする
- ドキュメントと設定だけの変更は、TDD の対象外

## 記録に書く計算の前に main に rebase する

記録（examples/*/verification.md、dvc.lock）に書く数値を出す計算を Slurm に流す前に、ブランチを main に rebase する。
計算のあとに rebase すると、記録の hash の付け替え（commit-refs remap）と dvc.lock の整え直し（dvc-relock）が要り、
場合によっては計算のやり直しになる。先に rebase すれば、マージ時の rebase で記録の hash が変わる場面を減らせる。

1. コードと入力の変更をすべてコミットする。記録に書く試行の入力（`trials/`）は `git add -f` でコミットする（スキル verification-record の「記録に書く試行」）。`git status --short` が空であることを確かめる
2. `git rebase main`。衝突があればここで解く。`dvc.lock` が衝突したら、rebase を終えてからスキル dvc-relock の前半で整える（rebase の途中では dvc を実行しない）。rebase の終わりに post-rewrite hook が記録の hash を付け替えたら、その差分をコミットしてから次へ進む（スキル verification-record の「rebase のあと」）
3. 段階（`dvc.yaml`）で実行するなら、先に `dvc checkout`（スキル dvc-relock の前半 1）。チェックが見る出力と計算の入力を揃えるため、チェックより前に行う
4. `git status --short` が空であることを確かめてから、`~/.claude/skills/language-checks/scripts/run-checks.sh -C <worktree> --base main` を通す（マージ時と同じ、main との分岐点からの差分の検査）
5. 記録用の計算を流す（スキル verification-record の「数値を出す前」）。記録に書く hash は、計算を流したときの HEAD（`dvc.lock` をコミットする前のコミット）、分岐点は `git merge-base HEAD main`

計算の待ち時間が長く、その間に main が進んだら、マージ前の rebase で hash が変わる。計算は測り直さず、付け替えは
commit-refs remap に任せる（ADR-0122。コード、入力、環境が変わらなければ測り直さない）。長い計算が複数のブランチで
並ぶときも、流す直前に rebase する。

ジョブは開始したときに worktree のコードや入力を読むので、開始の時刻が分からないあいだは、その worktree を変えない（rebase、
checkout、編集、`cargo build`）。待っている間に進めたい作業は別の worktree で行う。ジョブが終わったら、流した HEAD を
記録に書いたまま、`dvc.lock` と記録の差分を別のコミットにし、そのあと rebase する。main の変更がコード、入力、環境を変えて
いたら、付け替えではなく計算を再実行する（ADR-0122）。

## チェック

```sh
~/.claude/skills/language-checks/scripts/run-checks.sh
echo "rc=$?"
```

終了コードで判定する。`run-checks.sh | tail && git merge` のようにつなぐと、失敗してもマージに進む。`cargo test` は
`design-records check` と `linked_worktrees_share_references` も実行する。

## マージ

マージの前に、残さない試行を消す（`trials/`。スキル verification-record の「記録に書く試行」）。ADR を実装に合わせて直し
`status: accepted` にする、issue を閉じる（`status: closed`、`closed`、`resolution`）のも、このブランチの中で行う。

worktree から /git-merge を呼ぶときは、**ブランチ名を引数に渡す**: Skill ツールで `git-merge` を `<branch>` の引数つきで
呼ぶ。引数なしだと、Bash の cwd が呼び出しごとに main の worktree に戻るので、対象が main になり、
「既定ブランチ自身はマージの対象にしません」で止まる。

- 番号の重複は、/git-merge が rebase の直後に `design-records check` で見つけ、`design-records renumber`（issue は
  `--kind issue`）で振り直す。手で付け替えない。`renumber` は作業ツリーに追跡していない変更があると止まる
- rebase のあと、記録のコミットの hash は hook が付け替える（スキル verification-record）。付け替えた差分をコミットして
  から、チェックに戻る

## 後片付け

- 後片付け（worktree とブランチの削除、予約のタグの解放）は /git-merge が行う
- 核種を足すなどして作り直した `examples/vera/vera.msgpack` は、main の同じ場所にも置く。ほかの worktree は main の
  ものをコピーして使うので、置き換えないと古いライブラリで計算が始まる。追跡しない生成物なので、main の作業ツリーは
  汚れない。/git-merge は worktree を消すので、マージの前に（rebase とチェックのあと）コピーしておく:

  ```sh
  cp --reflink=auto examples/vera/vera.msgpack ~/github/xlc/examples/vera/vera.msgpack
  ```

  DVC のステージの出力を入れたブランチをマージしたあとは、main で `dvc checkout dvc.yaml` と `dvc checkout examples/vera/dvc.yaml`
  をして `dvc.lock` の出力に合わせる（ADR-0130、ADR-0159。手順はスキル dvc-relock の後半）
