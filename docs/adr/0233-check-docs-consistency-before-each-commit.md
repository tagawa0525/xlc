---
status: accepted
date: 2026-11-06
requires: []
supersedes: []
superseded-by:
issues: ["#0240"]
---

# ADR-0233: 文書の整合性の検査をコミットごとに pre-commit の hook で実行する

## 背景

文書の整合性の検査（crates/docs/tests/consistency.rs の `docs_are_consistent`。issue と ADR の番号の実在、issue の「関係する ADR:」と
ADR の `issues` の両方向の対応、リンク先のパスの有無など）は、`cargo test` の一部として `run-checks.sh` と、マージの
pre-merge-commit の hook（`run-checks.sh --merge`）だけが実行していた。`.pre-commit-config.yaml` の hook は design-records、
dvc-lock、no-inline-tests、verification-tables、python-ast、references-shell で、この検査は無かった（main e46b8639e）。

そのため、整合しない issue や ADR の変更はコミットを通り、ブランチのチェックかマージで初めて見つかる。2026-11-06 の
11:05〜13:45 JST に、担当の `run-checks.sh` で 2 回（#0236 と ADR-0219、ADR-0231、別のブランチで足した番号 0238 の issue と ADR-0222）、親のマージの hook で 1 回
（#0239 と ADR-0219。`git merge --abort` で取り消した）落ち、ほかに別のブランチにしか無い番号 0238 の issue を引いて「does not
exist」で 1 回落ちた。会話の記録（`~/.claude/projects/-home-tagawa-github-xlc/` の `.jsonl`）の tool_result で、
「lists ADR-NNNN, whose `issues` does not name #NNNN」の違反は 2026-10-29 から 11-06 までに 15 件の issue で出ている。
どれも直すのは数行で、見つかるのが遅いほど、直しのコミット、チェックのやり直し、マージの取り消しが増える。

測った時間（2026-11-06、r995、main e46b8639e から分けた worktree（xlc-docs は変えていない）。`/usr/bin/env time` の経過時間）:

| 実行                                                                    | 時間    |
| ----------------------------------------------------------------------- | ------- |
| `docs_are_consistent` だけ、1 コア（`taskset`）、空の target からの初回 | 3.5 秒  |
| 同じ、CPU を限らない初回                                                | 1.3 秒  |
| 同じ、1 コア、キャッシュあり                                            | 0.11 秒 |
| `cargo test -p xlc-docs` の全体、1 コア、キャッシュあり                 | 0.68 秒 |

xlc-docs は依存のクレートを持たない（crates/docs/Cargo.toml）ので、初回でもビルドは xlc-docs だけである。

## 検討した案

### 検査の範囲

- 案 A: `docs_are_consistent` だけを実行する（採る）。xlc-docs のほかのテストのうち、`production_sources_hold_no_test_code` と
  `record_tables_match_the_compare_outputs` は同じ検査を no-inline-tests と verification-tables の hook が持ち、
  `main_tracks_no_trials` は main でしか意味を持たず、`linked_worktrees_share_references` は worktree の `references` の symlink
  （git で追跡しないもの）を見る。残る 2 つは一時のリポジトリを作る検査のテストで、どれもコミットの中身に依らない
- 案 B: `cargo test -p xlc-docs` の全体。キャッシュありで 0.68 秒と 6 倍で、上の hook と重なる

### 実行するコミット

- 案 A: すべてのコミット（`always_run: true`。採る）。検査はリンク先として追跡しているすべてのファイルのパスを読み、
  消したファイルを指す参照も違反にする。prek は消したファイルを hook に渡さないので、`files` で絞ると、ファイルを消すだけの
  コミットや .md 以外のファイルを改名するコミットで走らない。キャッシュありで 0.11 秒なので、絞って稼ぐ時間は無い
- 案 B: 中身を読む拡張子（`TEXT_EXTENSIONS` の md、rs、py、toml、sh）のファイルが変わるコミットだけ。上の消すだけの
  コミットを見逃す
- 案 C: docs/ の .md が変わるコミットだけ。CLAUDE.md、examples/ の .md、.py、.rs、.toml、.sh、crates/ の .rs と .toml の issue と
  ADR の番号も検査の対象（crates/docs/src/check.rs の `scanned_for_refs`）なので、それらの変更を見逃す

### 排他の Slurm のジョブの間の CPU

hook の中の cargo は、claude-hooks の guard-exclusive-slurm（Bash のツールのコマンドを見る）に見えない。

- 案 A: ジョブの有無によらず、いつも 1 コアに限る（採る）。初回で 2.2 秒、キャッシュありで差は無い。CPU は、このプロセスに
  許された CPU（`/proc/self/status` の `Cpus_allowed_list`）の最初の 1 つにする。Slurm の `-c` のジョブの中のように
  許された CPU が限られていると、固定の番号では `taskset` が失敗する
- 案 B: `run-checks.sh` と同じく、`~/.claude/scripts/lib/slurm-exclusive.sh` を source し、排他のジョブがあるときだけ限る。
  判定は claude-hooks と揃うが、リポジトリの hook が `~/.claude` の生成物に依り、コミットのたびに squeue（上限 5 秒）を
  呼ぶ。限らずに稼げるのは初回の 2.2 秒だけである
- 案 C: 限らない。排他のジョブの時間の記録を乱しうる

## 決定と理由

各案の A を採る。`.pre-commit-config.yaml` に hook `docs-consistency` を足し、`always_run: true` と `pass_filenames: false` で、
`scripts/check-docs-consistency.sh` が、許された CPU の最初の 1 つに `taskset` で限って
`cargo test -q -p xlc-docs --test consistency -- --exact docs_are_consistent` を実行する。

違反をコミットの時点で止めれば、直しは同じコミットに入り、チェックやマージで見つかってやり直す周が無くなる。prek は
ステージしていない変更を退避してから hook を実行し、検査は git が渡す index（`git ls-files`）と作業ツリーを読むので、
ステージした中身を見る（グローバルの pre-commit は prek の後に Markdown の書式を自動で直してステージし直すので、その直しは
検査の後に入る）。`run-checks.sh` と pre-merge-commit の `cargo test` の検査は残す（`--no-verify` のコミットと、
別々のブランチの変更が合わさってはじめて崩れる参照を拾う）。

### 帰結

- 整合しない issue と ADR の変更、ほかのブランチにしか無い issue や ADR を引く変更は、コミットで止まる。直す側は、
  違反の行（`<パス>:<行>: …`）を hook の出力で受ける
- 新しい worktree の最初のコミットは、xlc-docs のビルドで 3.5 秒ほど長くなる。2 回目からは 0.1 秒ほど
- 排他のジョブの間も、hook は 1 コアしか使わない
- rebase の後の片付けのコミットが止まりうること、CPU を限る規則が `run-checks.sh` と違うこと、マージを締めるコミットで検査が
  2 度走ることは #0240 に残した

### 確認

- 2026-11-06（20c956ff1 の上の、ブランチに残さない一時のコミット）: #0239 の「関係する ADR:」を ADR-0222 にして ADR-0222 の
  `issues` を変えないコミットは「lists ADR-0222, whose `issues` does not name #0239」で止まった（`git commit -a` でも同じ）。
  ADR-0222 の `issues` に #0239 を足してステージしないまま残したコミットも止まり（ステージした中身を見る）、ステージすると通った
