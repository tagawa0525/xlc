---
status: accepted
date: 2026-10-31
requires: []
supersedes: []
superseded-by:
issues: ["#0184", "#0210"]
---

# ADR-0223: 生の dvc checkout を xlc のプロジェクトの hook で止め、scripts/dvc-checkout.py を使わせる

## 背景

生の `dvc checkout` は、`dvc.lock` に hash がまだ無いステージ（流していない参照解）の出力が戻らないため、main の上でも
いつも終了コード 255 で終わる（#0184）。hash があるのにキャッシュに無い出力（本当の失敗）も同じ「Checkout failed」の一覧に
並ぶので、終了コードでは見分けられず、一覧を目で確かめることになる。#0184 で `scripts/dvc-checkout.py`（未実行の出力を
数えるだけにし、hash があるのに戻らない出力があれば終了コード 1 で止まる）を作り、/dvc-relock と /topic-start、
`scripts/worktree-restore-outputs.sh` から呼ぶようにした。

それでも、/retro（35 回目）が数えたところ、親のセッションは直近の区間で 3 回中 2 回、過去 8 セッションで約 200 回、生の
`dvc checkout` を Bash ツールで実行し、「Checkout failed」の 6 件を目で確かめていた。スキルの文言はスクリプトを指して
いたが、守られなかった。ユーザーは「hookを追加したらいいと思います。」（ユーザー 2026-10-31）と決めた。

判断の軸:

- 生の `dvc checkout` を、覚えていなくても実行できないこと（止めたときに使うものが分かること）
- 1 つのリポジトリだけの手順は、そのリポジトリに置くこと
- `grep 'dvc checkout'` や文書の中の語など、dvc を実行しないコマンドを止めないこと

## 検討した案

- 案 A（採用）: xlc の `.claude/settings.json`（git で追跡）に PreToolUse の Bash の hook を足し、`.claude/hooks/` の Python の
  スクリプトがコマンドを読んで、dvc の checkout を直接呼ぶものを終了コード 2 で止める。xlc にはプロジェクトの hook がまだ無く、
  これが最初の例になる
- 案 B: nixfiles の claude-hooks（全リポジトリで動く共通の hook）に「`scripts/dvc-checkout.py` のあるリポジトリでは止める」規則を
  足す。共通の hook に 1 つのリポジトリの事情（スクリプトの名前と使い方）を入れることになり、xlc でスクリプトの形を変えると
  別のリポジトリを直すことになる
- 案 C: スキルの文言だけで守る。/dvc-relock と /topic-start は既にスクリプトを指していたが、上の背景のとおり止まらなかった

## 決定と理由

案 A にする。

- 置き場: `.claude/settings.json` の `hooks.PreToolUse` に `matcher: "Bash"` の command の hook を置き、`.claude/hooks/block-raw-dvc-checkout.py`
  を `"${CLAUDE_PROJECT_DIR}"` から呼ぶ。スクリプトは標準ライブラリだけで書き、`nix develop` の外の PATH の python3 で動く
  （hook はシェルから実行され、devShell に入っていない）。home の settings の PreToolUse（`claude-hooks pre-tool-use`）とは
  別の設定の階層の hook で、両方が走る（Claude Code は設定の階層の hook を置き換えずに足し合わせる）
- 判定: stdin の JSON の `tool_input.command` を、シェルの引用（`'…'`、`"…"`、`\`）、コメント、ヒアドキュメントの本文を
  読み分けて単純なコマンドに分け（`;`、`&&`、`||`、`|`、`&`、改行、`(…)`、`$(…)`）、各コマンドの先頭の代入と包むコマンド
  （`env`、`time`、`nice`、`timeout`、`taskset`、`nix develop`/`nix shell` の `-c`/`--command`、`bash -c` などの文字列）を
  外した先が `dvc`（または `python -m dvc`）で、最初のサブコマンドが `checkout` なら止める。引用の中の語（`grep 'dvc checkout'`、
  `echo "dvc checkout"`）、コメント、ヒアドキュメントの本文は実行されないので止めない
- 拒否の文言: 理由（未実行の出力でいつも 255 になり、本当の失敗が埋もれる）と、代わりに使う
  `nix develop -c scripts/dvc-checkout.py [<dvc.yaml か .dvc の指定>...]` を stderr に出す
- 読めないコマンド（閉じていない引用）は、止めずに終了コード 1 で知らせる（bash もそのコマンドを実行できないため）

理由:

- 案 C は既に試して守られなかった。止めるには、コマンドを実行する前に機械で判定する hook が要る
- 案 B と比べ、スクリプトの名前と使い方はこのリポジトリの事情で、xlc の中に置けば `scripts/dvc-checkout.py` と同じブランチで
  一緒に変えられる
- 語の検索（`dvc checkout` の文字列を含むか）では、スキルを直す作業の `grep 'dvc checkout'` やコミットの本文を止めるので、
  シェルの引用を読み分ける

### 帰結

- 生の `dvc checkout` は、Claude Code の Bash ツールからは実行できない。スクリプトの中（`scripts/dvc-checkout.py` 自身など）の
  `dvc checkout` は Bash ツールのコマンドに現れないので止めない
- 判定はシェルの文法の一部だけを読む。`"…"` の中の `$(…)` は読むが、バッククォートの中、`slurm-run.sh` などに文字列として渡す
  コマンド、`eval` は読まない（止め損なうが、誤って止めることはない）
- スキル（/openmc-reference）と docs の手順に残っていた生の `dvc checkout` のコマンドを、スクリプトの形に直した
- examples/ の記録の手順、dvc.yaml のコメント、スクリプトのエラーの文言に残る生の `dvc checkout` は #0210

### 確認

- `.claude/hooks/test_block_raw_dvc_checkout.py` が、止めるコマンドと通すコマンドの例、hook のプロセスの終了コードと stderr を
  確かめる（`scripts/run-python-tests.sh` が実行する）
- Claude Code の Bash ツールで実際に止まるかは、マージの後に親のセッションで確かめる
