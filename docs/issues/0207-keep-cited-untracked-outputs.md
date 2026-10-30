---
status: open
labels: [tooling, dvc]
opened: 2026-10-30
closed:
milestone:
blocked-by: []
resolution:
---

# 記録が md5 で引く追跡しない出力が worktree の後片付けで消える

## 問題

`examples/depletion/micro/*.msgpack` などは DVC の出力でなく、`.gitignore` で無視する生成物で、検証の記録と ADR が md5 で
引く。作ったブランチの worktree にだけあり、/git-merge の後片付け（`post-merge-cleanup.sh` の `git worktree remove`）が
無視されたファイルごと消す。DVC の出力はキャッシュから戻り、`examples/vera/p2/openmc-2K.msgpack` は
`scripts/worktree-restore-outputs.sh` が main からコピーするが、ほかの追跡しない出力には戻す経路が無い。

- 2026-10-30 に、ADR-0221 の依頼が「既存の `shifted-2A.msgpack`（md5 456654d5）をそのまま使う」と書いたが、どの worktree にも
  残っておらず、担当は置き換えた組成のテーブルを Slurm 618（20 分 9 秒）で作り直した（md5 は一致）
- ADR-0221 のマージ（xlc 9a69db2a9）の後片付けでも、記録と ADR が引く `table-nominal.msgpack`（md5 fde256d9）などの micro の
  テーブルが消え、main にも無い（2026-10-30 に `/home/tagawa/github` の下を `find` で探して 0 件）

関係する ADR: なし

## 完了の条件

記録か ADR が md5 で引く追跡しない出力が、手順を覚えていなくてもマージの後に残るか戻せる。候補（確かめていない）:

- 引く出力を DVC のステージの出力にする（キャッシュを worktree で共有し、`worktree-restore-outputs.sh` で戻る。ADR-0126、
  ADR-0130）。micro と branches の入力を dvc.yaml にする費用を数える
- `worktree-restore-outputs.sh` の逆向きのスクリプト（ブランチの記録が引く md5 の無視されたファイルを main にコピーする）を
  /git-merge のマージの前に呼ぶ
- `post-merge-cleanup.sh` が、記録が引く無視されたファイルを持つ worktree の削除を止める

## 経緯

- 2026-10-30: /retro（25 回目）で、スキル topic-start の「マージ」の節に、記録が引く md5 と無視されたファイルの md5 を
  シェルで突き合わせて main にコピーする手順を足す案を試して捨てた。ローカルの /code-review 2 周で、cwd が main に戻る
  （コマンドが空振りして何も出さない）、`:(glob)<dir>/*` が下のディレクトリを見ない、`git diff` の削除の行と文脈の行の md5 も
  拾う、DVC の出力を除かない、コピーが main の別の md5 のファイルを確かめる前に上書きする、/git-merge にコピーの段が無く
  後片付けで消える、`md5` を前に置かない md5（「置き換えた組成が da7265c2」）を拾わない、の指摘が出た。手順の文で補うと
  黙って失敗する形が残るので、仕組み（上の候補）で直す
