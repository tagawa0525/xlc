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

`examples/depletion/micro/*.msgpack` などは DVC の出力でなく、`examples/depletion/micro/.gitignore` の `*.msgpack` で無視する生成物で、検証の記録と ADR が md5 で
引く。作ったブランチの worktree にだけあり、/git-merge の後片付け（`post-merge-cleanup.sh` の `git worktree remove`）が
無視されたファイルごと消す。DVC の出力はキャッシュから戻り、`examples/vera/p2/openmc-2K.msgpack` は
`scripts/worktree-restore-outputs.sh` が main からコピーするが、ほかの追跡しない出力には戻す経路が無い。

- 2026-10-30 に、ADR-0221 の依頼が「既存の `shifted.out` と `shifted-2A.msgpack`（md5 456654d5）をそのまま使う」と書いたが、
  どの worktree にも残っていなかった。担当は履歴 `shifted-2A.msgpack` を `make_shifted.py`（Python、数秒）で作り直して md5 が
  一致し、依頼が挙げなかった置き換えた組成のテーブル `table-shifted.msgpack` も残っていなかったので、lattice を Slurm 618
  （20 分 9 秒。新しいテーブルは md5 da7265c2）で解き直した（examples/depletion/verification.md「ミクロ燃焼の線形の補正（2A）」の
  「計算」、担当の報告）。置き換えた組成の入力 `shifted.toml` は ADR-0221 で変わらず、前のテーブルが残っていれば 618 は
  要らなかったかは確かめていない（`bin/lattice` は変わった）
- ADR-0221 のマージ（xlc 9a69db2a9）の後片付けでも、記録と ADR が引く `table-nominal.msgpack`（md5 fde256d9）などの micro の
  テーブルが消え、main にも無い（2026-10-30 に `/home/tagawa/github` の下を `find` で探して 0 件）
- σ の履歴の補正のマージ（xlc fcb8db079。2026-10-30）では、`histories` の形が変わって main の
  `examples/depletion/history/` の生成物が新しい core で読めなくなったので、親が後片付けの前にブランチの worktree の
  `*.msgpack` と `*.log` の 19 件を main に `cp` で上書きした。記録（「履歴の残り R（2A の減速材の密度）」の「計算」）が
  「R は前の版の出力 3d73faee と同じ値」と引く前の形の `table-residual-nominal-histories.msgpack` は、この上書きと
  worktree の削除で残っていない（`/home/tagawa/github` の下を `find` で探して 222fabb0 の 1 件だけ）。`lattice histories`
  で数秒で作り直せるので損は小さいが、形が変わるたびに手でコピーし、前の版の md5 を失う形は同じ

出し直しの経路は、コミットした入力を checkout して再実行する形（ADR-0122 の案 E）で、この issue はその時間を省く置き場を
扱う。追跡しない参照の出力の置き場は #0018。

関係する ADR: なし

## 完了の条件

記録か ADR が引く（md5 か名前で）追跡しない出力が、手順を覚えていなくてもマージの後に残るか戻せる。候補（確かめていない）:

- 引く出力を DVC で追う（キャッシュを worktree で共有し、`worktree-restore-outputs.sh` で戻る。ADR-0126、ADR-0130）。
  DVC のキャッシュは中身の md5 で置くので、形が変わって出力を置き換えても前の版がキャッシュに残る。
  `dvc add` でファイルだけを追うか、既存の examples/depletion/dvc.yaml にステージを足すか（examples/depletion/branches には
  dvc.yaml が無い）の費用を数える
- `worktree-restore-outputs.sh` の逆向きのスクリプト（ブランチの記録が引く md5 の無視されたファイルを main にコピーする）を
  /git-merge のマージの前に呼ぶ
- マージの前に、記録が引く無視されたファイルを持つ worktree の削除を止める（`post-merge-cleanup.sh` はすべてのリポジトリで
  共有するので、xlc の側の検査として置く）

## 経緯

- 2026-10-30: /retro（25 回目）で、スキル topic-start の「マージ」の節に、記録が引く md5 と無視されたファイルの md5 を
  シェルで突き合わせて main にコピーする手順を足す案を試して捨てた。ローカルの /code-review 2 周で、cwd が main に戻る
  （コマンドが空振りして何も出さない）、`:(glob)<dir>/*` が下のディレクトリを見ない、`git diff` の削除の行と文脈の行の md5 も
  拾う、DVC の出力を除かない、コピーが main の別の md5 のファイルを確かめる前に上書きする、/git-merge にコピーの段が無く
  後片付けで消える、`md5` を前に置かない md5（「置き換えた組成が da7265c2」）を拾わない、の指摘が出た。手順の文で補うと
  黙って失敗する形が残るので、仕組み（上の候補）で直す
