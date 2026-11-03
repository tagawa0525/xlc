---
status: open
labels: [tooling, dvc]
opened: 2026-11-03
closed:
milestone:
blocked-by: []
resolution:
---

# dvc-relock のマージの直前の段を、手順の入口と許可に通す

## 問題

ブランチがステージを消すか outs の名前を変えると、main の作業ツリーに残る古い出力がマージの結果で未追跡になり、
pre-merge-commit の hook に止まる。dvc-relock の「マージが止まったら」は、nixfiles #569（2026-11-03）の hook が分けて出す
「マージするブランチが .gitignore から外した生成物」の一覧を受け、main の HEAD の記録（`dvc status` の `changed outs`、
キャッシュ）で消してよいかを決める。前の版の /code-review（2026-11-03 の 1〜3 周目と記録の周、HEAD 9a911d605）が
挙げた点のうち、次は節の作り直しと hook で解いた:

- `modified` の出力も消す対象に読めた（消すと戻らない）: 節が `dvc status` の `modified:`、`not in cache:`、記録が無い、を
  消さないものとして分けた（`Output.status` は dvc 3.67.1 の dvc/output.py）
- hook のパスを `.gitignore` の消えた行と手で照合し、ディレクトリの outs で迷う: hook が HEAD の `.gitignore` で無視されて
  いたファイルを別の一覧に出す（ディレクトリのパターンの中のファイルも入る。`git check-ignore --stdin` で試した）。hook が別の
  一覧を出せないとき（判定の失敗、古い hook）は、`git merge --abort` の後の `git check-ignore` で同じ判定をする
  （`git status --short --ignored` はディレクトリの中のファイルを `!! <ディレクトリ>/` とまとめて出すので使わない）
- hook の一覧に前からある未追跡のファイルや変更が混じる: hook が生成物を別の一覧に分ける
- 記録の無いファイルを消しうる（例の `dvc.yaml` 全体の `dvc status` は、出力でないパスに何も出さない）: パスごとの
  `dvc status <パス>` にし、出力でなければ失敗することで分けた（examples/core で試した）。戻すときもパスだけを指定する
- `dvc status` と /topic-start の「main の作業ツリーでは dvc を実行しない」、/git-merge の「hook の拒否はブランチで直す」との
  食い違い: 節が、ブランチでは直せない main の作業ツリーの追跡しないファイルを扱う例外と書いた
- スキルの入口の文と description が段を含まず、節が前半の `###` だった: 前半と後半の間の `##` にし、入口と description に足した

残る点（直していない。手順のとおりに実行して失敗する誤りではない）:

- /topic-start（委任を受けた担当も読む）は dvc-relock の後半だけを挙げ、マージが止まったときの段に触れない（CLAUDE.md の
  作業の流れからはたどれる）
- dvc-relock の allowed-tools に `dvc status examples/*` は足したが、main の worktree での `git merge --abort` と `rm` は
  無い（`rm` を許可に入れると、ほかの場面の削除も確認なしに通る）
- dvc.yaml を手で直して `.gitignore` の行が残ったとき、持ち主の無い出力は無視されたまま main に残る（`.gitignore` が変わらない
  ので hook は止まらない）
- `cache: false` に変えて git で追跡するようにした出力は、main の無視された版がブランチの追跡する版と違っても、git がマージで黙って
  上書きする（マージの結果と作業ツリーが同じになるので hook は止まらない）。前の版の節は前もって `.gitignore` の差を読んで
  これを見ていたが、hook の一覧を受ける形にしたので見なくなった（作り直した版のレビューの指摘）
- `dvc status` の状態の読み分けを手で行う。パスを受けて消してよいかを出すスクリプトにする案がある（作り直した版のレビューの
  指摘。試していない）。パスごとの `dvc status` は、ディレクトリの出力の中のファイルが多いと同じ答えを何度も求める
- `cache: false` の上書きは、/git-merge のマージを `git merge --no-overwrite-ignore`（無視されたファイルを上書きせずに止まる）で
  行えば止められる候補がある（2 周目のレビューの指摘。確かめていない。/git-merge は nixfiles のスキル）
- main の作業ツリーの古い出力を消す判断と、/topic-start と /git-merge の規則の例外は、スキルの節にだけあり ADR に無い
  （2 周目のレビューの指摘）
- 残すと決めた出力（`modified:` でユーザーが残すものなど）について、hook の示す「ブランチの .gitignore に行を戻す」を節が
  扱わない。allowed-tools に `git check-ignore`、`examples/` の外のパス（`bin/` など）の `dvc status`、`cd <main の worktree> && …`
  の形も無い（3 周目のレビューの指摘）

記録の周（2026-11-03、HEAD f29077a67）が挙げた blocking 2 件（`md5sum` の方法に `cache: false` の除外が無い、/git-merge の
扱いの誤り）は、親の判断で直した（`md5sum` の方法は外した）。残る点: /git-merge（nixfiles のスキル）は作業ツリーの変更による
拒否に触れず、この節への入口が無い。

捨てた方式: 前もって `.gitignore` の差を手で読んで消す方式（2 周目までの版）は、相対パス、`.dvc` ファイル、直下のパターン、
消して足し直した行で漏れ、`cache: false` の出力は消す要が無いのに消して戻せないことがある（2 周目と 3 周目の指摘）ので捨てた。
hook の一覧を `.gitignore` の差と手で照合し、`dvc status` で戻せるかを見る方式（4710ed575 までの版）は、記録の周の上の
blocking（`modified` を消しうる、照合で迷う、規則との食い違い）で捨て、hook の印と `dvc status` の状態ごとの扱いに替えた。
main とブランチの dvc.lock の outs を比べる案は、`.dvc` ファイルとディレクトリの outs を別に扱う手間と、hook の一覧で足りる
かの見込みから試していない（確かめていない）。

関係する ADR: なし

## 完了の条件

- /topic-start からマージが止まったときの段にたどれ、その段の操作（`git merge --abort`、消すこと）の許可を決める
- 手で直した dvc.yaml、`cache: false` の上書き、読み分けのスクリプト、ADR に書くかを、直すか、直さない理由を書いて閉じる

## 経緯

- 2026-11-03: #0093 のマージ（dcc9036a2）が、examples/core の持ち主の無くなった古い出力 6 本で pre-merge の hook に止まった。
  /retro の docs/dvc-relock-removed-stage-outs で dvc-relock の節を hook の一覧を使う形に直し、残りをここに書いた
- 2026-11-03: nixfiles #569 で hook が `.gitignore` から外れた生成物を分けて示すようになったので、節を hook の一覧を受けて
  DVC の記録で決める形に作り直し、解いた点と残る点を分けた
