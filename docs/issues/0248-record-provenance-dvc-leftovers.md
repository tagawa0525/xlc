---
status: open
labels: [tooling, dvc]
opened: 2026-11-06
closed:
milestone:
blocked-by: []
resolution:
---

# record_provenance.py の .dvc の読み方の残り

## 問題

`examples/record_provenance.py`（ADR-0176 と ADR-0222 の帰結のツール）に `.dvc` の md5 の履歴を出す機能（`parse_dvc`、
`dvc_history`、`show_dvc`）を足したブランチ `feat/record-provenance-dvc` のレビューの 1〜3 周目と記録の周（2026-11-06）で、
直さずに残した指摘:

- `dvc_history` は `.dvc` を触ったコミットの版をすべて `git show` して `parse_dvc` で読んでから履歴をさかのぼる。`.dvc` を消した
  コミットより前に読めない版があると、さかのぼりが止まるはずの所より前の版で全体が `FAILED` になり、版の数だけ `git show` を
  実行する
- `show_dvc` は `uncommitted` と `unknown` でも終了コード 0 を返す（dvc.lock のステージの `locked_since` も同じ）
- 履歴は `git log -- <パス>` の新しい順に並べて同じ md5 が続く間を 1 行にまとめるので、マージの両側が `.dvc` を変えると、
  2 本のブランチのコミットが日付の順に混じり、md5 を入れたコミットを取り違えうる（`locked_since` も同じ）
- `dvc_history` と `locked_since` が、`git log` と `git show <コミット>:./<名前>` で版をさかのぼる同じ処理を別々に持つ
- `uncommitted` のとき、作業ツリーの md5 の行の下にコミット済みの履歴を出す（`locked_since` は `uncommitted` のときコミットを
  出さない）
- md5 を先頭 8 桁で出すので、ディレクトリの出力の `.dir` の印が出ない（dvc.lock のステージの出力も同じ）
- （2 周目）`git log -- <名前>` と `git show <コミット>:./<名前>` は名前を pathspec として読むので、`[`、`*`、`?` を含む名前は
  別のファイルに当たる（`--literal-pathspecs` が候補。xlc の生成物の名前には無い）
- （2 周目）スキル verification-record の「Y の行のコミットが X か X の子孫」は下限だけを確かめ、記録が古すぎる X を書いても通る。
  同じ md5 が何行かに出るときに、どの行と照らすかも書いていない
- （2 周目）`parse_dvc` と `parse_stage` は手で書いた YAML の読み手で、dvc の書き方（折り返し、引用符）ごとに扱いが要る。
  PyYAML（default の devShell にある）の `yaml.safe_load` で読む案は、ADR-0176 の「標準ライブラリで済むこと」と合わないので、
  採るなら ADR-0176 を見直す
- （3 周目）`path:` の次の行に引用符で囲んだパスがあると、`parse_dvc` は引用符ごと出す（dvc 3.67.1 は引用符の名前を同じ行に
  書く。手で直した .dvc か別の版で起きうる）。空の `path:` の後に何行も続く形と、この形にはテストが無い
- （3 周目）ディレクトリ（`.dvc` の設定のディレクトリなど）を渡すと `IsADirectoryError` の traceback で止まり、残りの引数を読まない
- （3 周目）`outs:` の中の 1 桁目の注釈の行で `parse_dvc` が読むのをやめる（dvc は注釈を書かない）
- （3 周目）テストの `contextlib` と `io` の import が関数の中に 3 か所ある
- （記録の周）スキル verification-record の .dvc の項は、d03a12687 で足した「md5 が前の版と同じなら行は X より前のコミットを
  指す」と、同じ段落の「作り直した生成物は md5 が版ごとに違う」「出るコミットは…生成物を作ったコミットの後にある」が食い違う。
  record_provenance.py の docstring の「生成物を作ったコミット（その前）ではない」も同じ。段落の 1 行が折り返されていない

関係する ADR: なし

## 完了の条件

上の点を直すか、直さない理由を書いて閉じる。

- （例外の周、親の周）`parse_dvc` は `-` だけの行（キーを次の行から書く項目）を区切りと認めず、続く行を前の項目に上書きして out を黙って落とす（dvc は書かない形）。`path:` の値の後ろの注釈（`# …`）を剥がさない。`state` の判定が (path, md5) だけで docstring の「outs が HEAD と違えば uncommitted」より狭い。`.dvc` でない引数を混ぜたときと `show_dvc` の「does not exist」のテストが無い。`\r\n` の `.dvc` を「no outs」と報告して原因を示さない
- （例外の周の 2 周目。blocking。806ac120a で直した）スキル verification-record の .dvc の項で、md5 が前の版と同じときの文（4233c858f「記録には Y を入れた行のコミット（X より前）と、X で作り直した事実の出どころ（X のジョブのログ）の両方を書く」）に、再現できる誤りが残る: ジョブのログにはコミットが入らない（~/.claude/scripts/slurm-run.sh。nixfiles の PR #627 で `commit=` をジョブのコメントに足す）、「Y の行のコミットが X か X の子孫」の検査はこの場合に必ず失敗するのに分岐を書いていない、同じ節の 4（計算の欄には実行したコミット X）と違う第 3 のコミットを記録に足す。前の 2 回の直し（351b8ead1、4233c858f）も同じ文で (a) が出たので、/git-merge の「直した指摘が次の周にも出たとき」に当たり、親が止めた。候補: 文を「md5 が同じなら行は X より前のコミットを指す（中身は同じ）」の事実だけにし、X の出どころは同じ節の 4 の X と Slurm の番号に任せる（記録の今の書き方。examples/depletion/verification.md の 151〜153 行）。修正方法は親に任された（ユーザー 2026-11-07 の判断の要約）ので、806ac120a で候補のとおり文を事実だけにした
- （例外の周の 2 周目）`git()` が終了コードが 0 でないときに理由を問わず None を返し、`dvc_history` が `git show` のどんな失敗も「.dvc を消したコミット」と同じに扱って履歴を黙って打ち切る（index.lock、リポジトリの外のパス）。例外の周の指摘の箇条書きが「完了の条件」の下に付いていて「問題」の節に無い
- （例外の周の 3 周目）スキル verification-record の .dvc の項は、806ac120a で「この出力だけでは、X で作り直したかは決められない」を消したので、md5 が同じときの行（X より前のコミット）と、X で作り直していない記録の行が見分けられないことを書いていない
