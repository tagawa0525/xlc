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
`dvc_history`、`show_dvc`）を足したブランチ `feat/record-provenance-dvc` のレビューの 1 周目と 2 周目（2026-11-06）で、
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

関係する ADR: なし

## 完了の条件

上の点を直すか、直さない理由を書いて閉じる。
