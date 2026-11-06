---
status: open
labels: [depletion, tooling]
opened: 2026-11-06
closed:
milestone:
blocked-by: []
resolution:
---

# core の燃焼の 1 項目を変えた比較のツールの重複と、比べる 2 本の取り違え

## 問題

examples/depletion/compare_core_xenon_boron.py（ADR-0219 の `xenon-boron` の joint と nested を比べる）は、compare_core_guess.py（ADR-0231 の
`initial-guess`）の `compare` を呼ぶが、`report`、`parse`、`main` は名前と `format` と run の鍵のほかは同じコピーで、ログの読み方や md5 の記録を
変えるときは 2 か所を直す。examples/verification_tables.py の `core_xenon_boron_runs` も `core_guess_runs` のコピー（差の表は `core_pair_differences` を
共有する）。コピーは 2 つめで、CLAUDE.md の「重複は 3 回まで許容」の中。

燃焼の履歴は `xenon-boron` と `initial-guess` を書かない（ADR-0231 と同じく、どちらも同じ収束の判定を満たす解で履歴の形を変えない）ので、比較の
ツールは引数の順だけで 2 本を見分ける。joint と nested を入れ替えて渡すと、`SHARED` の項目の確認を通り、差の符号と固有値計算の比が逆の要約になる。

関係する ADR: なし

## 完了の条件

1 項目を変えた 2 本の core の燃焼を比べるツールを 1 つにする（変えた項目の名前、2 本の名前、`format` を引数か表で渡す）か、3 つめのコピーが
要るときにまとめる。取り違えは、ログに入力の項目を出す（`depletion` の始めの行など）か、履歴に書くかで見分けられるようにするかを決める。

## 経緯

- 2026-11-06: perf/0156-xenon-boron-joint のレビューの 1 周目で挙がった（コピーの重複と、履歴が `xenon_boron` を持たないので引数の順を確かめられない
  こと）。再現する誤りではないので、このブランチでは直さない
- 2026-11-06: 親の周（1 周目、ブランチ全体）で残した指摘: (1) 記録と compare-core-small-xenon-boron.json が md5 で引く core-small-equilibrium-1.log と
  core-small-equilibrium-nested-1.log は branches/.gitignore の `*.log` で無視され、`dvc add` も無視されたファイルを受け付けないので追っていない（compare-core-small-guess.json の
  ログと同じ扱い。ログの扱いは、ADR-0222 の対象を比較の JSON が引く生成物に広げるブランチ docs/0222-compare-json-generated の issue の対象）。(2) nested の履歴 core-small-equilibrium-nested-1.msgpack（md5 0c3828ae）は、compare-core-small-guess.json が
  ADR-0231 の core-small-equilibrium-1.msgpack に記録した md5 と同じで（同時の反復の前の形を bit 単位で再現する）、joint の実行がその名前の履歴を md5 7620f11e に
  上書きしたことを記録が書いていない。(3) ADR-0219 の「確認」が、決定の「端で解く」「効きが分かるまで探す」「有限でない濃度はエラー」のテストの名前を挙げていない。
  (4) examples/verification_tables.py の GIT_TRACKED の注が trial-0231 の理由（#0235）だけで、trial-0156 の理由（#0156 の試行）を書いていない
