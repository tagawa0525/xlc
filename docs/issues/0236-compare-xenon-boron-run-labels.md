---
type: Issue
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
共有する）。compare_core_acceleration.py が 3 つめのコピー（下）で、CLAUDE.md の「重複は 3 回まで許容」を超えたので、まとめる時に当たる。

燃焼の履歴は `xenon-boron` と `initial-guess` を書かない（ADR-0231 と同じく、どちらも同じ収束の判定を満たす解で履歴の形を変えない）ので、比較の
ツールは引数の順だけで 2 本を見分ける。joint と nested を入れ替えて渡すと、`SHARED` の項目の確認を通り、差の符号と固有値計算の比が逆の要約になる。

examples/depletion の比較のツールの `rho_pcm`（compare_core_cycle.py、compare_core_micro.py、compare_core_partial_rods.py、compare_core_b10.py、
compare_core_restart.py）と、燃焼度の一致の許容 `BURNUP_TOLERANCE`（compare_core_depletion.py、compare_depletion.py、compare_depletion_assembly.py、
compare_core_micro.py、compare_core_partial_rods.py）と `same_burnup` は、3 つを超える複製で、examples/depletion/core_history.py にまとめられる。
compare_core_acceleration.py の `main`、`parse`、`report` の前半と verification_tables.py の `core_acceleration_runs` は、compare_core_guess.py と
compare_core_xenon_boron.py に続く 3 つ目のコピーで、まとめる時に当たる。compare_core_threads.py も燃焼の順（直方体の並びの順か色の順か）を引数の
位置で決め、ログで確かめない（core のログは `node-order` を出さない）ので、上と同じく取り違えを見分けられない。

関係する ADR: なし

## 完了の条件

1 項目を変えた 2 本の core の燃焼を比べるツール（compare_core_guess.py、compare_core_xenon_boron.py、compare_core_acceleration.py）を 1 つにする
（変えた項目の名前、2 本の名前、`format` を引数か表で渡す）。`rho_pcm`、`BURNUP_TOLERANCE`、`same_burnup` を core_history.py にまとめる。取り違えは、ログに入力の項目を出す（`depletion` の始めの行など）か、履歴に書くかで見分けられるようにするかを決める。

## 経緯

- 2026-11-06: perf/0156-xenon-boron-joint のレビューの 1 周目で挙がった（コピーの重複と、履歴が `xenon_boron` を持たないので引数の順を確かめられない
  こと）。再現する誤りではないので、このブランチでは直さない
- 2026-11-06: 親の周（1 周目、ブランチ全体）で挙がった 4 件（小さな炉心のログを `dvc add` していない、nested の履歴の md5 が ADR-0231 の記録の履歴と同じ、
  ADR-0219 の「確認」のテストの名前、verification_tables.py の GIT_TRACKED の注）は、親が 2 周目の後にマージの前に直した
