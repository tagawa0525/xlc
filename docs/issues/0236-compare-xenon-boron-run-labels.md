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

比較のツール（compare_core_pair.py、compare_core_threads.py）は、2 本の燃焼の役を、今ディスクにある core の入力の項目
（`initial-guess`、`xenon-boron`、`acceleration`、`node-order`）の値で決める。履歴とログはこれらの項目を書かないので、燃焼の後に
入力を書き換えたときや、入力を変えずに別の版の bin/core で解いた出力が残っているときは、実行に使った項目と違う役で比べうる。例:
`branches/core-small-search-1.msgpack` は今 CMFD の燃焼で、記録の「固有値計算の初期値（小さな炉心）」の平らな初期値の燃焼は CMFD を
足す前の bin/core で解いた。どちらの入力も `acceleration` を書かないので、平らな初期値の古い出力が残っていれば、記録の `initial-guess`
のコマンドは止まらずに違う解き方の 2 本を比べる（その出力は手元に無く、実行しては確かめていない）。

重複の共通化と、入力による役の確認は済んだ（経緯の 2026-11-19）。それより前の版の本文は main の d17396374 にある。

関係する ADR: なし

## 完了の条件

core が、比べる項目（または入力の md5）をログか履歴に書き、比較のツールがそれと入力を突き合わせて止まる（crates/core の変更）。

## 経緯

- 2026-11-06: perf/0156-xenon-boron-joint のレビューの 1 周目で挙がった（コピーの重複と、履歴が `xenon_boron` を持たないので引数の順を確かめられない
  こと）。再現する誤りではないので、このブランチでは直さない
- 2026-11-06: 親の周（1 周目、ブランチ全体）で挙がった 4 件（小さな炉心のログを `dvc add` していない、nested の履歴の md5 が ADR-0231 の記録の履歴と同じ、
  ADR-0219 の「確認」のテストの名前、verification_tables.py の GIT_TRACKED の注）は、親が 2 周目の後にマージの前に直した
- 2026-11-19: 3 つのツールを compare_core_pair.py（項目を最初の引数で受け、要約の `format` と `runs` の鍵を表に持つ）にまとめ、
  verification_tables.py の 3 つの表も `core_pair_runs` にまとめた。`rho_pcm`、`BURNUP_TOLERANCE`、`same_burnup` は core の比較のツールの分を
  core_history.py にまとめた（history/compare_core_history.py と compare_core_energy.py の分も）。compare_depletion.py と compare_depletion_assembly.py は格子計算と OpenMC の燃焼を比べるツールで、core の
  履歴を読まないので、自分の許容のまま残した。取り違えは、ツールが燃焼ごとに core の入力を受け、入力が名前を書く履歴と入力の名前の `.log` を読み、2 つの入力が項目（と `history`、`threads`）のほかで違えば止まり、
  項目（`initial-guess`、`xenon-boron`、`acceleration`、`node-order`）の値で 2 本の役を確かめる形で止めた（compare_core_threads.py は入力の
  `threads` とログのスレッドの数も比べる）。同じ入力を 2 つの版の bin/core で解いた燃焼は入力で見分けられないので、比べる形をやめた
  （compare_core_threads.py の ADR-0234 の前の bin/core の扱いを含む）。その形で書いた要約（小さな炉心の #0244 の前の CMFD、Problem 9 の
  形の試行）は、書いたツールの名前を記録に書いて残した
