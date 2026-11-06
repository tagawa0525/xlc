---
status: open
labels: [core, maintenance]
opened: 2026-11-06
closed:
milestone:
blocked-by: []
resolution:
---

# core の CMFD のレビューで直さなかった指摘

## 問題

ADR-0232 の実装（ブランチ perf/0076-core-cmfd）のローカルのレビューで、/git-merge の「直す指摘の基準」の blocking に当たらず直さなかった指摘:

1 周目（c062017c、2026-11-06）:

- CMFD は既定で、補正を作るときにノードの中性子束が正でなければ止まる。加速しない反復は途中の負の値を通し、収束した時点だけを確かめていた。
  止まるのは crates/core のテストの強い吸収体（Σt h = 1000）のノードで、ほかの例（小さな炉心、全炉心の BOC）では起きていない。負のノードを
  補正から外すか、そのノードだけ比を 1 にする形などは測っていない（#0242 と同じ枝）
- `stepped_core_converges_with_refinement` の吸収体で埋めた比較は加速しない反復で解く（CMFD のこの問題でのエラーは `cmfd_needs_positive_node_fluxes`
  が確かめる）。上の扱いを直したら既定に戻す
- 変化の見積もり `Remaining` が `crates/core/src/nem.rs`、`crates/cmfd/src/lib.rs`、`crates/lattice/src/moc.rs` に 3 つある。`xlc-cmfd` の
  ものを公開して core と lattice が使う形（#0241 で lattice を移すときにまとめる）
- `xlc_cmfd::Coarse::new` は `labels`、`matrices`、`fission`、`neighbors`、`flux` の長さを確かめない（今の呼び出しは core だけで長さは揃う）
- `Core::coarse` は外側の反復ごとにノードごとの `Vec` を作り、`Coarse::new` が平らな配列にコピーする。D̃ と ADF の平均は物質を置き換えるまで
  変わらないが、毎回求め直す（全炉心の BOC で 1 回 12 ms、外側の反復 1 回のスイープの約半分）
- `Spectrum::chi_sum` は χ の和で、呼び出し側が求める。core の呼び出しはセルに 1 つの項だけを渡す（lattice の移行に備えた形）
- `small_core_cases.py` の `acceleration` の行は `radial-nodes = 1` の行の後に置き換えで足し、その行が無ければ黙って足さない（テスト
  `test_unaccelerated_cases_differ_only_by_the_acceleration` は止める）
- ADR-0232 の検討した案の回数と時間は、コミットしていない試行（環境変数で定数と形を切り替えた版）の値（決めた形の値は「確認」の記録の値）
- `compare_core_acceleration.py` の `per_solve` は、`compare_core_guess.run_of` が読んだのと同じ span をもう一度読む

2 周目（339e17ce、2026-11-06）:

- `docs/plans/architecture.md` の共通のクレートの一覧に `xlc-cmfd` が無い（依頼の変える範囲の外）
- `verification_tables.py` の `core_acceleration_runs` は 1 回の時間を小数 2 桁で描き、小さな炉心の最小（0.003 s など）が 0.00 になる（比較のツールの
  `--report` は 3 桁）
- `compare_core_acceleration.py` の `main`、`parse`、`report` の前半と、`core_acceleration_runs` は、初期値と Xe とホウ素の比較のものの 3 つ目の写し
- `xlc-numeric` にモジュールを足すと、それを使わない xs と lattice の実行ファイルのバイトが変わった原因を確かめていない（ADR-0232 の案 S2）。
  共通のクレートに足すたびにステージが changed になるなら、構造の側で直す候補がある

3 周目（416d97a9、2026-11-06）:

- examples/depletion/verification.md の初期値の節に足した段落の「追跡する。。」の句点の重なりと、段落の末尾の句点の抜け
- 同じ段落の「この節のログは追っていない」が探索の履歴にだけ掛かって読める（平衡 Xe のログ md5 02a10843 も追っていない。nested-1 の `.log.dvc` は別の実行）
- 初期値の節の平らな初期値の履歴（`compare-core-small-guess.json` の md5 058973f1 と d6762f89）は `.dvc` が無い（このブランチより前から）
- 平衡 Xe の燃焼の差の原因の候補に、固有値計算の数の違い（267 回と 255 回）を並べていない
- ADR-0232 の `issues` から #0241〜#0243 を外したが、docs/adr/README.md は「帰結」に書いた残りの issue を `issues` に足してもよいとする（どちらも規則の中）

関係する ADR: なし

## 完了の条件

それぞれ直すか、直さない理由を書いて閉じる。

## 経緯

- 2026-11-06: ADR-0232 のブランチの 1 周目のレビューから
- 2026-11-06: 2 周目のレビューから
- 2026-11-06: 3 周目のレビューから
