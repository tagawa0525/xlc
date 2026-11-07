---
status: open
labels: [core, lattice, depletion]
opened: 2026-11-07
closed:
milestone: M4
blocked-by: []
resolution:
---

# ノードの出力の重み（#0154）のレビューの残り

## 問題

ブランチ feat/0154-node-power-by-fission-energy のレビューの、直さなかった指摘（再現を書けないもの、重複と効率）:

- crates/core/src/case.rs の `node_constants`: 内挿した Σf が 0 以下の群は核分裂あたりのエネルギーを 0 とし、補正で Σf が正になっても κΣf は 0 の
  まま。補正した κΣf の符号も確かめない（通常の入力で内挿した Σf が 0 以下になる例は見つけていない）
- crates/lattice/src/history.rs の `history_branches`: チェーンと履歴の照合は核種の名前と順だけで、同じ核種の並びで Q の違うチェーン（作り直した
  チェーン）を通す。履歴にチェーンの Q か識別を持たせれば止められる
- crates/lattice/src/history.rs: `kappa_fission` は `micro` と同じ SPH 因子と少数群の中性子束を状態ごとに求め直す（重複と計算の繰り返し）
- crates/core/src/nem.rs: `fission_powers` は `fission_rates` のループのコピー（断面積とエラーだけが違う）
- crates/core/src/rods.rs の `mix`: κΣf の有無を数えてから NaN を返す取り出しで混ぜる（すべてか、どれも持たないかの検査を先に書けば NaN が要らない）
- crates/core/src/depletion.rs: 燃焼の前の κΣf の検査がすべての状態を見る（テーブルの検査が、すべての状態が持つか、どれも持たないかを保証する）
- examples/depletion/compare_core_energy.py: 集合体の燃焼度を重金属の質量から求め直す（履歴の点の `assembly_burnups` と同じ定義）。引数の空の
  パスの検査（`all(str(p) for p in paths)`）は `Path('')` が `.` になるので働かない（compare_core_guess.py と同じ形）
- `kappa_fission` を持たない燃焼の履歴からのテーブル（examples/depletion/branches の `table-2A-*` のうち `table-2A-micro-boron-rod-axis.msgpack` を除くもの、
  `table-2O-all.msgpack`、examples/depletion/history の `table-*`）は、core の燃焼で使うとエラーで止まる。それらを読む記録を再実行するときに作り直す
  （CLAUDE.md の後回しの規則）。examples/depletion/p9-timing の `core-p9` のステージは、作り直した `table-2A-micro-boron-rod-axis.msgpack` と `bin/core`
  が deps で changed のまま
- examples/verification_tables.py の `core_energy_points`: 核分裂率の燃焼の出力の割合が 0 の集合体を、κΣf の燃焼の割合によらず表から落とす
  （食い違いを TableError にしない）
- 「restart と装荷替え（小さな炉心）」の記録が md5 で引く core の履歴（`branches/core-small-restart-full.msgpack` ほか）は `dvc add` で追跡して
  いないので、同じケースを解き直すと上書きされて戻せない（Slurm 721 が #0154 のブランチの worktree で上書きした）
- examples/depletion/verification.md の「ノードの出力の重み（小さな炉心）」（例外の周の 2 周目の指摘。直していない）:
  - 最初の `scripts/dvc-checkout.py` と `C=` の行が `&&` の連鎖の外にある
  - lattice の手順が dvc add で追うテーブルをその場で上書きし、作り直したテーブルが記録の md5 と同じかを確かめない
  - コマンドは `trials` に比較の JSON を書くので、記録が引く `branches/compare-core-small-energy.json`（入力のパスが `branches/energy`）を作り直す手順が無い
  - 「計算」の Slurm 723 の項は、汚れた作業ツリーの変更を記録の 1 ファイルとしか書いていない（ほかに docs/issues/0266 と docs/plans/m4-depletion.md も
    コミットしていなかった）。Slurm 722 の項は、テーブルの値の比較も流したのに「燃焼から比較まで」と書く。723 は出力のそろった worktree で走り、
    DVC の出力を戻す経路を通っていないのに、流していない部分を build と lattice だけとしている。どちらも、md5 と `cases` を比べたのは 2 回目の出力だけ
- docs/issues/0266 は、直している最中の記録の節の特定の版（0edff5d2e）を案内する（節は今のファイルにあり、測ったコミットと日付で引ける）
- docs/plans/m4-depletion.md は、#0154 の完了を項目 7 の済んだ部分でなく、番号のない段落に書き、マージのコミットを持たない

関係する ADR: なし

## 完了の条件

点ごとに、直すか直さないかを決めて、直すものを直す。

## 経緯

- 2026-11-07: #0154 のブランチの 1〜3 周目と記録の周の後の例外の周のレビューで出た、直さなかった指摘を集めた
