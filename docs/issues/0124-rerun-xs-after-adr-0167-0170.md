---
status: closed
labels: [verification]
opened: 2026-10-14
closed: 2026-11-19
milestone:
blocked-by: []
resolution: "wontfix: 前提が変わった。ADR-0167 と ADR-0170 の後のステージは、その後のライブラリの作り直しで何度も再実行した。残る記録の本文は #0224、OpenMC の 1B は #0120 が扱う"
---

# ADR-0167 と ADR-0170 の後に xs のステージとその下流を再実行する

## 問題

高速化の後の記録の測り直し（ブランチ docs/remeasure-after-speedups、Slurm のジョブ 332）は、分岐点 c59097f0 のコードで
走った。そのあと main に ADR-0167（MLBW のエネルギーに依存する散乱半径）と ADR-0170（J・K の表の θ の下限）が入り、
`bin/xs` が変わった。このため、マージの後の `dvc status` で次が changed になる（`bin/lattice` と `bin/core` は分岐点と同じ md5）:

- examples/vera/dvc.yaml の `xs-range`、`xs-nuclide@*`（60 ステージ）、`xs`（ジョブ 332 で合わせて約 26 分）
- examples/core、m0、reflector の `xs`
- examples/depletion の `chain`（`bin/chain` が main で変わった）

ライブラリが変われば、その下流（`lattice-1B`、examples/vera/p3、p5 のすべてのステージ、core、m0、reflector の格子と炉心）も走る。

ステージでない記録の数値（examples/vera/verification.md の「OpenMC との比較」の 1C〜2Q、examples/mox の k と実効断面積）も、
同じ分岐点のライブラリのもの。

ADR-0167 は NRO=1 の MLBW の評価だけに効き、VERA の核種はそれまでの xs で処理できていたので、ADR-0167 では
ライブラリは変わらない見込み（確かめていない）。ADR-0170 の影響は見積もっていない。

関係する ADR: なし

## 完了の条件

- `dvc repro -P` を Slurm で排他で流し、`dvc.lock` をコミットする。このコマンドは #0120 の OpenMC のステージも実行するので、
  #0120 と同じジョブで閉じる
- ライブラリと出力が変わったかを記録ごとに書き、変わった記録は数値と差の表を書き直す。ステージでない 1C〜2Q と MOX も、
  ライブラリが変わったら再実行する

## 経緯

- 2026-10-14: 再測定のブランチの rebase で、xs のステージが古くなったことが分かり、#0120 から分けた
- 2026-11-19: issue の棚卸しで閉じた。xs とその下流は、その後のライブラリの作り直し（ADR-0177 など）と新しい既定の後の再実行で何度も再実行した。表の外の本文（vera の「OpenMC との比較」と mox を含む）は #0224 が、OpenMC の 1B の参照解は #0120 が受け持つ
