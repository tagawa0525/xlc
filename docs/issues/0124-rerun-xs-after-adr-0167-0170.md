---
status: open
labels: [verification]
opened: 2026-10-14
closed:
milestone:
blocked-by: []
resolution:
---

# ADR-0167 と ADR-0170 の後に xs の段階とその下流を再実行する

## 問題

高速化の後の記録の測り直し（ブランチ docs/remeasure-after-speedups、Slurm のジョブ 332）は、分岐点 4b644566 のコードで
走った。そのあと main に ADR-0167（MLBW のエネルギーに依存する散乱半径）と ADR-0170（J・K の表の θ の下限）が入り、
`bin/xs` が変わった。このため、マージの後の `dvc status` で次が changed になる（`bin/lattice` と `bin/core` は分岐点と同じ md5）:

- examples/vera/dvc.yaml の `xs-range`、`xs-nuclide@*`（60 段階）、`xs`（ジョブ 332 で合わせて約 26 分）
- examples/core、m0、reflector の `xs`
- examples/depletion の `chain`（`bin/chain` が main で変わった）

ライブラリが変われば、その下流（`lattice-1B`、examples/vera/p3、p5 のすべての段階、core、m0、reflector の格子と炉心）も走る。

段階でない記録の数値（examples/vera/verification.md の「OpenMC との比較」の 1C〜2Q、examples/mox の k と実効断面積）も、
同じ分岐点のライブラリのもの。examples/m1 の記録（e4dbc0de）は、このブランチで再測定していない。

ADR-0167 は NRO=1 の MLBW の評価だけに効き、VERA の核種はそれまでの xs で処理できていたので、ADR-0167 では
ライブラリは変わらない見込み（確かめていない）。ADR-0170 の影響は見積もっていない。

関係する ADR: なし

## 完了の条件

- `dvc repro -P` を Slurm で排他で流し、`dvc.lock` をコミットする。このコマンドは #0120 の OpenMC の段階も実行するので、
  #0120 と同じジョブで閉じる
- ライブラリと出力が変わったかを記録ごとに書き、変わった記録は数値と差の表を書き直す。段階でない 1C〜2Q と MOX も、
  ライブラリが変わったら再実行する。M1 も今のコードで再測定する

## 経緯

- 2026-10-14: 再測定のブランチの rebase で、xs の段階が古くなったことが分かり、#0120 から分けた
