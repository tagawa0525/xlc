---
status: open
labels: [vera, openmc]
opened: 2026-10-14
closed:
milestone: M3
blocked-by: []
resolution:
---

# 1B の OpenMC の参照解を 2026-10-12 のライブラリと標準のタリーで作り直す

## 問題

examples/vera/dvc.yaml の `openmc-1B` と `openmc-1B-extract` は、2026-10-11 の 2ff665e4（565 K のライブラリで relock）の
deps で lock されたままで、`dvc status` で changed になる。そのあとの 2026-10-12 の変更が deps を変えた:

- 319a7948: OpenMC のライブラリ（`openmc-library`）に Ti と 565 K の構造材の核種を足した（#0073）
- 06434194: 参照解に標準のタリーの組（ADR-0152）を仕込んだ（`openmc_model.py`、`openmc_tallies.py`）
- e905948a、45e8c77f、ac3fc088、7340855b: run の補助（`examples/openmc_run.py`）と model の整理

`dvc repro -P` は、この 2 段階を走らせ、1B の参照解を今のライブラリとスクリプトで作り直す（16 スレッドで 40〜60 分）。
高速化の後の記録の測り直し（ブランチ docs/remeasure-after-speedups）は、参照解を替えないためにこの 2 段階を走らせず、
古い lock のまま残した。下流の `compare-1B` も、上流が changed なので走らせていない（xlc の側の `lattice-1B` は走らせた）。

標準のタリーの組の決定の範囲でなく、その帰結の走らせ直しで、関係する issue は #0073。

関係する ADR: なし

## 完了の条件

- `dvc repro examples/vera/dvc.yaml:compare-1B` で `openmc-1B`、`openmc-1B-extract`、`compare-1B` を走らせ、`dvc.lock` を
  コミットする（Slurm で排他、スキル openmc-reference）
- examples/vera/verification.md の「参照解」の 1B の行と、「OpenMC との比較」の 1B の差を、新しい参照解で書き直す。
  2026-10-11 の参照解（1.181783 ± 0.000063）との差が統計誤差の中かを書く

## 経緯

- 2026-10-14: 高速化の後の測り直しで `dvc status` が 2 段階を changed と出し、分けた
