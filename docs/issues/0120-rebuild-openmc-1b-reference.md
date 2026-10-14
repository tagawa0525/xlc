---
status: open
labels: [vera, openmc]
opened: 2026-10-14
closed:
milestone: M3
blocked-by: ["#0124"]
resolution:
---

# 1B の OpenMC の参照解を 2026-10-12 のライブラリと標準のタリーで作り直す

## 問題

examples/vera/dvc.yaml の `openmc-1B` と `openmc-1B-extract` は、2026-10-11 の 9c121906（565 K のライブラリで relock）の
deps で lock されたままで、`dvc status` で changed になる。そのあとの 2026-10-12 の変更が deps を変え、ADR-0136 の段階の
分割で `openmc-1B` の cmd と deps（`openmc-1B-model` の model.xml、ライブラリ、入力の TOML）も変わった:

- 8ffd3162: OpenMC のライブラリ（`openmc-library`）に Ti と 565 K の構造材の核種を足した（#0073）
- 3343ca05: 参照解に標準のタリーの組（ADR-0152）を仕込んだ（`openmc_model.py`、`openmc_tallies.py`）
- f86b2f2a、5b90a5d7、1fd290c3、2982b580: run の補助（`examples/openmc_run.py`）と model の整理

`dvc repro -P` は、この 2 段階を実行し、1B の参照解を今のライブラリとスクリプトで作り直す（`openmc-1B` は 16 スレッドで約 18 分。
ADR-0136）。ADR-0136 は OpenMC の段階を `openmc-1B-model`（model.xml）と `openmc-1B` に分け、標準のタリーで model.xml が
変わったため実行の段階は古いまま残し、作り直しは別に行うとした（#0018 の残りの項目にも同じ記述がある）。この issue はその作り直しを受け持つ。
高速化の後の記録の測り直し（ブランチ docs/remeasure-after-speedups）は、参照解を替えないためにこの 2 段階を実行せず、
古い lock のまま残した。下流の `compare-1B` も、上流が changed なので実行していない（xlc の側の `lattice-1B` は実行した）。

標準のタリーの組の決定の範囲でなく、その帰結の再実行で、関係する issue は #0073。

関係する ADR: なし

xs の段階とその下流（`lattice-1B` を含む）も ADR-0167 と ADR-0170 の後で古く、#0124 で再実行する。#0124 を先に閉じれば、
ここで走るのは OpenMC の段階と `compare-1B` だけになる。
同じく main の OpenMC のスクリプトの変更（openmc_core_2d.py、openmc_model.py など）で、`openmc-5A-2D-extract` の deps も古い。
抽出だけの段階（statepoint は `references/openmc/vera-5A-2D/`）なので、同じ機会に走らせ直し、5-2D の OpenMC との比較の値が
変わらないことを確かめる。

## 完了の条件

- `dvc repro examples/vera/dvc.yaml` で古い段階（`openmc-1B`、`openmc-1B-extract`、`compare-1B`、`openmc-5A-2D-extract`）を
  実行し、`dvc.lock` をコミットする（Slurm で排他、スキル openmc-reference）
- 5-2D の OpenMC との比較の値が変わらないことを確かめる
- examples/vera/verification.md の「参照解」の 1B の行と、「OpenMC との比較」の 1B の差を、新しい参照解で書き直す。
  2026-10-11 の参照解（examples/vera/verification.md の「参照解」の 1B、9c121906、1.181783 ± 0.000063）との差が統計誤差の中かを書く

## 経緯

- 2026-10-14: 高速化の後の再測定で `dvc status` が 2 段階を changed と出し、分けた
