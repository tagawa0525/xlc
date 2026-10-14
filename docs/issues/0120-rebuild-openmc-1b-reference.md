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

main の ADR-0167（MLBW のエネルギーに依存する散乱半径）で `bin/xs` が変わり、`xs` の段階も古い（再測定のジョブは ADR-0167 の
マージの前の分岐点 4b644566 で走った。examples/core、m0、reflector の `xs` とその下流も同じく古い）。完了の条件のコマンドは
`xs`（約 26 分）と `lattice-1B` も実行する。ADR-0167 は NRO=1 の MLBW の評価だけに効き、VERA の核種はそれまでの xs で処理
できていたので、vera.msgpack と 1B の xlc の k は変わらない見込み（確かめていない）。変わったら、参照解の変化と分けて書く。
同じく main の OpenMC のスクリプトの変更（openmc_core_2d.py、openmc_model.py など）で、`openmc-5A-2D-extract` の deps も古い。
抽出だけの段階（statepoint は `references/openmc/vera-5A-2D/`）なので、同じ機会に走らせ直し、5-2D の OpenMC との比較の値が
変わらないことを確かめる。

## 完了の条件

- `dvc repro examples/vera/dvc.yaml` で古い段階（`openmc-1B`、`openmc-1B-extract`、`compare-1B`、`openmc-5A-2D-extract`、
  古ければ `xs` と `lattice-1B`）を実行し、`dvc.lock` をコミットする（Slurm で排他、スキル openmc-reference）
- 5-2D の OpenMC との比較の値が変わらないことを確かめる
- examples/core、m0、reflector の古い `xs` とその下流も `dvc repro` で揃え、出力が変わるかを記録に書く
- examples/vera/verification.md の「参照解」の 1B の行と、「OpenMC との比較」の 1B の差を、新しい参照解で書き直す。
  2026-10-11 の参照解（1.181783 ± 0.000063）との差が統計誤差の中かを書く

## 経緯

- 2026-10-14: 高速化の後の再測定で `dvc status` が 2 段階を changed と出し、分けた
