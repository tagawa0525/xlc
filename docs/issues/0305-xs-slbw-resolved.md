---
type: Issue
status: open
labels: [xs, depletion]
opened: 2026-11-20
closed:
milestone:
blocked-by: []
resolution:
---

# xs が読めない MF2 の形式（ENDF/B-VII.0 の NIS > 1、VIII.1 の KBK = 1）が、燃焼チェーンの中性子のファイルに残る

## 問題

- 中性子のファイルの全部を 1 つずつ `xs range` に通すと（1 核種だけの入力。2026-11-20、xlc a15e0a0f）、止まるファイルは次の 3 つで、
  分離共鳴の SLBW（LRF=1）で止まるファイルは無い（ADR-0265）。ファイルの数は VII.0 が 393、VII.1 が 423、VIII.0 が 557、VIII.1 が 558
  - VII.0 の天然の V と Zn: 複数の同位体（NIS = 2、5）の MF2 で、「NIS=2 isotopes are not supported」
  - VIII.1 の Sr-88: 分離共鳴のバックグラウンドの R 行列（KBK = 1）で、「background R-matrix (KBK=1) ... are not supported」
  - 「no group bounds」で止まるファイル（VII.0 の 5、VII.1 の 8、VIII.0 の 9、VIII.1 の 9）は、1 核種だけの入力で分離共鳴の上端が熱の
    打ち切りより低いためで、読めないのではない
  - 再現: 1 核種だけを書いた入力（`[nuclide.X-1]` に `endf` と `temperatures`）で `xs range` を実行する
- `xs range` は MF2 の分離共鳴の上端までしか読まない。`xs nuclide` が読むもの（MF1〜MF8）の読み込みと形式の検査（`xs check`、ADR-0266）に
  全部のファイルを通しても、上の 3 つは同じ文で止まる。ほかの形式で止まるファイルと形式ごとの issue は、examples/depletion/verification.md の
  「中性子のファイルの読み込みの走査」。上の 3 つは #0292 の 9 核種に当たらない
- 3 つのファイルは、燃焼のライブラリの 65 核種に入らない

関係する ADR: ADR-0264、ADR-0265

## 完了の条件

- [x] VII.0 と VII.1 の版ごとのライブラリ（examples/depletion/endfb70、endfb71）の 65 核種が、分離共鳴の SLBW で止まらない（ADR-0265。
  ほかの形式で止まる核種は #0307）
- [ ] 燃焼チェーンの全核種で、上の形式のどれをどう扱うかを決める（#0303 の完了の条件の 2 つ目と合わせる）

## 経緯

- 2026-11-20: 版ごとのピンセルの燃焼の準備（#0303）で、VII.0 と VII.1 のライブラリの `xs range` が SLBW で止まり、起こした
- 2026-11-20: SLBW を ENDF-102 の式で読むようにし（ADR-0265）、中性子のファイルの全部を通し直して、残る形式を今の姿にした
