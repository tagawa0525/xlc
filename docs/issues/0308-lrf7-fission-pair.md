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

# R 行列の制限形（LRF=7）の核分裂の粒子対（MT=18）を読めない

## 問題

- ENDF/B-VIII.1 の Pu-239（n-094_Pu_239.endf）の分離共鳴は LRF=7 で、捕獲（MT=102）、弾性散乱（MT=2）に加えて核分裂（MT=18、
  PNT=-1）の粒子対を持つ。xs（crates/xs/src/resonance/r_matrix_limited.rs）は中性子と γ の粒子対だけを扱い（ADR-0055）、
  `xs nuclide` が「particle pair 3: MA=1 MT=18 is not supported」で止まる（Slurm 989 のログ
  /home/tagawa/github/slurm-logs/depletion-endfb81-989.out、main の f6360d711）。このため VIII.1 の燃焼のライブラリが作れず、
  版ごとのピンセルの燃焼の比較（#0303）が VIII.1 で止まっている
- 4 つの版（VII.0、VII.1、VIII.0、VIII.1）の中性子の全評価のうち、LRF=7 の評価は 28（VII.0 は 0）。その粒子対と制御の値を
  ファイルから数えると（2026-11-20）、核分裂の粒子対を持つのは VIII.1 の Pu-239 だけ。ほかに xs が LRF=7 で止まる評価は次の
  とおりで、この issue では直さない
  - Cl-35（VII.1、VIII.0、VIII.1）と Ca-40（VIII.0、VIII.1）: 荷電粒子の粒子対（MT=600、800。#0264 に挙げてある）。VII.1 の Cl-35 は
    シフト因子（SHF=-1）も持つ
  - VIII.1 の Cu-63: 弾性散乱の粒子対がシフト因子（SHF=1）を持つ
  - VIII.1 の Sr-88: spin group がバックグラウンドの R 行列（KBK=1、2）と位相シフトの指定（KPS=2）を持つ（粒子対ではなく spin group の形式）
- 版ごとのピンセルの燃焼の VIII.1 の 65 核種（ブランチ feat/0303-pin-depletion-by-endf-version の examples/depletion/endfb81/xs.toml）
  には、上の 4 核種は入っていない

関係する ADR: なし

## 完了の条件

- [ ] LRF=7 の核分裂の粒子対を ENDF-102 の定義で読み、NJOY RECONR の点ごとの断面積と比べる
- [ ] VIII.1 の Pu-239 の `xs nuclide` が終了コード 0 で終わる

## 経緯

- 2026-11-20: VIII.1 の燃焼のライブラリの作成（Slurm 989）が Pu-239 で止まり、起こした
