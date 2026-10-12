---
status: open
labels: [core, vera]
opened: 2026-10-12
closed:
milestone: M3
blocked-by: ["#0063"]
resolution:
---

# Problem 5 の 3D の全炉心の入力（定数、層、バンクの位置）

## 問題

Problem 5（零出力試験）の計算の仕組み（名前つきの状態を並列に解くこと、比較の道具）は ADR-0146 で入った。
臨界ホウ素濃度の探索は別のブランチ（feat/p5-boron-search）で入れる。3D の全炉心の入力は骨組みだけで（examples/vera/p5/core-hzp.toml に、バンクを
挿入しない場合の状態だけ）、まだ解けない。残っているのは次。

- 炉心の形と装荷: 5-2D の装荷図（ブランチ feat/p5-2d の examples/vera/p5/core-plant.toml、未マージ）を、軸方向の
  反射体の層を持つ 3D の炉心の形にする。サイクル（core-hzp-cycle.toml）と台帳も要る
- 軸方向の層: Problem 3 の層（ブランチ feat/p3-xlc の examples/vera/p3/core-3B-assemblies.toml の形、未マージ）に
  合わせ、Pyrex の本数と濃縮度ごとの設計に広げる。軸方向のメッシュは層の境界に乗せる（ADR-0105）
- 定数: 濃縮度（2.11、2.619、3.10 %）と Pyrex の本数と計装管の組ごとの分岐のテーブル（ホウ素、減速材、燃料の軸。
  ホウ素の格子点は #0038。場合のホウ素濃度 1170〜1291 ppm と臨界ホウ素濃度の探索の範囲を含むこと）。制御棒の軸（AIC と B4C）と、制御棒が入る軸方向の
  反射体の層の制御棒の軸（ADR-0127、#0030）。ライブラリの核種は #0073
- バンクの位置: 状態に書く形（ADR-0127、#0063）が入ったら、バンクを挿入する場合（1、3〜10、12〜19、22〜31）の
  状態を足し、compare_hzp.py の状態の検査（いまは挿入した場合をエラーにする）をバンクの位置を比べる形にする
- DVC の段階: core（`core core-hzp.toml > core-hzp.out`）と比較（`python3 compare_hzp.py <spec> core-hzp.toml
  core-hzp.out > compare-hzp.out`）を examples/vera/p5/dvc.yaml に足す。dvc.yaml は feat/p5-2d が作るので、
  そのマージの後に足す（いま作ると、両方のブランチが同じファイルを足すうえ、入力のない段階が `dvc repro -P` を
  止める）

関係する ADR: ADR-0146

## 完了の条件

core-hzp.toml の全 32 の場合（と、場合 1、2 の臨界ホウ素濃度の探索）を解き、compare_hzp.py の結果を
examples/vera/verification.md に書く。

## 経緯

- 2026-10-12: 計算の仕組みと比較の道具を入れたブランチ feat/p5-hzp-tests で起こした
