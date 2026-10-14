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

Problem 5（零出力試験）の計算の仕組み（名前つきの状態を並列に解くこと、比較のツール、臨界ホウ素濃度の探索）は
ADR-0146、ADR-0145 で入った。3D の全炉心の入力は骨組みだけで（examples/vera/p5/core-hzp.toml に、バンクを
挿入しない場合の状態だけ）、まだ解けない。残っているのは次。

- 炉心の形と装荷: 5-2D の装荷図（ブランチ feat/p5-2d の examples/vera/p5/core-plant.toml、未マージ）を、軸方向の
  反射体の層を持つ 3D の炉心の形にする。サイクル（core-hzp-cycle.toml）と台帳も要る
- 軸方向の層: Problem 3 の層（ブランチ feat/p3-xlc の examples/vera/p3/core-3B-assemblies.toml の形、未マージ）に
  合わせ、Pyrex の本数と濃縮度ごとの設計に広げる。軸方向のメッシュは層の境界に乗せる（ADR-0105）
- 定数: 濃縮度（2.11、2.619、3.10 %）と Pyrex の本数と計装管の組ごとの分岐のテーブル（ホウ素、減速材、燃料の軸。
  ホウ素の格子点は #0038。場合のホウ素濃度 1170〜1291 ppm と探索の範囲 1100〜1400 ppm を含むこと）。制御棒の軸（AIC と B4C）と、制御棒が入る軸方向の
  反射体の層の制御棒の軸（ADR-0127、#0030）。ライブラリの核種は #0073
- バンクの位置: 状態に書く形（ADR-0127、#0063）が入ったら、バンクを挿入する場合（1、3〜10、12〜19、22〜31）の
  状態を足し、compare_hzp.py の状態の検査（いまは挿入した場合をエラーにする）をバンクの位置を比べる形にする
- DVC のステージ: core（`core core-hzp.toml > core-hzp.out`）と比較（`python3 compare_hzp.py <spec> core-hzp.toml
  core-hzp.out > compare-hzp.out`）を examples/vera/p5/dvc.yaml に足す。dvc.yaml は feat/p5-2d が作るので、
  そのマージの後に足す（いま作ると、両方のブランチが同じファイルを足すうえ、入力のない段階が `dvc repro -P` を
  止める）

関係する ADR: ADR-0145、ADR-0146

## 完了の条件

core-hzp.toml の全 32 の場合（と、場合 1、2 の臨界ホウ素濃度の探索）を解き、compare_hzp.py の結果を
examples/vera/verification.md に書く。

## 経緯

- 2026-10-12: 計算の仕組みと比較のツールを入れたブランチ feat/p5-hzp-tests で起こした
- 2026-10-13: 最小の入力の試行（445da315、feat/p3-xlc の xs.toml に 560、570 K を足したライブラリ）で、3.10% の燃料の層の格子 `31-fuel` と `31-fuel-aic` が分岐の 8 状態を通った（k∞ は 1100 ppm、560 K で 1.206516 と 0.866747、1400 ppm、570 K で 1.171749 と 0.845511）。残りの格子と炉心は、格子の自己遮蔽の高速化（perf/0039-shielding-time）のあとに再実行する
- 2026-10-14: 3.10% の最小の入力の格子 19 段階（単一集合体 7、反射体 6、制御棒のテーブル 2、内挿と 2D の小さな炉心 4）を
  7abbc031（分岐点 704b59e5）の実行ファイルで解いた。main（a66351f2）へ rebase して実行ファイルが変わった（ADR-0157、
  ADR-0158。bin/lattice の md5 bd5af095）ので、8cbedf55 の trials/p5-hzp/same_state.py で、`31-fuel` と `reflector-I-31-fuel`
  の最初の状態（1100 ppm、560 K、560 K）を新しい実行ファイルで解き直し、前の出力の同じ状態とビット単位で同じことを確かめた
  （lattice/31-fuel.msgpack の md5 3a2eef2a、lattice/reflector-I-31-fuel.msgpack 016d2181、
  lattice/reflector-I-31-fuel-fuel.msgpack 9809e8dd）。19 段階は流し直さず、dvc.lock に記録する。平板（`axial-31`、
  `axial-31-rcca`）は ADR-0158 で解き方が変わるので、新しい実行ファイルで解き直した
