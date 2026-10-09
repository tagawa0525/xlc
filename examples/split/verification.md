# 集合体を分けたノード（VERA 2B）

集合体を n × n の部分ノードに分ける（[ADR-0108](../../docs/adr/0108-split-assemblies-with-shared-constants.md)）と、
定数を作った体系（鏡面反射の集合体）では分けないときと同じ解になることを確かめる例。範囲と作業の順序は
[M3 の計画](../../docs/plans/m3-core-hzp.md)。

| ファイル                            | 内容                                                                                                        |
| ----------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| `assembly-1.toml`                   | examples/vera/assembly-2B.toml の集合体の、形状関数を付けた少数群定数（`fg-1.msgpack`）と格子計算のピン出力 |
| `assembly-2.toml`                   | 同じ集合体を 2 × 2 の部分ノードに分けた少数群定数（`fg-2.msgpack`）                                         |
| `plant.toml`                        | 集合体を2体並べた炉心の形（外周は鏡面反射）                                                                 |
| `assemblies-*.toml`、`cycle-*.toml` | 台帳とサイクル。右の集合体を 90° 回して置く                                                                 |
| `case-1.toml`、`case-2.toml`        | ノードを分けない（1 × 1）ケースと 2 × 2 のケース。再構成したピン出力を書く                                  |

ライブラリは examples/vera の `vera.msgpack`。生成物（少数群定数、ピン出力、ログ）は git で追跡しない。

```sh
cd examples/split
../../target/release/lattice assembly-1.toml > lattice-1.log 2>&1   # 約 4 分（32 スレッド）
../../target/release/lattice assembly-2.toml > lattice-2.log 2>&1
../../target/release/core case-1.toml > core-1.log 2>&1
../../target/release/core case-2.toml > core-2.log 2>&1
```

## 分けないときとの一致（ee9e1b8、2026-10-09）

| ケース | k-eff    |
| ------ | -------- |
| 1 × 1  | 1.185801 |
| 2 × 2  | 1.185801 |

炉心全体のピン出力（17 行 × 34 列。2 × 2 では中央の行と列の半分のセルを1つに戻したもの）の差は、相対で最大 2.3e-8
（NEM の収束の程度）。回した集合体と回していない集合体が同じ列のノードに並ぶ（切ったセルの幅の丸めが違う）場合も
扱えた。

1 × 1 の core の左の集合体のピン出力と、格子計算の `pin-powers` の差は最大 0.17%。格子計算のピン出力がトラックから
推定した体積で平均していることによると見込む（#0041）。格子計算は2回で 508 秒かかった（#0039）。
