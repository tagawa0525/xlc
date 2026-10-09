# M2 集合体計算

## 目的

M1 までの格子計算は、無限希釈の断面積でピンセルを1つ解くだけである。M2 では、17×17 の PWR 燃料集合体を
2D MOC で解き、共鳴の自己遮蔽を入れる。集合体の均質化した2群定数と、ピンごとの出力分布を出す。

## 検証の対象

- C5G7（2D）: 7群の巨視的断面積が与えられたベンチマーク。自己遮蔽を含まないので、MOC の形状・トラッキング・
  真空境界・反復解法だけを検証できる
- VERA Problem 1B〜1D（ピンセル）と 2B〜2D（17×17 集合体、毒物なし）: 減速材と被覆管は 600 K、
  燃料は 600, 900, 1200 K。この温度はすべて M1 のライブラリにある。565 K の 1A, 2A は、
  H in H2O の S(α,β) に 565 K がないので対象にしない
- 参照解は、同じ核データから手元で作った OpenMC で求める（コード同士の検証。[ADR-0031](../adr/0031-verify-against-openmc-with-same-data.md)）。VERA の公開値（KENO-VI、
  ENDF/B-VII.0）とも比べるが、核データの差を含むので参考とする

## 決めたこと

### 材料と核データ

- 材料の組成: [ADR-0054](../adr/0054-model-vera-materials.md)
- R-matrix limited 形式（Fe-54, 57）: [ADR-0055](../adr/0055-reconstruct-r-matrix-limited-subset.md)
- MF6 の LAW=2, 3: [ADR-0056](../adr/0056-treat-mf6-law2-law3-as-two-body.md)（M3 で [ADR-0082](../adr/0082-leave-out-mf6-angular-data.md) に置き換え）
- 出射分布の補間則: [ADR-0057](../adr/0057-interpolate-outgoing-distributions-as-endf.md)
- 非分離共鳴の LRF=1 と Fe-58 の l=3: [ADR-0058](../adr/0058-read-lrf1-as-lrf2-with-l3.md)

### 形状とトラッキング

- 形状: [ADR-0059](../adr/0059-use-square-lattice-geometry.md)（M3 で [ADR-0083](../adr/0083-generalize-to-rectilinear-lattice.md) に置き換え）
- トラッキング: [ADR-0060](../adr/0060-track-whole-assembly-cyclically.md)、境界条件: [ADR-0093](../adr/0093-choose-lattice-boundary-conditions.md)

### 反復の加速と並列化

- CMFD による加速: [ADR-0061](../adr/0061-accelerate-moc-with-cmfd.md)（M3 で [ADR-0083](../adr/0083-generalize-to-rectilinear-lattice.md) に置き換え）
- optimally diffusive CMFD: [ADR-0062](../adr/0062-use-optimally-diffusive-cmfd.md)
- CMFD の後の更新と境界の初期値: [ADR-0063](../adr/0063-update-angular-flux-after-cmfd.md)
- スイープの並列化と再現性: [ADR-0064](../adr/0064-sum-sweeps-per-azimuth.md)
- 減衰の計算: [ADR-0066](../adr/0066-compute-attenuation-per-track.md)（[ADR-0065](../adr/0065-tabulate-attenuation.md) を置き換え）
- スレッド数の既定: [ADR-0067](../adr/0067-default-threads-to-physical-cores.md)

### 自己遮蔽

- ピンごとの減速計算（採らなかった方式を含む）: [ADR-0068](../adr/0068-solve-pinwise-slowing-down.md)（[ADR-0023](../adr/0023-shield-by-bondarenko-and-equivalence.md) を置き換え）
- 減速計算の範囲: [ADR-0069](../adr/0069-cap-cell-slowing-down-at-20-kev.md)
- 散乱核の標的の扱い: [ADR-0070](../adr/0070-treat-targets-at-rest.md)
- SPH 因子: [ADR-0071](../adr/0071-correct-with-sph-factors.md)
- ライブラリの σ0 ごとの重み: [ADR-0072](../adr/0072-weight-shielded-data-by-flux-calculator.md)（[ADR-0051](../adr/0051-weight-shielded-data-by-nr.md) を置き換え）
- ライブラリのポイントワイズの断面積、分離共鳴の上端、ポテンシャル散乱: [ADR-0073](../adr/0073-store-pointwise-data-in-library.md)
- 範囲の外の遮蔽: [ADR-0074](../adr/0074-shield-outside-range-by-equivalence.md)

### 散乱の非等方性

- MOC の非等方散乱: [ADR-0075](../adr/0075-solve-moc-with-anisotropic-scattering.md)（[ADR-0040](../adr/0040-scatter-isotropically-in-m0.md) を置き換え）
- ライブラリの散乱行列のモーメント: [ADR-0076](../adr/0076-store-legendre-scattering-moments.md)
- MF6 の角度のデータ: M2 の計画の扱いは実装の前に変えた（[ADR-0082](../adr/0082-leave-out-mf6-angular-data.md)）

### 離散化の設定

- OpenMC と比べる集合体は方位角 30・間隔 0.03 cm、ほかは VERA の入力の設定: [ADR-0111](../adr/0111-use-30-azimuths-for-vera-assemblies.md)（[ADR-0077](../adr/0077-keep-vera-discretization.md) を置き換え）

### 出力

- 集合体を均質化した2群定数（ピンセルと同じ形式）と、ピンごとの核分裂率の分布を出す。
  ADF とピン出力の再構成用のデータは、炉心側で使う M3 で加える

## 進め方

6 の集合体の UFG の参照の道具（#0019）を除いて済（最後のマージは fd9a7a7）。

1. 一般の正方格子の形状、集合体全体のトラッキング、真空境界と、巨視的断面積の入力。C5G7 で検証する
2. ファイルの約束事への移行（入力を TOML、データファイルを MessagePack、エラーを thiserror と anyhow に。
   [ADR-0010](../adr/0010-write-inputs-in-toml.md)、[ADR-0011](../adr/0011-write-data-files-in-messagepack.md)、[ADR-0016](../adr/0016-type-errors-with-thiserror-and-anyhow.md)）。コードが増える前に済ませる
3. CMFD による加速とスレッドによる並列化
4. Zircaloy-4 の核種（Sn, Fe, Cr）の核データ処理（R-matrix limited 形式、MF6 の LAW=2, 3）
5. OpenMC の準備と、VERA Problem 1B〜1D, 2B〜2D の参照解
6. 自己遮蔽（ピンごとの UFG の減速計算、SPH 因子、集合体の UFG の参照解）。VERA Problem 1B〜1D で検証する。
   等価なセルの誤差は、集合体全体を UFG の MOC で直接解く参照の道具で、案内管や集合体の端のそばのピンで測る。
   OpenMC とは、中性子束の重みの群定数と直接比べる
7. 17×17 集合体の入力、ピンごとの出力分布、集合体の2群定数。VERA Problem 2B〜2D で検証する
8. 散乱の非等方性（ライブラリの散乱行列の P_L のモーメントと、MOC の非等方の散乱源）。VERA Problem 2B〜2D の
   水の穴のまわりのピンの核分裂率で検証する
