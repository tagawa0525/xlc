# M0 縦切り計画

## 目的

`xs` → `lattice` → `core` の3段階を最小構成で端から端までつなぎ、段階間のデータ形式と
各段階の計算の骨格を確定する（[ADR-0030](../adr/0030-build-a-vertical-slice-first.md)）。物理の精度は M1 以降で上げる。ここで決めた簡略化は、
後のマイルストーンで置き換える前提のものである。

## 対象

- 核種は H-1, O-16, U-235, U-238 とする。被覆管（Zr）は含めず、ピンセルは燃料と水の2領域とする
- ENDF/B-VIII.0 を使う（[ADR-0020](../adr/0020-use-endf-b-viii-0.md)）

## xs: 核データ処理

- 共鳴と温度: [ADR-0038](../adr/0038-use-zero-kelvin-cross-sections-in-m0.md)（M1 で置き換え）。分離共鳴の形式: [ADR-0021](../adr/0021-follow-njoy-processing-flow.md)
- 群構造と重み関数: [ADR-0022](../adr/0022-use-wims-69-groups.md)、[ADR-0039](../adr/0039-weight-with-njoy-iwt4.md)
- 散乱: [ADR-0040](../adr/0040-scatter-isotropically-in-m0.md)（M2 で [ADR-0075](../adr/0075-solve-moc-with-anisotropic-scattering.md) に置き換え）、[ADR-0041](../adr/0041-transfer-by-two-body-kinematics.md)（M3 で [ADR-0081](../adr/0081-use-mf4-angular-distributions.md) に置き換え）
- 核分裂: [ADR-0042](../adr/0042-collapse-chi-by-fission-source.md)

## lattice: 格子計算

- ピンセルの MOC と cyclic tracking: [ADR-0043](../adr/0043-track-pin-cells-cyclically.md)
- 散乱は等方とする（[ADR-0040](../adr/0040-scatter-isotropically-in-m0.md)）。2群への縮約と拡散係数: [ADR-0044](../adr/0044-collapse-to-two-groups.md)
- 検証: 均質なセルでの解析解（無限媒質の多群固有値）との一致と、トラック間隔・角度の収束性で行う

## core: 炉心計算

- 3D の2群 NEM: [ADR-0045](../adr/0045-expand-nem-to-fourth-order.md)。CMFD 加速は M3 で入れる（[ADR-0028](../adr/0028-accelerate-core-with-nonlinear-cmfd.md)）
- 境界条件は、反射とゼロ中性子束（解析解と比べるため）とする。真空境界（入射部分中性子流ゼロ）は
  M3 で入れる
- パッケージ名 `xlc-core`: [ADR-0001](../adr/0001-split-into-three-stage-crates.md)
- 検証: 均質な直方体炉心の解析解（2群拡散方程式、B² = Σ(π/L)²）との一致で行う

## 段階間のデータ形式

自前のテキスト形式とする（[ADR-0009](../adr/0009-write-text-intermediate-files.md)。M2 でデータファイルは [ADR-0011](../adr/0011-write-data-files-in-messagepack.md)、入力ファイルは [ADR-0010](../adr/0010-write-inputs-in-toml.md) に改めた）。

- 多群ライブラリ（`xs` → `lattice`）: 群境界、核種ごとの σt, νσf, σf, χ, 散乱行列
- 少数群定数（`lattice` → `core`）: 2群の Σt, Σa, νΣf, Σf, χ, 散乱行列, D

## 完了条件

ENDF から作った多群ライブラリで、ピンセルの k∞ と均質な炉心の k_eff を、
3段階のコマンドを順に実行して求められること。各段階の検証が上記のとおり通ること。

状態: 済（最後のマージは 85b5e79）。
