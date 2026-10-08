# M3 炉心 HZP

## 目的

M2 までで、燃料と案内管だけの 17×17 集合体を格子計算で解ける。M3 では、零出力（HZP）の炉心を 3D NEM で解くのに
要る機能を、xs・lattice・core にわたって揃える。対象は VERA Problem 3〜5 の炉心（Pyrex、制御棒、支持格子、
ノズル、炉心板、バッフル、反射体）である。検証は機能ごとに付けるが、進め方の軸は機能の実装とする。

## 範囲

- 物質は VERA Problem 3〜5 に要るものに限る: Pyrex、AIC（Ag-In-Cd）、B4C、SS304、Inconel、支持格子、
  シンブルプラグ、計装用シンブル。可燃性毒物（Gd、IFBA、WABA）は燃焼が本題なので M4 で扱う
- 新しい核種は C、Si、P、Ti、Mn、Ni、Ag、Cd、In（同位体ごと）。形式を増やす作業は見込まない（[ADR-0021](../adr/0021-follow-njoy-processing-flow.md)）
- 炉心の外の構造: [ADR-0078](../adr/0078-put-baffle-in-reflector.md)
- 天然の C: [ADR-0079](../adr/0079-split-natural-carbon.md)
- 熱水力のフィードバックは扱わない（[ADR-0029](../adr/0029-defer-thermal-hydraulic-feedback.md)）。燃料・被覆管・減速材の温度は入力で与える

## 決めたこと

### 核データ

- 565 K の熱中性子散乱と核種ごとの温度: [ADR-0080](../adr/0080-interpolate-sab-in-temperature.md)（[ADR-0047](../adr/0047-choose-temperatures-per-nuclide.md) を置き換え）
- 群の間の移行と MF4 の角度分布: [ADR-0081](../adr/0081-use-mf4-angular-distributions.md)（[ADR-0041](../adr/0041-transfer-by-two-body-kinematics.md) を置き換え）
- MF6 の角度のデータ: [ADR-0082](../adr/0082-leave-out-mf6-angular-data.md)（[ADR-0056](../adr/0056-treat-mf6-law2-law3-as-two-body.md) を置き換え）

### 格子計算

- セルの形と長方形の格子、CMFD の粗メッシュ: [ADR-0083](../adr/0083-generalize-to-rectilinear-lattice.md)（[ADR-0059](../adr/0059-use-square-lattice-geometry.md)、[ADR-0061](../adr/0061-accelerate-moc-with-cmfd.md) を置き換え）
- 出力に、集合体の面ごとの不連続因子（ADF）、ピン出力の再構成用の形状関数、反射体の少数群定数を加える。
  反射体の定数の作り方: [ADR-0084](../adr/0084-make-reflector-constants-with-fuel.md)
- 状態の分岐とテーブル: [ADR-0085](../adr/0085-branch-on-core-states.md)
- 構造材の塊の自己遮蔽: [ADR-0086](../adr/0086-shield-structural-lumps.md)
- 影の計算の全断面積: [ADR-0087](../adr/0087-use-dilute-totals-in-shadows.md)
- 等価なセルの枠のリング: [ADR-0088](../adr/0088-add-frame-ring-to-equivalent-cell.md)
- 2項有理近似の重み: [ADR-0089](../adr/0089-weight-two-terms-with-potential-scattering.md)
- AIC と B4C は光学的に非常に厚い。そのような領域では k がトラック間隔に対して不規則に動く（[ADR-0077](../adr/0077-keep-vera-discretization.md)）。
  制御棒のセルを入れたところで影響を測り、必要ならトラックの置き方を見直す

### 炉心計算

- 集合体ごとに違う少数群定数（装荷図）と、軸方向の層（燃料、支持格子、プレナム、ノズル、炉心板、上下の水）を
  入力で与える。ノードの定数は状態（ホウ素濃度、温度、制御棒）でテーブルを内挿する
- NEM に ADF を入れ、CMFD で固有値の反復を加速する（[ADR-0027](../adr/0027-solve-core-with-nem.md)、[ADR-0028](../adr/0028-accelerate-core-with-nonlinear-cmfd.md)。粗メッシュのソルバは [ADR-0003](../adr/0003-extract-shared-computation-on-second-use.md)）
- 制御棒の先端がノードの途中にあるときの扱い（カスピング）と、ピン出力の再構成の方法は、実装の前に文献で
  比べて決める
- 全炉心を解く: [ADR-0090](../adr/0090-solve-full-core.md)

### 検証

- 新しい物質の格子計算: VERA 2E〜2K と 2Q（KENO の値）と、同じ核データの OpenMC（連続エネルギー）
- ADF、反射体の定数、NEM: 2D の 4-2D（カラーセット）、REF1-2D（反射体）、5-2D（1/4 炉心）。参照解は OpenMC（2D）と
  KENO の値
- 3D の炉心: Problem 3、4 は OpenMC（3D）と KENO の値。Problem 5（全炉心の零出力試験）は、仕様書の実測値（臨界
  ホウ素濃度、制御棒バンクの価値、等温温度係数）と KENO の値で比べる。3D の全炉心の OpenMC は計算時間が大きいので
  作らない

## 進め方

1. データファイルの型と MessagePack の読み書き、Gauss-Legendre と Legendre 多項式を、形式のクレートと共通の
   計算のクレートに移す（段階の間を渡る型が増える前に行う。[ADR-0002](../adr/0002-keep-stage-crates-independent.md)、[ADR-0003](../adr/0003-extract-shared-computation-on-second-use.md)）
2. 新しい核種の処理と、MF4・MF6 の角度分布の散乱行列への反映
3. 熱中性子の散乱行列の複数の温度と、lattice での温度の内挿
4. セルの形の追加（中実の棒、支持格子の枠、長方形の領域、場所で違う格子の幅）
5. 格子計算の出力（ADF、ピン出力の形状関数、反射体の定数）と状態の分岐のテーブル
6. core: 装荷図と軸方向の層、テーブルの内挿、ADF を入れた NEM、CMFD（粗メッシュのソルバを共通のクレートに切り出す）
7. core: 制御棒の部分挿入、ピン出力の再構成
8. Problem 3〜5 の入力と、炉心の結果のまとめ
