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
- 運動学のしきい値より下から始まる反応: [ADR-0100](../adr/0100-move-tables-to-kinematic-thresholds.md)

### 格子計算

- セルの形と長方形の格子、CMFD の粗メッシュ: [ADR-0083](../adr/0083-generalize-to-rectilinear-lattice.md)（[ADR-0059](../adr/0059-use-square-lattice-geometry.md)、[ADR-0061](../adr/0061-accelerate-moc-with-cmfd.md) を置き換え）
- 出力に、集合体の面ごとの不連続因子（ADF）、ピン出力の再構成用の形状関数、反射体の少数群定数を加える。
  反射体の定数の作り方: [ADR-0084](../adr/0084-make-reflector-constants-with-fuel.md)
- 状態の分岐とテーブル: [ADR-0085](../adr/0085-branch-on-core-states.md)
- 構造材の塊の自己遮蔽: [ADR-0086](../adr/0086-shield-structural-lumps.md)
- 影の計算の全断面積: [ADR-0087](../adr/0087-use-dilute-totals-in-shadows.md)
- 等価なセルの支持格子のリング: [ADR-0088](../adr/0088-add-frame-ring-to-equivalent-cell.md)
- 2項有理近似の重み: [ADR-0089](../adr/0089-weight-two-terms-with-potential-scattering.md)
- AIC と B4C は光学的に非常に厚い。そのような領域では k がトラック間隔に対して不規則に動く（[ADR-0111](../adr/0111-use-30-azimuths-for-vera-assemblies.md)）。
  制御棒のセルを入れたところで影響を測り、必要ならトラックの置き方を見直す

### 炉心計算

- 集合体ごとに違う少数群定数（装荷図）と、軸方向の層（燃料、支持格子、プレナム、ノズル、炉心板、上下の水）を
  入力で与える。ノードの定数は状態（ホウ素濃度、温度、制御棒）でテーブルを内挿する
- NEM に ADF を入れ、CMFD で固有値の反復を加速する（[ADR-0027](../adr/0027-solve-core-with-nem.md)、[ADR-0028](../adr/0028-accelerate-core-with-nonlinear-cmfd.md)。粗メッシュのソルバは [ADR-0003](../adr/0003-extract-shared-computation-on-second-use.md)）。面の条件の形と外周の真空の境界は [ADR-0101](../adr/0101-impose-transport-current-ratio-at-vacuum-faces.md)（[ADR-0099](../adr/0099-couple-core-nodes-by-get-with-vacuum-boundary.md) を置き換え）
- 制御棒の先端がノードの途中にあるときの扱い（カスピング）: [ADR-0127](../adr/0127-weight-partially-rodded-nodes-by-axial-fine-mesh-flux.md)（proposed）。ピン出力の再構成: [ADR-0096](../adr/0096-reconstruct-pin-power-analytically.md)
- 全炉心を解く: [ADR-0090](../adr/0090-solve-full-core.md)
- Problem 5 の計算: 臨界ホウ素濃度は入力の範囲を囲い込んだまま Illinois 法で探す（[ADR-0145](../adr/0145-search-critical-boron-by-illinois-in-input-bracket.md)）。
  1 つのケースに名前つきの状態を複数書いて並列に解き、バンクの価値、DBW、ITC は比較のツールが状態の k の差から求める（[ADR-0146](../adr/0146-solve-named-states-of-a-case-in-parallel.md)）

### 検証

- 新しい物質の格子計算: VERA 2E〜2K と 2Q（KENO の値）と、同じ核データの OpenMC（連続エネルギー）
- ADF、反射体の定数、NEM: 2D の 4-2D（カラーセット）、REF1-2D（反射体。形状を表すのに #0067 が要る）、5-2D（1/4 炉心）。参照解は OpenMC（2D）と
  KENO の値
- 3D の炉心: Problem 3、4 は OpenMC（3D）と KENO の値。Problem 5（全炉心の零出力試験）は、仕様書の実測値（臨界
  ホウ素濃度、制御棒バンクの価値、等温温度係数）と KENO の値で比べる。3D の全炉心の OpenMC は計算時間が大きいので
  作らない

## 進め方

項目ごとに状態（済・作業中・未着手）と、済んだ部分のマージのコミットを書く。細かい作業は milestone: M3 の issue で追う。

1. データファイルの型と MessagePack の読み書き、Gauss-Legendre と Legendre 多項式を、形式のクレートと共通の
   計算のクレートに移す（段階の間を渡る型が増える前に行う。[ADR-0002](../adr/0002-keep-stage-crates-independent.md)、[ADR-0003](../adr/0003-extract-shared-computation-on-second-use.md)）
   — 済（2891d70）
2. 新しい核種の処理と、MF4 の角度分布の散乱行列への反映（MF6 の角度のデータは入れない。[ADR-0082](../adr/0082-leave-out-mf6-angular-data.md)）
   — 済（MF4 は 97eea1e、MF6 を入れない決定は af65bdc、分離共鳴のない核種のライブラリは 44070e2。新しい核種（C、Si、P、
   Ti、Mn、Ni、Ag、Cd、In の 29 核種）のライブラリと NJOY との比較は cce0b6e と 007c227（examples/m3-nuclides/verification.md）。
   減速計算の範囲より上の散乱源を NJOY にそろえたのは 770d09d、運動学のしきい値から反応を始めるのは e694dca（[ADR-0100](../adr/0100-move-tables-to-kinematic-thresholds.md)）。
   VERA の k はこの2つで 6 桁まで変わらない（63fcd3d）。しきい値の値を 0 にしたのは 92ef88e（#0024）。NJOY との差のうち、拡げた格子が小さい反応で粗いこと（#0022）は
   未解決。新しい物質の格子計算（2E〜2K）はまだ）
3. 熱中性子の散乱行列の複数の温度と、lattice での温度の内挿 — 済（92b1e89）
4. セルの形の追加（中実の棒、支持格子の枠、長方形の領域、場所で違う格子の幅） — 作業中（枠は 07eb7fb、長方形の格子は
   42c0e67、構造材の自己遮蔽は d91afe5 と 2f2a71b。中実の棒は同心円のセルで表せるが、シンブルプラグの物質（SS304）は
   新しい核種が要り、まだ確かめていない）
5. 格子計算の出力（ADF、ピン出力の形状関数、反射体の定数）と状態の分岐のテーブル — 作業中（反射境界の ADF は c7a1ee2、
   アルベド境界は 21e52b5、真空の辺のある格子の ADF（1ノードの NEM。ADR-0094）は 7bce5d3、真空の辺の比 J/φ
   （ADR-0101）は 8cd9422。格子のノードごとの均質化と I 型の反射体の定数（ADR-0095。f*と、反射体の Koebke の等価な D（ADR-0104）、
   core の f* の面）は f6839b5 で、core が MOC の k を +0.0006 pcm で再現する（examples/reflector/verification.md、2fc6a8fe）。燃料が真空に
   面する 17×17 の集合体の ADF（#0027）、L 型（#0029。いまの格子の形では作れない。#0068）の反射体、非反射の格子の CDF（#0020）は未。軸方向の反射体の定数（ADR-0128。集合体の高さ全体の平板、燃料の媒質、`mix`、少数群定数の z の辺）は済み（VERA 3A の平板で core が MOC の k を再現する。examples/reflector/verification.md）。上の層の制御棒の軸と Problem 3 の比較は #0030。
   角型の反射体の格子（燃料に角の点だけで接するノードにも反射体の印を付ける）は済み（#0029。examples/reflector/verification.md）。鏡面反射の格子の形状関数と
   CDF（ADR-0096）は 0998ec8。分岐のテーブル（ADR-0097、ADR-0106）は ecec076 の VERA 2B での確かめまで済み（ホウ素の
   格子点は #0038））
6. core: 装荷図と軸方向の層、テーブルの内挿、ADF を入れた NEM、CMFD（粗メッシュのソルバを共通のクレートに切り出す）
   — 作業中（ADF は f62bbed、真空の境界は bf9ed9d と 8cc6a90、格子計算の比 J/φ を課す真空の面は 2a46331
   （ADR-0101）、両側が反射の1ノードの方向の発散の修正は b5f574b（#0025）、核分裂のないノードの ANM
   （ADR-0102）は 3700729、装荷図と軸方向の層（ADR-0105。炉心の形、台帳、サイクル、ケースの4つの入力）は
   7bda1ce。反射体の定数を装荷図の回転で東西南北に置くことの確認（4 つの向きで k とピン出力が
   同じ）は済み（#0029）。集合体を分割するノード（ADR-0108。#0033）は鏡面反射の燃料集合体まで済み、反射体の定数の分割は未。炉心の外形に沿った境界と、反射体どうしの面に各側の ADF を使う扱い（ADR-0132）は
   済み（#0034。5-2D での確認は #0029、反射体の列の上下の層は #0072）。CMFD は要否を Problem 3〜5 で決めるまで止めている（#0076）。
   テーブルの内挿（ケースの状態。ADR-0097、ADR-0146）は済）
7. core: 制御棒の部分挿入、ピン出力の再構成 — 作業中（ピン出力の再構成は ADR-0096 で、検証は #0043。制御棒の部分挿入は
   ADR-0127 で方式を決め、実装は #0063）
8. Problem 3〜5 の入力と、炉心の結果のまとめ — 作業中（Problem 5 の計算の仕組み: 臨界ホウ素濃度の探索（ADR-0145）、名前つきの状態を並列に解くケースと、実測値と KENO-VI の値と比べるツール examples/vera/p5/compare_hzp.py
   （ADR-0146）は済。3D の全炉心の入力（炉心の形、軸方向の層、定数、バンクの位置）は #0092。5-2D の全炉心の入力と段階、
   KENO-VI との比較は 5da6c0fb（examples/vera/verification.md の「Problem 5-2D の炉心の計算」。k は +235〜+294 pcm、集合体出力の
   RMS 4〜5% で中心が低く外周が高い。原因の切り分けは #0086）。OpenMC の参照解との比較は #0073 の本番のあと。L 型の位置の
   近似の大きさは #0029）
