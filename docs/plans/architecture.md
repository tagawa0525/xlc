# PWR 炉心解析コード 全体計画

## 目的

ENDF 形式の評価済み核データから出発し、PWR 炉心全体の解析までを自前で行うコード体系を作る。

## 全体構成

計算を3つのステージ（`xs`: 核データ処理、`lattice`: 格子計算、`core`: 炉心計算）に分け、ステージごとのクレートとして
実装する。ステージ間はファイルでデータを受け渡す。ステージが共有するデータファイルの形式と計算は、共通のクレート
（`xlc-format`、`xlc-numeric`、`xlc-depletion`）に置く。

- ステージのクレートとファイルでの受け渡し: [ADR-0001](../adr/0001-split-into-three-stage-crates.md)
- ステージのクレートの間の依存: [ADR-0166](../adr/0166-build-chain-in-its-own-crate.md)（[ADR-0162](../adr/0162-add-depletion-crate.md) を置き換え）
- 共有する計算の切り出し: [ADR-0003](../adr/0003-extract-shared-computation-on-second-use.md)
- ステージごとのソルバの設定: [ADR-0004](../adr/0004-keep-solver-settings-per-stage.md)
- 依存クレート: [ADR-0006](../adr/0006-use-de-facto-crates.md)（[ADR-0005](../adr/0005-prefer-std-and-own-numerics.md) を置き換え）
- 線形代数: [ADR-0007](../adr/0007-use-faer-for-linear-algebra.md)、小さな密行列: [ADR-0008](../adr/0008-keep-own-small-dense-solvers.md)

## ファイルの約束事

- 入力ファイル: [ADR-0010](../adr/0010-write-inputs-in-toml.md)
- データファイル: [ADR-0011](../adr/0011-write-data-files-in-messagepack.md)（[ADR-0009](../adr/0009-write-text-intermediate-files.md) を置き換え）
- 後方互換と版: [ADR-0012](../adr/0012-drop-backward-compatibility.md)
- 読み込みの厳密さ: [ADR-0013](../adr/0013-read-files-strictly.md)
- 単位と群の並び: [ADR-0014](../adr/0014-fix-units-and-group-order.md)
- 核種の名前: [ADR-0015](../adr/0015-name-nuclides-like-u-235.md)
- エラー: [ADR-0016](../adr/0016-type-errors-with-thiserror-and-anyhow.md)、panic の禁止: [ADR-0017](../adr/0017-forbid-panics-in-production.md)
- 出力とログ: [ADR-0018](../adr/0018-log-with-tracing-to-stderr.md)
- Python と開発環境: [ADR-0019](../adr/0019-use-python-only-for-verification.md)

## 手法の選択

### xs: 核データ処理

- 対象データ: [ADR-0020](../adr/0020-use-endf-b-viii-0.md)
- 処理の流れと共鳴の形式: [ADR-0021](../adr/0021-follow-njoy-processing-flow.md)
- 群構造: [ADR-0022](../adr/0022-use-wims-69-groups.md)
- 自己遮蔽: [ADR-0068](../adr/0068-solve-pinwise-slowing-down.md)（分離共鳴の範囲）、[ADR-0074](../adr/0074-shield-outside-range-by-equivalence.md)（範囲の外）。当初の方式は [ADR-0023](../adr/0023-shield-by-bondarenko-and-equivalence.md)

### lattice: 格子計算

- 輸送計算: [ADR-0024](../adr/0024-solve-lattice-with-2d-moc.md)
- 燃焼計算: [ADR-0025](../adr/0025-deplete-with-cram.md)
- 出力: [ADR-0026](../adr/0026-output-two-group-tables-from-lattice.md)

### core: 炉心計算

- 手法: [ADR-0027](../adr/0027-solve-core-with-nem.md)
- 固有値の加速: [ADR-0028](../adr/0028-accelerate-core-with-nonlinear-cmfd.md)
- 燃焼度と核種の追跡: [ADR-0035](../adr/0035-track-history-by-micro-depletion.md)
- 熱水力のフィードバック: [ADR-0029](../adr/0029-defer-thermal-hydraulic-feedback.md)

## 進め方

最初に、3つのステージを最小構成で端から端までつなぐ（縦切り。[ADR-0030](../adr/0030-build-a-vertical-slice-first.md)）。その後で各ステージを作り込む。

| マイルストーン  | 内容                                                                                                 | 検証                                        | 状態                                                          |
| --------------- | ---------------------------------------------------------------------------------------------------- | ------------------------------------------- | ------------------------------------------------------------- |
| M0 縦切り       | 少数核種（H-1, O-16, U-235, U-238）の無限希釈断面積 → ピンセル MOC → 均質な炉心の NEM                | 解析解、モンテカルロ計算                    | 済（85b5e79）                                                 |
| M1 核データ処理 | 共鳴の再構成・ドップラー拡がり・熱中性子散乱・f-table を備え、PWR に必要な全核種を処理する           | NJOY の出力との比較                         | 済（M1 の対象核種。4f04365。M3 の新しい核種は M3 で処理する） |
| M2 集合体       | 17×17 集合体の MOC と自己遮蔽                                                                        | C5G7, VERA Problem 1–2                      | 集合体の UFG の参照のツール（#0019）を除いて済（d17a530）     |
| M3 炉心 HZP     | 3D NEM、CMFD、ADF、ピン出力の再構成                                                                  | VERA Problem 3–5                            | 作業中（docs/plans/m3-core-hzp.md）                           |
| M4 燃焼         | 格子計算の燃焼と分岐計算（可燃性毒物を含む）、炉心での燃焼度追跡とミクロ燃焼（燃焼履歴の補正を含む） | VERA Problem 6 以降、格子燃焼のベンチマーク | 計画の下書き（docs/plans/m4-depletion.md）                    |

「燃焼まで」（M4）を最初の到達目標とする。熱水力のフィードバックを含む BEAVRS との比較は、その後に検討する。

格子計算の精度の目標は、同じ核データで作ったモンテカルロ計算の参照解に対して、燃焼を通じて k の差が ±100 pcm
（許容の上限は ±150 pcm）、ピン出力の差が rms 0.2% 以内とする（[ADR-0032](../adr/0032-target-lattice-accuracy.md)）。達成の状況は検証の例
（examples/vera/verification.md など）に記録する。

参照解は、同じ核データから手元で作った OpenMC で求める（[ADR-0124](../adr/0124-pin-reference-assets-with-nix.md)）。

## 今後の拡張

M4 の後に、動特性と、取り出した燃料の冷却計算（ORIGEN と同じく、放射能と崩壊熱を求める）を加える。

### M4 の分岐テーブルに入れるデータ

- 動特性と冷却計算が使うデータ: [ADR-0033](../adr/0033-store-kinetics-and-cooling-data-in-m4-tables.md)

### M4 で作るもの

- 燃焼チェーン: [ADR-0034](../adr/0034-build-full-depletion-chain.md)
- 燃焼履歴: [ADR-0035](../adr/0035-track-history-by-micro-depletion.md)
- core のノードごとの状態: [ADR-0036](../adr/0036-keep-core-state-per-assembly.md)
- 可燃性毒物（Gd 入りの燃料、IFBA、Pyrex、WABA）を検証の対象に含める。実際の PWR の集合体のほとんどに入っており、
  Gd は燃料の外側から順に燃えるので、燃料のリングの分割と燃焼のステップの設計に効くため
- 分岐計算は燃焼点の状態から始まるので、独立したケースの集まりではない。lattice の1回の実行の中で燃焼と分岐を
  行うか（CASMO と同じ）、燃焼点の状態をファイルで渡すかを M4 の計画で決める（燃焼のソルバの置き場所は [ADR-0003](../adr/0003-extract-shared-computation-on-second-use.md)）

### ソルバの形

- core の演算子と固有値の反復: [ADR-0037](../adr/0037-separate-operators-from-iterations.md)
- 燃焼のソルバ: [ADR-0003](../adr/0003-extract-shared-computation-on-second-use.md)

### 必要になってから加えるもの

いずれも上の設計の上に後から載せられる。

- 臨界ホウ素濃度の探索、制御棒の価値、反応度係数、停止余裕
- 随伴中性子束: 実効遅発中性子割合と即発中性子寿命、核データの共分散と組み合わせた感度・不確かさの解析に使う
- 検出器の応答: 実機のデータは炉内の検出器の反応率で与えられるので、比べるのに要る
- ガンマ線の輸送と発熱: 熱水力のフィードバックを入れる段階でピン出力に効く
