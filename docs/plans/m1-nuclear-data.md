# M1 核データ処理の作り込み

## 目的

M0 の `xs` は、0 K・無限希釈・自由ガスの断面積しか作れない。M1 では、PWR の燃料・被覆管・冷却材の核種について、
温度と自己遮蔽に依存する多群ライブラリを作れるようにする。格子計算側で自己遮蔽を使うのは M2 とし、
M1 は核データ処理だけを対象にする。

## 対象核種

U-234, U-235, U-236, U-238, O-16, O-17, H-1（H in H2O）, B-10, B-11, Zr-90, Zr-91, Zr-92, Zr-94, Zr-96

被覆管はジルカロイを Zr だけで近似する。Sn, Fe, Cr は必要になった時点で追加する（M2 で追加した。[ADR-0054](../adr/0054-model-vera-materials.md)）。

## 決めたこと

- 共鳴の形式（MLBW の追加）: [ADR-0021](../adr/0021-follow-njoy-processing-flow.md)
- ドップラー拡がり: [ADR-0046](../adr/0046-broaden-with-sigma1.md)
- 核種ごとの温度: [ADR-0047](../adr/0047-choose-temperatures-per-nuclide.md)
- 非分離共鳴の自己遮蔽: [ADR-0048](../adr/0048-shield-unresolved-like-unresr.md)
- 幅のゆらぎの期待値: [ADR-0049](../adr/0049-integrate-width-fluctuations-adaptively.md)
- σ0 ごとの群定数とライブラリの中身: [ADR-0050](../adr/0050-tabulate-sigma0-per-nuclide-and-temperature.md)
- σ0 ごとの群定数の重み: [ADR-0051](../adr/0051-weight-shielded-data-by-nr.md)（M2 で [ADR-0072](../adr/0072-weight-shielded-data-by-flux-calculator.md) に置き換え）
- 熱中性子散乱: [ADR-0052](../adr/0052-scatter-h-in-h2o-by-sab.md)
- `xs` の入力: [ADR-0053](../adr/0053-read-xs-input-from-file.md)

## 検証

各段階を NJOY（BROADR, UNRESR, THERMR, GROUPR）の同じ条件の出力と比べる。M0 と同じく、
`references/` の NJOY で参照解を作り、`crates/xs/examples/` の比較プログラムで比べる。

## 進め方

1. MLBW 形式
2. ドップラー拡がり
3. 狭い共鳴の近似による σ0 依存の群定数
4. 非分離共鳴の自己遮蔽
5. S(α,β) による熱中性子散乱
6. ライブラリの形式の更新と、M1 の対象核種によるライブラリの作成
