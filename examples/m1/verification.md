# M1: 温度と自己遮蔽に依存する核データ

ENDF/B-VIII.0 から、M1 の対象核種（U-234, U-235, U-236, U-238, O-16, O-17, H-1, B-10, B-11, Zr）の
温度と背景断面積 σ0 に依存する多群ライブラリ（`xs`）を作り、ピンセルの2群定数（`lattice`）を計算する例。
手法と近似は M1 で決めた次の ADR のとおり（後に置き換えたものは、各 ADR の `superseded-by` が
新しい決定を指す）。範囲と作業の順序は [M1 の計画](../../docs/plans/m1-nuclear-data.md)。

- 共鳴の形式（MLBW の追加）: [ADR-0021](../../docs/adr/0021-follow-njoy-processing-flow.md)
- ドップラー拡がりと温度: [ADR-0046](../../docs/adr/0046-broaden-with-sigma1.md)、[ADR-0047](../../docs/adr/0047-choose-temperatures-per-nuclide.md)
- 非分離共鳴: [ADR-0048](../../docs/adr/0048-shield-unresolved-like-unresr.md)、[ADR-0049](../../docs/adr/0049-integrate-width-fluctuations-adaptively.md)
- σ0 ごとの群定数と重み: [ADR-0050](../../docs/adr/0050-tabulate-sigma0-per-nuclide-and-temperature.md)、[ADR-0051](../../docs/adr/0051-weight-shielded-data-by-nr.md)
- 熱中性子散乱: [ADR-0052](../../docs/adr/0052-scatter-h-in-h2o-by-sab.md)
- `xs` の入力: [ADR-0053](../../docs/adr/0053-read-xs-input-from-file.md)

## 実行

ENDF/B-VIII.0 の中性子と熱中性子散乱のサブライブラリを `references/endf/ENDF-B-VIII.0_neutrons/` と
`references/endf/ENDF-B-VIII.0_thermal_scatt/` に置いて（別の場所なら `xs.toml` の `endf-directory` と
`thermal-directory` を書き換えて）実行する。

```sh
examples/m1/run.sh
```

生成物（`m1.msgpack`, `pincell-fg.msgpack`）はこのディレクトリに書き出す。git では追跡しない。

| ファイル       | 内容                                                                               |
| -------------- | ---------------------------------------------------------------------------------- |
| `xs.toml`      | 多群ライブラリの入力（14 核種、核種ごとの温度、σ0 9 点、H in H2O の S(α,β)）       |
| `pincell.toml` | ピンセルの入力（M0 と同じ 3.1 wt% UO2 と水、被覆管なし。燃料 900 K、減速材 600 K） |
| `run.sh`       | 2段階を順に実行するスクリプト                                                      |

## 結果

| 段階      | 結果                                                                | 計算時間                   |
| --------- | ------------------------------------------------------------------- | -------------------------- |
| `xs`      | 14 核種の WIMS 69群ライブラリ（燃料・被覆管 4 温度、冷却材 2 温度） | 10.4 min（核種ごとに並列） |
| `lattice` | ピンセルの k∞ = 0.7257（無限希釈）                                  | 0.8 s                      |

ピンセルの k∞ の温度による変化（同じライブラリ、無限希釈）:

| 燃料    | 減速材  | k∞     |
| ------- | ------- | ------ |
| 293.6 K | 293.6 K | 0.7284 |
| 900 K   | 293.6 K | 0.7278 |
| 293.6 K | 600 K   | 0.7262 |
| 900 K   | 600 K   | 0.7257 |

上の 2 つの表は、[ADR-0133](../../docs/adr/0133-scatter-h-in-h2o-up-to-the-evaluated-sab-limit.md) の前（H in H2O の S(α,β) の上限 4 eV、表の外は 0、その上は静止した
標的）のライブラリの値で、測り直していない（#0070）。

同じ 293.6 K の M0（0.7292）とは、水素の熱中性子散乱（M0 は自由ガス、M1 は S(α,β)）だけが違う。
格子計算は M2 まで無限希釈の値を使うので、燃料温度の効果（ドップラー効果）は自己遮蔽のある場合よりずっと小さい。
水の密度は入力のまま（0.74 g/cm3）で、減速材の温度は主に熱中性子散乱に効く。

## 検証

各段階は NJOY の同じ条件の出力と比べた（`crates/xs/examples/` の比較プログラム）。

| 段階                             | NJOY との一致                                                   |
| -------------------------------- | --------------------------------------------------------------- |
| ドップラー拡がり（BROADR）       | 2e-6                                                            |
| 非分離共鳴の平均断面積（RECONR） | 非分離共鳴の範囲で 2e-4                                         |
| 非分離共鳴の自己遮蔽（UNRESR）   | 4e-4                                                            |
| σ0 ごとの群平均断面積（GROUPR）  | 1e-4（U-238 の群 13 の捕獲だけ 6e-4）                           |
| S(α,β) の熱中性子散乱（THERMR）  | 断面積 2.1e-3、出射分布 1.1%                                    |
| 解析式の出射分布（GROUPR）       | O-17 の (n,2n) などの散乱行列は断面積の差（1.4e-3）の範囲で一致 |

ずれの原因はそれぞれ調べてあり、NJOY 側の近似（BROADR の区間探索、GROUPR の台形公式、UNRESR の求積表、
THERMR の S(α,β) の補間）による。U-236 の核分裂スペクトル χ は、NJOY が 1e-5 eV での θ のスペクトルを
15 MeV までの全入射エネルギーに使うため、高エネルギーの群で十数 % ずれる（こちらは ENDF の θ(E) に従う）。

### S(α,β) の熱中性子散乱（THERMR）

上の表の S(α,β) の行は、H in H2O の 293.6 K を、S(α,β) の上限を評価済みデータの上限 10.00008 eV、表の外を THERMR と同じ
短時間衝突近似にして（[ADR-0133](../../docs/adr/0133-scatter-h-in-h2o-up-to-the-evaluated-sab-limit.md)）、THERMR（許容誤差 0.001、上限 10.00008 eV）と GROUPR の
MF6 MT222 と比べたもの（`compare_njoy_thermr`、references/njoy-runs/thermr/H1-10eV-293.6。c9246058、2026-10-11）。上限より
下の群で、行の和の差は最大 2.1e-3（3.3〜4 eV）、行の形の差（差の絶対値の和 / 行の和）は最大 1.1e-2（0.996〜1.02 eV）。

```sh
cargo run --release -p xlc-xs --example compare_njoy_thermr -- references/njoy-runs/thermr/H1-10eV-600/tape20 \
  references/njoy-runs/thermr/H1-10eV-600/tape30 references/njoy-runs/thermr/H1-10eV-600/tape26 125 600 10.00008
```

600 K（references/njoy-runs/thermr/H1-10eV-600）では、行の和の差は 4〜9.877 eV で −2.9e-3、3.3〜4 eV で −1.8e-3（xlc の方が
小さい）、0.97〜1.12 eV で +1.1〜+1.5e-3。4〜10 eV の差の大半は補間の違いによる: THERMR は ln S を 3 点の 2 次式で補間し
（thermr.f90 の `terpq`）、xlc は ENDF の補間則（ln S を α、β について線形）に従う（ADR-0133）。同じ表と同じ短時間衝突近似で
両方の補間の散乱断面積を求めると、THERMR の補間の方が 3.6 eV で 0.17%、4.5 eV で 0.21%、7 eV で 0.25%、9.8 eV で 0.21% 大きく、
1.1 eV で 0.04%、0.05 eV で −0.15% 違う（`sab_interpolation.py`、14c70e60、2026-10-11）。THERMR の許容誤差を 1e-4 にしても、
THERMR の散乱断面積は 0.05〜9.8 eV で 0.02% 以内しか変わらない（H1-10eV-600-tight の PENDF。GROUPR は止まった）ので、
許容誤差による差ではない。1 eV より下の差（±1.5e-3）は補間の違いだけでは説明できず、確かめていない。

```sh
python3 examples/m1/sab_interpolation.py references/endf/ENDF-B-VIII.0_thermal_scatt/tsl-HinH2O.endf 600 3.6 4.5 7 9.8
```

## 計算時間について

`xs` の計算時間の大半は U-236 の非分離共鳴の自己遮蔽（1 温度あたり約 3 分）である。50 keV 以上では
競合反応の幅もゆらぐので、幅のゆらぎの期待値が3重の適応積分になり、1 系列あたりの被積分関数の評価が
約 1000 万回になる（競合幅のない系列では約 15 万回）。
