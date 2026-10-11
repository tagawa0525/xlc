# M0: 3段階の縦切り

ENDF/B-VIII.0 から、多群ライブラリ（`xs`）、ピンセルの2群定数（`lattice`）、炉心の k_eff（`core`）までを
順に計算する例。手法と近似は M0 で決めた次の ADR のとおり（後に置き換えたものは、各 ADR の `superseded-by` が
新しい決定を指す）。範囲と完了条件は [M0 の計画](../../docs/plans/m0-vertical-slice.md)。

- 進め方とクレート: [ADR-0030](../../docs/adr/0030-build-a-vertical-slice-first.md)、[ADR-0001](../../docs/adr/0001-split-into-three-stage-crates.md)（`xlc-core` の名前を含む）
- 核データと処理の流れ: [ADR-0020](../../docs/adr/0020-use-endf-b-viii-0.md)、[ADR-0021](../../docs/adr/0021-follow-njoy-processing-flow.md)（M0 は Reich-Moore だけ）、[ADR-0022](../../docs/adr/0022-use-wims-69-groups.md)、[ADR-0039](../../docs/adr/0039-weight-with-njoy-iwt4.md)
- 断面積と散乱: [ADR-0038](../../docs/adr/0038-use-zero-kelvin-cross-sections-in-m0.md)、[ADR-0040](../../docs/adr/0040-scatter-isotropically-in-m0.md)、[ADR-0041](../../docs/adr/0041-transfer-by-two-body-kinematics.md)、[ADR-0042](../../docs/adr/0042-collapse-chi-by-fission-source.md)
- 格子計算: [ADR-0024](../../docs/adr/0024-solve-lattice-with-2d-moc.md)、[ADR-0043](../../docs/adr/0043-track-pin-cells-cyclically.md)、[ADR-0044](../../docs/adr/0044-collapse-to-two-groups.md)
- 炉心計算: [ADR-0027](../../docs/adr/0027-solve-core-with-nem.md)、[ADR-0045](../../docs/adr/0045-expand-nem-to-fourth-order.md)
- 段階間のファイル: [ADR-0009](../../docs/adr/0009-write-text-intermediate-files.md)

## 実行

ENDF/B-VIII.0 の中性子サブライブラリを `references/endf/ENDF-B-VIII.0_neutrons/` に置いて（`xs.toml` と `dvc.yaml` が
この場所を指す）、DVC で実行する（ADR-0123）。入力かコードが変わった段階だけが走る。

```sh
dvc repro examples/m0/dvc.yaml
```

生成物（`wims69.msgpack`, `pincell-fg.msgpack`, `core.out`）はこのディレクトリに書き出す。git では追跡せず、DVC の
キャッシュに置く。各段階の入力・コード・出力の md5 は `dvc.lock` に書かれる。

| ファイル          | 内容                                                                                                  |
| ----------------- | ----------------------------------------------------------------------------------------------------- |
| `xs.toml`         | 多群ライブラリの入力（4 核種、293.6 K、無限希釈だけ）                                                 |
| `pincell.toml`    | ピンセルの入力（3.1 wt% UO2 と水、被覆管なし）                                                        |
| `core.toml`       | 炉心の計算のケース（ピンセルの2群定数による均質な裸の直方体 340 × 340 × 366 cm）。入力の形は ADR-0105 |
| `cycle.toml`      | サイクルの装荷（256 体の ID の並び）                                                                  |
| `plant.toml`      | 炉心の形（16 × 16 体の位置、ピッチ 21.25 cm、高さ 366 cm）                                            |
| `assemblies.toml` | 台帳（ピンセルの2群定数の集合体 256 体）                                                              |
| `dvc.yaml`        | 3段階（xs → lattice → core）の入力、コード、出力とコマンド                                            |

## 結果

| 段階      | 結果                                                      | 計算時間               |
| --------- | --------------------------------------------------------- | ---------------------- |
| `xs`      | H-1, O-16, U-235, U-238 の WIMS 69群ライブラリ（293.6 K） | 49 s（核種ごとに並列） |
| `lattice` | ピンセルの k∞ = 0.7292                                    | 0.6 s                  |
| `core`    | 炉心の k_eff = 0.7259（核分裂中性子あたりの漏れ 0.0043）  | 13 s                   |

各段階の検証は、次のとおり。

- `xs`: ポイントワイズの断面積は NJOY RECONR と許容誤差 0.1% 以内で一致する。群平均断面積は NJOY GROUPR と 0.11% 以内、
  ν̄ は 2e-5 以内、χ は 3e-4 以内で一致する（`crates/xs/examples/` の比較プログラム）
- `lattice`: 均質なセルの k が無限媒質の解析解と 1e-6 で一致し、1群に縮約すると νΣf/Σa がセルの k∞ に一致する
- `core`: 均質な裸の直方体の k が2群拡散の解析解に収束する（各方向8ノードで 1.1e-4）

## M0 の近似による限界

k の値は実機（k∞ 約 1.3）よりかなり低い。主な原因は次のとおりで、M1 以降で置き換える。

- 無限希釈の断面積を使っており、U-238 の共鳴の自己遮蔽がない（共鳴吸収を大きく過大評価する）
- 熱中性子散乱は自由ガス模型（S(α,β) は M1 で使えるようにしたが、この例では使わない）
- 散乱は等方（P0）で輸送補正がなく、拡散係数は D = 1/(3Σt)（高速群で小さめになる）
- 被覆管がない
- 炉心の計算に加速がない（CMFD は M3）。大型炉心では反復に時間がかかる

## 装荷図の入力（ADR-0105）に書き換えた炉心（d4355e5、2026-10-09）

炉心の入力を、ケース（`core.toml`）、サイクル（`cycle.toml`）、炉心の形（`plant.toml`）、台帳（`assemblies.toml`）の
4つに分けた。体系は前と同じ均質な裸の直方体で、径方向は集合体の位置ごとに1ノード（16 × 16、21.25 cm 角）、
軸方向は 18 等分。燃焼を扱わないので、集合体の ID は通し番号（P001〜P256）にした。

| 入力                                     | k_eff    | 核分裂中性子あたりの漏れ | 計算時間 |
| ---------------------------------------- | -------- | ------------------------ | -------- |
| 4つの入力（このコミット）                | 0.723150 | 0.0110                   | 9.3 s    |
| 前の1つの入力（main の 6d77883 の core） | 0.723150 | 0.0110                   | —        |

同じ少数群定数（このコミットの lattice で作った `pincell-fg.msgpack`）で、前の入力と k と漏れが一致する。上の「結果」の
表（k_eff = 0.7259、漏れ 0.0043）は M0 のときの lattice の定数によるもので、その後の lattice の変更（拡散係数の
定義など）で定数が変わった。差の内訳はこの節では調べていない。
