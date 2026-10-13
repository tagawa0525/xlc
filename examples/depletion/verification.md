# 燃焼のソルバ（CRAM）と燃焼チェーン

`xlc_depletion::deplete`（[ADR-0161](../../docs/adr/0161-solve-depletion-with-ipf-cram.md)）を、解析解のある小さなチェーンと、OpenMC の
`openmc.deplete.cram` と比べる（[M4 の計画](../../docs/plans/m4-depletion.md) の進め方の 4）。全核種の燃焼チェーン
（[ADR-0164](../../docs/adr/0164-build-full-chain-from-endf.md)）を、OpenMC の `Chain.from_endf` で同じ ENDF から作ったチェーンと比べる（進め方の 2。「燃焼チェーン」）。

| ファイル                                           | 内容                                                                                                   |
| -------------------------------------------------- | ------------------------------------------------------------------------------------------------------ |
| `cases.toml`                                       | OpenMC と比べる燃焼行列（例示の PWR の燃料のチェーン、複素の固有値の環、硬いチェーン、300 核種の乱数） |
| `compare_openmc.py`                                | 同じ行列を OpenMC の `CRAM16`、`CRAM48` と xlc に与えて差を出す（references の devShell）              |
| `../../crates/depletion/examples/deplete_text.rs`  | 行列を文字で受け取って `deplete` を呼ぶ（`compare_openmc.py` が呼ぶ）                                  |
| `chain.toml`、`dvc.yaml`                           | 全核種のチェーン `chain.msgpack`（`bin/chain`）と OpenMC のチェーン `openmc-chain.xml` を作る段階      |
| `compare_chain.py`                                 | 2 つのチェーンの核種、半減期、崩壊の分岐、反応、収率を比べる                                           |
| `chain-cases.toml`、`compare_chain_depletion.py`   | 2 つのチェーンに同じ 1 群の反応率を与え、1 ステップを CRAM48 で解いて比べる                            |
| `../../crates/depletion/examples/deplete_chain.rs` | チェーンと反応率から燃焼行列を組み立てて `deplete` を呼ぶ（`compare_chain_depletion.py` が呼ぶ）       |

誤差と差は、絶対値の最大を初期の数密度の最大で割ったもの（ADR-0161 の「確認」）。

## 解析解との比較

`cargo test -p xlc-depletion`（crates/depletion/src/cram/tests.rs）。閾値は ADR-0161 の「確認」。

### 計算

bb018540（分岐点は main の 321efe1c）、2026-10-14。`cargo test -p xlc-depletion --release -- --nocapture` の `ERROR` の行の、
テストと次数ごとの最大。

| チェーン                                          | 次数 16 | 次数 48 |
| ------------------------------------------------- | ------- | ------- |
| 1 核種の崩壊（x ∈ [1e-6, 1e5]）                   | 1.1e-15 | 2.2e-15 |
| 直列の 3 核種（‖A t‖ = 1e9 の硬いチェーンを含む） | 1.9e-15 | 1.9e-15 |
| 分岐                                              | 2.2e-16 | 7.2e-16 |
| 複素の固有値の環（kt = 0.1〜100）                 | 6.1e-12 | 1.0e-15 |
| 親の 1e-9 倍の娘（相対誤差）                      | 1.7e-15 | 4.1e-16 |

次数 16 は、固有値が負の実軸から離れた環でだけ誤差が開く。

## OpenMC との比較

```sh
cargo build --release -p xlc-depletion --example deplete_text
nix develop .#references -c python3 examples/depletion/compare_openmc.py \
    examples/depletion/cases.toml target/release/examples/deplete_text
```

### 計算

1ceb8b8f（分岐点は main の 321efe1c）、2026-10-14。Slurm のジョブ 314（2 CPU、排他でない）。

| 行列                                                          | 次数 16 の差 | 次数 48 の差 | OpenMC の 16 と 48 の差 |
| ------------------------------------------------------------- | ------------ | ------------ | ----------------------- |
| pwr-fuel（14 核種、Pu-239 の α 崩壊の環、30 日）              | 1.8e-18      | 3.2e-16      | 3.2e-16                 |
| cycle（固有値が複素の環、kt = 10）                            | 5.6e-16      | 2.8e-16      | 6.1e-12                 |
| stiff（‖A t‖ = 1e9 のチェーン）                               | 4.4e-16      | 3.9e-16      | 1.0e-15                 |
| random-300（300 核種、環を含む乱数の行列、5 列に 100 個の娘） | 1.6e-16      | 4.1e-16      | 5.1e-16                 |

すべて閾値 1e-13 の内にある。値が初期の最大の 1e-12 倍より大きい核種の相対差の最大は、random-300 で 1.8e-9（16）と
2.3e-10（48）、ほかは 6.8e-13 以下。random-300 の相対差は、初期の最大の 1e-12 倍に近い小さな核種で、絶対の差（丸めの水準）を
小さな値で割ったもの。OpenMC の 16 と 48 の差は環の行列でだけ 6.1e-12 に開き、解析解の環の次数 16 の誤差（6.1e-12）と合う。

## 燃焼チェーン

```sh
dvc repro examples/depletion/dvc.yaml          # chain.msgpack と openmc-chain.xml
cargo build --release -p xlc-depletion --example deplete_chain
cd examples/depletion
nix develop ../..#references -c python3 compare_chain.py chain.msgpack openmc-chain.xml
nix develop ../..#references -c python3 compare_chain_depletion.py chain-cases.toml chain.msgpack openmc-chain.xml \
    ../../target/release/examples/deplete_chain
```

### 計算

72a9a554（分岐点は main の 1f50ae7a）、2026-10-14。チェーンの段階は同じ `bin/chain`（md5 db8179c6）と OpenMC 0.15.2 の
`Chain.from_endf`（ENDF/B-VIII.0 の decay 3,821、nfy 31、neutrons 557）。`bin/chain` は 1 スレッドで約 1.2 秒、`Chain.from_endf` は
約 10 秒（手元。排他のジョブではない）。

### チェーンの比較

| 項目                       | 結果                                                                                                      |
| -------------------------- | --------------------------------------------------------------------------------------------------------- |
| 核種                       | どちらも 3,820（同じ集合）。不安定 3,561、反応を持つ 555、核分裂する 88                                   |
| 半減期                     | 不安定の 3,561 核種で一致（相対差 0）。半減期が 0 の不安定な 14 核種（Xe-136 など）はどちらも安定とみなす |
| 崩壊の分岐（型と娘が同じ） | 5,056。分岐比の差の最大 2.1e-7（Bi-196m2）、1e-12 を超えるもの 377                                        |
| 崩壊の分岐（娘が違う）     | SF 127、娘が decay にない 8                                                                               |
| 反応（生成物が同じ）       | 2,430                                                                                                     |
| 反応（生成物が違う）       | xlc が MF8 の状態に分けるもの 104、生成物が decay にない 1（Pu-245 (n,p)）                                |
| 核分裂の収率（自分の収率） | 31 核種の 52 のエネルギーの点で、0 でない生成物の数と値が一致（差 0）                                     |
| 核分裂の収率（借りる）     | 57 核種。借りる先が同じなのは 4（Am-242 の Am-242m など）                                                 |

違いはすべて ADR-0164 の規則による:

- 分岐比の差: 和の 1 からのずれ（最大 9e-7）を、xlc は比で割り直し、OpenMC は最大の分岐に足す
- SF: xlc は親を減らし生成物を持たない。OpenMC は娘を親そのものにするので、燃焼行列で SF の分岐が親を減らさない
- 娘が decay にない 8 件（Cf-239、Cf-256、Es-240、Es-243、Es-258、Rf-253、Ds-279m、Ni-48）: xlc は娘を生まず、OpenMC は
  `replace_missing` で β の向きを推して別の核種にする。He-5 の α の娘（中性子）は、どちらも追わない
- MF8 の状態: xlc は MF8 が LMF = 9、10 で挙げる状態ごとのチャネルに分け（Am-241 (n,γ) の Am-242m、Np-237 (n,2n) の Np-236m、
  Am-243 (n,γ) の Am-244m など）、OpenMC の `from_endf` は基底状態だけにする。xlc は 12 件で、MF8 の励起エネルギーと decay の
  準安定状態の励起エネルギーが 1e-3 を超えて違う（Am-244m が 72.9 keV と 86.1 keV、Np-236m が 185 keV と 56 keV など）ので、
  近い順の対応で決め、2 件（Cr-51 (n,α) の Ti-48 の 3.63 MeV、Ir-191 (n,2n) の Ir-190 の 36.15 keV）は decay に対応する
  準安定状態がなく基底状態に落とす（`chain` の警告に全件が出る）
- Pu-245 (n,p) の Np-245: xlc は残りの核を生まず（H-1 だけ）、OpenMC は β− の向きに進めて Pu-245 にする
- 借りる収率: xlc は同じ原子番号の質量数の最も近い核種（U-239 は U-238、Np-239 は Np-238、Pu-243 は Pu-242）、OpenMC は
  isotone か U-235（U-239 は Pu-241、Np-239 は Pu-240、Pu-243 は Cm-245）
- decay にない収率の生成物（U-235 で 112 の状態）は、どちらも準安定状態を基底状態へ、基底状態を β− の向きに移し、値が一致する

### 同じ反応率での 1 ステップ

反応率と収率のエネルギー（0.0253 eV）は chain-cases.toml。全核種のチェーン（3,820 核種）の燃焼行列を、xlc は `Chain::matrix` と
`deplete`（次数 48）で、OpenMC は `Chain.form_matrix` と `CRAM48` で解く。差は、絶対値の最大を初期の数密度の最大で割ったものと、
初期の最大の 1e-12 倍より大きな核種の相対差。

| 場合                    | 差 / 初期の最大 | 経路の核種の相対差の最大                                         | 相対差の最大（全核種） |
| ----------------------- | --------------- | ---------------------------------------------------------------- | ---------------------- |
| U-238 の捕獲、30 日     | 2.8e-15         | 2.8e-15（U-239、Np-239、Pu-239）                                 | 1.4e-11（U-235m）      |
| U-235 の核分裂、10 時間 | 6.5e-15         | 9.6e-15（I-135、Xe-135、Pm-149、Sm-149、Cs-135、Xe-136、Sm-150） | 1.0e-7（Ge-81m）       |
| U-235 の核分裂、30 日   | 1.4e-13         | 9.3e-15（同上）                                                  | 1.0e-7（Ge-81m）       |
| Pu-239 の核分裂、30 日  | 4.5e-13         | 1.1e-15（Pu-239 と同上）                                         | 1.0e-7（Ge-81m）       |

経路の核種は丸めの水準で一致し、核種の並び、崩壊と反応と収率の向きが OpenMC と同じ行列になる。Ge-81m の 1.0e-7 は、先行核
Ga-81 の分岐比の和が 1 − 1e-7 で、xlc は比で割り直し（Ge-81m への分岐が 1e-7 増える）、OpenMC は差を最大の分岐（Ge-81 への
分岐）に足すため。U-235m の 1.4e-11（値は初期の最大の 1.3e-12 倍）の原因は確かめていない。この 4 つの場合では、SF と MF8 の
状態の違いは効かない（U-238 の SF は崩壊の 5.5e-7 で、崩壊そのものが 30 日で 1e-11 の水準）。
