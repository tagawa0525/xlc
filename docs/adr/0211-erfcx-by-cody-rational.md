---
status: proposed
date: 2026-10-28
requires: [ADR-0203, ADR-0209]
supersedes: []
superseded-by:
issues: ["#0181"]
---

# ADR-0211: erfcx の x ≥ 2.5 を Cody の有理近似で求める

## 背景

ADR-0209 の既定では、自己遮蔽の時間の大半が自由ガスの核の D の散乱源で、1C の 0 MWd/kg のピンセルの 1 スレッドのプロファイルの 39% が
erfcx（x ≥ 2.5 の 40 段の連分数）だった（#0181 の「問題」）。2A の燃焼のジョブは ADR-0209 の前の 18 分 58 秒から 1 時間 47 分 59 秒になった。

erfcx は xs の `special::erfcx` と、その複製の lattice の `ufg::free_gas::erfcx` にある（ステージのクレートは互いに依存しない。ADR-0203）。
xs では `resonance/fluctuation.rs`（非分離共鳴の幅の揺らぎの積分）、`multigroup/free_gas.rs`（自由ガスの核）が erfcx を、
`multigroup/outgoing.rs`（Maxwell 分布の累積 `gamma_p_3_2`）が erf を通して x ≥ 2.5 の erfcx を使う。

判断の軸:

- 正確さ: 倍精度の多倍長の参照との相対誤差（ulp）で、今の連分数と同じかそれより小さい
- 速さ: 1 スレッドの `free_gas` の時間
- 保守: 依存は標準ライブラリだけ。係数を手で書くなら、出典と、係数を確かめた方法を残す

## 検討した案

精度は、x の 2.5〜4、4〜6 の一様な乱数 4 万点ずつ、6〜50 の 2 万点、50〜1e8 の対数で一様な 5,000 点と区間の端（計 10 万 5 千点）で、
mpmath 1.4.1 の 60 桁の erfc(x)·exp(x²) との差を、参照値の倍精度の ulp で測った（最大 / 二乗平均）。速さは、2.5 ≤ x < 50 の 10 万点を
200 周して 1 回あたりの時間（Zen 5、`target-cpu=znver5`、1 コア）。どちらも手元の使い捨ての計測（2026-10-28）。

| 案                                    | 2.5〜4    | 4〜6      | 6〜50     | 50〜1e8   | 1 回    |
| ------------------------------------- | --------- | --------- | --------- | --------- | ------- |
| 今の連分数（40 段）                   | 2.53/0.71 | 2.88/0.85 | 2.97/0.78 | 2.59/0.80 | 68.5 ns |
| Cody の有理近似（CALERF のまま）      | 5.48/1.27 | 1.48/0.52 | 1.43/0.49 | 1.40/0.49 | 3.2 ns  |
| Cody、2.5〜4 を補償 Horner 法（採用） | 1.83/0.53 | 1.48/0.52 | 1.43/0.49 | 1.40/0.49 | 6.9 ns  |
| Faddeeva package の erfcx             | 2.98/0.81 | 3.98/1.05 | 2.36/0.59 | 2.68/0.69 | 3.0 ns  |

- 案 A: 連分数の段の数を x で減らす。今の連分数と 200 段の連分数がビットで同じになる段は、x = 2.5 で 38、4 で 21、6 で 14。連分数を呼ぶ x の
  78% が 2.5 ≤ x < 6（#0181）なので、減らせる段は少ない。採らない
- 案 B: Cody の有理 Chebyshev 近似（W. J. Cody, "Rational Chebyshev approximations for the error function", Math. Comp. 23 (1969)
  631–637）。Netlib の specfun の CALERF（1990-03-19 版、パブリックドメイン）が倍精度の係数を持つ。erfcx を 0.46875 ≤ x ≤ 4 で 8 次 / 8 次の
  有理式、x > 4 で (1/√π − w R(w))/x（w = 1/x²、R は 4 次 / 4 次）で求める。有理式そのものの誤差は、mpmath で係数の 10 進の値のまま
  評価して、2.5 ≤ x ≤ 4 で相対 8.1e−19（0.007 ulp）、4 ≤ x ≤ 4e4 で 5.9e−19。CALERF の Horner 法のままでは、2.5〜4 で丸めが積もって
  最大 5.5 ulp になり、連分数より粗い（FMA の Horner 法でも 3.9 ulp）
- 案 C: 案 B の 2.5〜4 の分子と分母を、補償 Horner 法（Graillat, Langlois, Louvet, "Compensated Horner scheme", 2005）で求める。
  各段の積と和の丸めを誤差なしの変換（積は融合積和、和は Knuth の TwoSum）で取り出して足し戻すので、係数と x が正の多項式は丸め 1 回
  ほどの誤差になる。x > 4 は R の項が 1/√π に比べて小さく、丸めが薄まるので、FMA の Horner 法のままで連分数より細かい
- 案 D: Faddeeva package（S. G. Johnson、MIT ライセンス）の erfcx。x < 50 は y = 4/(4 + x) の 100 区間の 6 次の Chebyshev 多項式
  （係数 700 個）、x ≥ 50 は連分数の 5 段。速いが、4〜6 で最大 4.0 ulp と連分数より粗く、表が大きく、ライセンスの表示が要る。採らない

## 決定と理由

案 C にする。2.5 ≤ x ≤ 4 は Cody の 8 次 / 8 次の有理式を補償 Horner 法で、x > 4 は Cody の 1/x² の有理式を FMA の Horner 法で求め、
x < 2.5 の級数は変えない。最大の誤差がすべての範囲で連分数より小さく（1.83 ulp 以下。連分数は 2.97 ulp）、1 回の時間は 1/10。

- 積和は `f64::mul_add` で書く。融合積和は 1 回の丸めなので、FMA の命令の無い CPU（ライブラリの fma）でも同じビットになる
  （`target-cpu` を付けずにコンパイルした計測でも、テストの 34 点のビットが同じ）
- 係数は CALERF の倍精度の DATA 文の 10 進の値を丸めずにコピーし、上の mpmath の評価で確かめた。コードのコメントに出典と確かめ方を書く
- xs と lattice に同じコードを置く（ADR-0203 の複製の理由のまま）。xs の値も丸めの桁で変わるので、ライブラリを作り直す。作り直すのは
  1C と 2A の燃焼の再測定に要る examples/depletion（65 核種）と examples/depletion/p2（90 核種）だけで、ほかの例のライブラリは `bin/xs` が
  変わって changed になるが、結果が要るときに作り直す（ADR-0130）

## 確認

- `erfcx_is_within_two_ulp_of_multiprecision_values`（xs の `special::tests`、lattice の `ufg::free_gas::tests`）: 区間の端と、上の
  計測で誤差の大きかった x を含む 30 点で、mpmath の 60 桁の値との差が 2 ulp 以下（今の連分数は x = 4.463 で 2.88 ulp なので落ちる）
- `erfcx_gives_the_same_bits_as_its_copy`（同じ 2 つのテスト）: 2 つのクレートに同じ 34 点のビットの表を置き、xs と lattice の erfcx が
  同じ値を返す
- 1 スレッドの `free_gas`（1C の 0 MWd/kg のピンセル `issue-0133/pin-1C-0.toml` に ADR-0209 の既定を明示したもの、手元、2026-10-28）:
  9.69 → 5.93 秒。k の標準出力は 6 桁で同じ（1.172387）

## 帰結

- ライブラリの値と lattice の自由ガスの核が丸めの桁で変わる。examples/depletion と examples/depletion/p2 のライブラリを作り直し、1C と 2A の
  燃焼を再測定する（examples/depletion/verification.md）
- 自由ガスの散乱源の残りの手間は #0181 に書く
