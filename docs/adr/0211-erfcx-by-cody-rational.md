---
status: proposed
date: 2026-10-28
requires: [ADR-0203, ADR-0209]
supersedes: []
superseded-by:
issues: ["#0181"]
---

# ADR-0211: erf と erfcx の x ≥ 0.46875 を Cody の有理近似で求める

## 背景

ADR-0209 の既定では、自己遮蔽の時間の大半が自由ガスの核の D の散乱源で、1C の 0 MWd/kg のピンセルの 1 スレッドのプロファイルの 39% が
erfcx（x ≥ 2.5 の 40 段の連分数）だった（#0181 の「問題」）。2A の燃焼のジョブは ADR-0209 の前の 18 分 58 秒から 1 時間 47 分 59 秒になった。
x < 2.5 では、erf を級数で（x = 2 で 29 回、2.4 で 35 回の割り算）、erfcx と erfc をその 1 − erf から求めていた。

erf と erfcx は xs の `special` と、その複製の lattice の `ufg::free_gas` にある（ステージのクレートは互いに依存しない。ADR-0203）。
lattice の自由ガスの核は erf と erfcx を呼ぶ。xs では `resonance/fluctuation.rs`（非分離共鳴の幅の揺らぎの積分）、`multigroup/free_gas.rs`
（自由ガスの核）が erfcx と erf を、`multigroup/outgoing.rs`（Maxwell 分布の累積 `gamma_p_3_2`）が erf を使う。xs の erfc を呼ぶのはテストだけ。

判断の軸:

- 正確さ: 倍精度の多倍長の参照との相対誤差（ulp）で、今の計算と同じかそれより小さい
- 速さ: 1 スレッドの `free_gas` の時間
- 保守: 依存は標準ライブラリだけ。係数を手で書くなら、出典と、係数を確かめた方法を残す

## 検討した案

精度は、x の 0.46875〜1.5、1.5〜2.5 の一様な乱数 2 万点ずつ、2.5〜4、4〜6 の 4 万点ずつ、6〜50 の 2 万点、50〜1e8 の対数で一様な 5,000 点と
区間の端で、mpmath 1.4.1 の 60 桁の erfc(x)·exp(x²) との差を、参照値の倍精度の ulp で測った（最大 / 二乗平均。x < 2.5 は最大だけ）。
時間は、同じ点を 200 周した 1 回あたりの時間。どれも手元の使い捨ての計測（Zen 5、`target-cpu=znver5`、1 コア、2026-10-28）で、時間は
ほかのジョブと CPU を共有した、時間を記録するジョブではない目安。

| erfcx の案                             | 0.47〜1.5 | 1.5〜2.5 | 2.5〜4    | 4〜6      | 6〜50     | 50〜1e8   |
| -------------------------------------- | --------- | -------- | --------- | --------- | --------- | --------- |
| 今の計算（級数、2.5 から連分数 40 段） | 131       | 16,444   | 2.53/0.71 | 2.88/0.85 | 2.97/0.78 | 2.59/0.80 |
| Cody の有理近似（CALERF のまま）       |           |          | 5.48/1.27 | 1.48/0.52 | 1.43/0.49 | 1.40/0.49 |
| Cody、4 以下を補償 Horner 法（採用）   | 2.05      | 2.36     | 1.83/0.53 | 1.48/0.52 | 1.43/0.49 | 1.40/0.49 |
| Faddeeva package の erfcx              |           |          | 2.98/0.81 | 3.98/1.05 | 2.36/0.59 | 2.68/0.69 |

erfcx の 1 回の時間は、2.5 ≤ x < 50 で、連分数 57〜69 ns、CALERF のまま 3.2 ns、Faddeeva 3.0 ns、採用した形 4.7 ns。0.46875 ≤ x < 2.5 で、
級数 37 ns、採用した形 8.0 ns。erf（0.46875〜2.5 の一様な乱数 4 万点、mpmath の erf との差）は、級数が最大 13.3 ulp（x = 2.42。掛ける
exp(−x²) に x² の丸めが x² 倍で効く）、34 ns で、1 − erfcx(x) exp(−x²) が最大 3.5 ulp、10 ns。erfc は、0.46875〜2.5 の 4 万点で 1 − 級数が
最大 16,712 ulp、0.46875〜26 の一様な乱数 3 万点で erfcx(x) exp(−x*x) が最大 510 ulp（x = 24.6）、exp(−x²) を下の分け方にすると最大 4.9 ulp。
負の x の erfcx（−26〜0 の 3 万点）は、2 exp(x*x) − erfcx(−x) が最大 502 ulp、分け方で最大 4.6 ulp。

- 案 A: 連分数の段の数を x で減らす。今の連分数と 200 段の連分数がビットで同じになる段は、x = 2.5 で 38、4 で 21、6 で 14。連分数を呼ぶ x の
  78% が 2.5 ≤ x < 6（#0181）なので、減らせる段は少ない。採らない
- 案 B: Cody の有理 Chebyshev 近似（W. J. Cody, "Rational Chebyshev approximations for the error function", Math. Comp. 23 (1969)
  631–637）。Netlib の specfun の CALERF（1990-03-19 版、パブリックドメイン）が倍精度の係数を持つ。erfcx を 0.46875 ≤ x ≤ 4 で 8 次 / 8 次の
  有理式、x > 4 で (1/√π − w R(w))/x（w = 1/x²、R は 5 次 / 5 次）で求める。有理式そのものの誤差は、mpmath で係数の 10 進の値のまま
  評価して、0.46875 ≤ x ≤ 4 で相対 1.2e−18（0.01 ulp）、4 ≤ x ≤ 4e4 で 5.9e−19。CALERF の Horner 法のままでは、2.5〜4 で丸めが積もって
  最大 5.5 ulp になり、連分数より粗い（FMA の Horner 法でも 3.9 ulp）
- 案 C: 案 B の 4 以下の分子と分母を、補償 Horner 法（Graillat, Langlois, Louvet, "Compensated Horner scheme", 2005）で求める。
  各段の積と和の丸めを誤差なしの変換（積は融合積和、和は Knuth の TwoSum）で取り出して足し戻すので、係数と x が正の多項式は丸め 1 回
  ほどの誤差になる。x > 4 は R の項が 1/√π に比べて小さく、丸めが薄まるので、FMA の Horner 法のままで連分数より細かい
- 案 D: Faddeeva package（S. G. Johnson、MIT ライセンス）の erfcx。x < 50 は y = 4/(4 + x) の 100 区間の 6 次の Chebyshev 多項式
  （係数 700 個）、x ≥ 50 は連分数の 5 段。速いが、4〜6 で最大 4.0 ulp と連分数より粗く、表が大きく、ライセンスの表示が要る。採らない
- 案 C の範囲を 2.5 ≤ x に限り、x < 2.5 の級数を残す: 初めにこの形で実装した（レビューの前のコミット dabd8ac8a）。級数の 1 − erf は erf が
  1 に近いほど桁落ちし、2.5 の直下で 16,444 ulp。自由ガスの核の遠い下方散乱は erfcx(u) − erfcx(x2) の打ち消し合いなので、この誤差が差に
  出うる（u と x2 が近いと、減少関数の差の符号が逆になりうる）。級数の割り算の数も連分数と同じほどで、遅いまま残る。採らない。採用した形も
  誤差が 1 ulp を超えうるので、隣り合う浮動小数点数の間では ulp の桁の増減が残る（区間の境の 0.46875 と 4 を含む。打ち消し合いの差に
  出る誤差は、どの x でも ulp の桁）

## 決定と理由

案 C を、CALERF と同じ 0.46875 から使う。

- erfcx: 0.46875 ≤ x ≤ 4 は Cody の 8 次 / 8 次の有理式を補償 Horner 法で、x > 4 は Cody の 1/x² の有理式を FMA の Horner 法で求め、
  x < 0.46875 は級数の erf から exp(x²)(1 − erf(x))（erf(x) < 0.5 なので 1 − erf は桁落ちしない。最大 0.89 ulp）。最大の誤差は x ≥ 2.5 で
  1.83 ulp 以下（連分数は 2.97 ulp）、0.46875 ≤ x < 2.5 で 2.36 ulp 以下（級数は 16,444 ulp）で、1 回の時間は 1/4.6〜1/12
- erf: x ≥ 0.46875 は 1 − erfcx(x) exp(−x²)（lattice は x ≥ 6 で 1）、x < 0.46875 は級数。最大 3.5 ulp（級数は 13.3 ulp）で、1/3.4 の時間
- erfc（xs）: x ≥ 0.46875 は erfcx(x) exp(−x²)、x < 0.46875 は 1 − 級数。exp(−x²) は CALERF と同じく exp(−s²) exp(−(x − s)(x + s))
  （s は x を 1/16 の倍数に切り捨てた値で、s² は丸めなしで表せる）に分け、x*x の丸めを指数に入れない。最大 4.9 ulp（x ≤ 26）
- 負の x の erfcx: erfc(−x) = 2 − erfc(x) から 2 exp(x²) − erfcx(−x)。exp(x²) は同じ分け方。最大 4.6 ulp
- erf は exp(−x*x) のまま（分けない）。x*x の丸めの効きは erfc(x) 倍に縮み、xs と lattice の erf のビットが同じになる
- 分子と分母は 1 つのループで進める。独立な 2 つの依存の列が並び、別々のループの 1.7 倍の速さ（13.7 → 7.9 ns、0.46875 ≤ x < 2.5）で、
  ビットは同じ
- 積和は `f64::mul_add` で書く。融合積和は 1 回の丸めなので、FMA の命令の無い CPU（ライブラリの fma）でも同じビットになる
  （`target-cpu` を付けずにコンパイルした計測でも、テストの点のビットが同じ）
- 係数は CALERF の倍精度の DATA 文の 10 進の値を丸めずにコピーし、上の mpmath の評価で確かめた。コードのコメントに出典と確かめ方を書く
- xs と lattice に同じコードを置く（ADR-0203 の複製の理由のまま）。xs の値も変わるので、ライブラリを作り直す。作り直すのは 1C と 2A の
  燃焼の再測定に要る examples/depletion（65 核種）と examples/depletion/p2（90 核種）だけで、ほかの例のライブラリは `bin/xs` が変わって
  changed になるが、結果が要るときに作り直す（ADR-0130）

## 確認

テストは xs の `special::tests` と lattice の `ufg::free_gas::tests` の両方にあり、点の表は xs の `special/erfcx_cases.rs` に置いて、2 つの
テストが include! で同じ表を読む。多倍長の参照値は `special/erfcx_cases.py`（mpmath 1.4.1、60 桁）が書き出す。

- `erfcx_is_within_two_and_a_half_ulp_of_multiprecision_values`: 区間の端と、上の計測で誤差の大きかった x を含む 41 点で、mpmath の
  60 桁の値との差が 2.5 ulp 以下（連分数は x = 36.8 で 2.97 ulp、級数は x = 2.475 で 16,444 ulp なので落ちる）
- `erf_is_within_four_ulp_of_multiprecision_values`: 10 点で 4 ulp 以下（級数は x = 2.42 で 13.3 ulp なので落ちる）
- `erfc_is_within_five_ulp_of_multiprecision_values`（xs だけ）: 11 点で 5 ulp 以下（exp(−x*x) では x = 24.6 で 510 ulp なので落ちる）
- `erfcx_of_negative_x_is_within_five_ulp_of_multiprecision_values`: 8 点で 5 ulp 以下（2 exp(x*x) − erfcx(−x) では x = −5.94 で 30 ulp）
- `erf_and_erfcx_of_nan_are_nan`: NaN は NaN を返す（lattice の erf は x ≥ 6 だけを 1 にする）
- ulp は、多倍長の値の側の隣の浮動小数点数までの間隔で数える（値が 2 の冪の直下のとき、上の間隔で数えると誤差を半分に数える）
- `erfcx_gives_the_same_bits_as_its_copy`: 40 点のビットの表に、2 つのクレートの erfcx が同じ値を返す。点は x = 0 と x ≥ 0.46875 で、
  四則と融合積和だけで求まる（libm の exp を通る x < 0.46875 は、libm の版でビットが変わりうるので置かない）
- 1 スレッドの `free_gas`（1C の 0 MWd/kg のピンセル `issue-0133/pin-1C-0.toml` に ADR-0209 の既定を明示したもの、手元、2026-10-28。
  ほかのジョブと CPU を共有した、時間を記録するジョブではない交互の 2 回）: main の 00d789a1f で 9.74〜9.75 秒、238ab21ff（erf を
  erfcx から求めるコミット。後の負の x と erfc の変更は核の呼ぶ値を変えない）で 4.84〜4.86 秒。k の標準出力は 6 桁で同じ（1.172387）

## 帰結

- ライブラリの値と lattice の自由ガスの核が変わる。examples/depletion と examples/depletion/p2 のライブラリを作り直し、1C と 2A の
  燃焼を再測定する（examples/depletion/verification.md）
- 自由ガスの散乱源の残りの手間は #0181 に書く
