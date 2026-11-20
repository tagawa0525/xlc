---
type: ADR
status: accepted
date: 2026-11-20
requires: [ADR-0055]
supersedes: []
superseded-by:
issues: ["#0264"]
---

# ADR-0273: 荷電粒子のチャネルの透過因子を、転回点の外の Steed の CF2 と、障壁の下の Riccati の式の ln P の積分で求める

## 背景

4 つの版（ENDF/B-VII.0、VII.1、VIII.0、VIII.1）の評価をすべて xs で処理する（ユーザー 2026-11-20）。分離共鳴の R 行列の制限形
（LRF=7。[ADR-0055](0055-reconstruct-r-matrix-limited-subset.md)）のうち、荷電粒子の粒子対を持つ 5 評価（VII.1、VIII.0、VIII.1 の
Cl-35 の (n,p) の MT600、VIII.0 と VIII.1 の Ca-40 の (n,p) の MT600 と (n,α) の MT800）が「a charged particle is not supported」で
止まっていた（#0264）。荷電粒子のチャネルの透過因子は Coulomb の波動関数の P_l(η, ρ) = ρ / (F_l² + G_l²)（ENDF-102 2023 年版の
式 D.80〜D.84）で、障壁の下では e^{−2πη} の桁まで小さくなる。共鳴エネルギーでの評価の範囲は、Cl-35 の陽子で η 0.88〜3.2、ρ 0.82〜2.9、
Ca-40 の陽子で η 3.0〜17.9、ρ 0.18〜1.06、α で η 5.7〜9.4、ρ 2.3〜3.9（l は 0〜2。2026-11-20 に評価のファイルから求めた）。Ca-40 の
(n,p) は閾値 541.9 keV が分離共鳴の範囲（〜1.5 MeV）の中にあり、閾値に近づくと η → ∞ になる。

判断の軸:

- 正確さ: 近似に落とさない。障壁の深い下でもアンダーフローとオーバーフローを起こさず、倍精度の丸めの大きさで正しい値を出す
- 実行速度: 1 つのエネルギーで spin group とチャネルごとに何度も呼ぶ
- 保守: 依存を足さない。式の出典をたどれる

## 検討した案

- 案 A: NJOY2016 samm（`pspcou`）と同じ方法。ρ < 1.02 では級数（`coulx`。η > 10ρ かつ η > 5 では η の大きいときの漸近の式 `bigeta`、
  η ≥ 5 では漸近の式と Taylor 展開、許容 1e-6）、ほかは Steed の方法（`coulfg`）で、転回点の内側で G が 1e6 を超えると JWKB の近似。
  透過因子が 1e-35 以下なら 0 にし、共鳴の幅を換算幅にするときは 0 の透過因子を 1 にする（`betset`）。近似で、式の切り替えで値が
  跳ぶ: Ca-40 の (n,p) の NJOY の断面積は、η = 10ρ になる 822768 eV で 5.47e-9 b から 6.72e-4 b に不連続に変わり、α の l = 0 の JWKB の
  透過因子は厳密な値より 1.3〜2.0% 大きい（examples/depletion/verification.md の「版ごとのライブラリの分離共鳴の R 行列の制限形の荷電粒子の
  粒子対」）
- 案 B: Steed の CF1（F'/F）と CF2（H'/H）だけで F、G、F'、G' を求める。障壁の下では H'/H = (S + iP)/ρ の虚部が実部より e^{−2πη} の桁で
  小さく、CF2 の丸めの誤差（実部の大きさ × 機械精度）に埋もれる
- 案 C: 原点のまわりの級数で F と G を求める（`coulx` の一般の形）。G の級数は対数と digamma の項の打ち消しを持ち、ρη が大きいと
  桁落ちする。F と G を別々に持つと、G² が倍精度を超える（η ≳ 110）ので指数の桁を別に持つ必要がある
- 案 D: 転回点の外で CF2 から Y = H'/H を求め、Riccati の式 Y' = V − Y²（V = 2η/ρ + l(l+1)/ρ² − 1）を ρ の内向きに積分する。Wronskian
  F'G − FG' = 1 から ρY = S + iP なので、P を直接得る。P の代わりに ln P を積分すればアンダーフローしない（前の担当の案）
  - D1: 埋め込みの Runge-Kutta（Dormand-Prince の 5(4) 次）で刻みを決める。相対の許容 1e-12 で、η 3.3、ρ 0.8 の l = 0 に 336 ステップ、η 17.9、
    ρ 0.18 の l = 2 に 1867 ステップ、η 100、ρ 0.03 の l = 0 に 6026 ステップかかった（P の誤差は 1e-12 の桁。x = ln ρ で積分しても 324、1666、5448 ステップ。
    手元の試作、2026-11-20）。ln P の大きさ（〜2πη）に対して絶対の誤差を小さくする必要があり、5 次の方法では刻みが小さい
  - D2: x = ln ρ で、S（= ρ Re Y）と ln P の Taylor 級数を漸化式で求め（Riccati の式は 2 次の多項式なので係数の漸化式が閉じる）、末尾の係数から
    刻みを決める。20 次で同じ 3 点に 12、26、87 ステップで、P の相対の誤差は 1e-16〜7e-13（手元の試作、2026-11-20）
  - D3: WKB の解を引いた残りだけを積分する。残りは小さいが、転回点のすぐ内側で WKB の解の微分が発散するので、転回点の近くを別に積分する
    必要がある
- 案 E: 多倍長の計算（mpmath のような）を使う。依存が増え、1 つの値に 0.02〜35 秒かかった（mpmath 1.4.1、40 桁。35 秒は η 1000、ρ 0.003）

## 決定と理由

案 D2 にする（`crates/xs/src/resonance/coulomb.rs`）。

- 転回点 ρ_t = η + √(η² + l(l+1)) と 2 の大きい方を ρ₀ とし、ρ ≥ ρ₀ では CF2（DLMF 33.8.2 の連分数を修正 Lentz 法で）から P = ρ Im Y を
  直接求める。ρ₀ の下限の 2 は、転回点の小さいチャネル（l = 0、η = 0.01 で ρ_t = 0.02）で CF2 の項が増える（ρ = 0.02 で 2515 項、ρ = 2 で
  52 項）のを避けるため。CF2 の先頭の項を分子に置くと η が 1e5 を超えたとき初期値で桁があふれるので、分母の連分数に Lentz 法を当てる
- ρ < ρ₀ では、ρ₀ の Y から x = ln ρ で dS/dx = S − S² + P² + l(l+1) + 2ηρ − ρ²、d(ln P)/dx = 1 − 2S を内向きに積分する。案 D1 よりステップの数が
  1 桁少なく、誤差は ln P の丸めの大きさ（|ln P| × 2e-15 以内）で、D3 のような領域の分け方が要らない。次数は、Ca-40 の 0 K の再構成の
  時間が 12 次で 13.9 秒、16 次で 8.8 秒、20 次で 7.3 秒、28 次で 6.9 秒（手元、2026-11-20）だったので 20 にした
- 障壁の下（ρ < ρ_t）で S < 0 なら、内向きに S < 0 が続き（S = 0 では dS/dx = P² + l(l+1) + 2ηρ − ρ² > 0）、ln P は減り続ける。ln P が
  −746（最小の正の非正規化数 e^{−744.4} より下）を下回れば、P は倍精度で 0 なので、そこで 0 を返す。近似でなく、閾値のすぐ上
  （Ca-40 の陽子で η ≈ 3000/√(E − E_thr [eV])）の積分を打ち切れる
- η = Z_a Z_b α μ c / (ħ k)（式 D.80）。α は NJOY2016 phys の `finstri` と同じく EV、CLIGHT、HBAR から求める

荷電粒子の粒子対の扱い（[ADR-0055](0055-reconstruct-r-matrix-limited-subset.md) の範囲に足す）:

- ZA > 0、MA > 0 で MT が 600〜849（陽子、重陽子、三重陽子、³He、α の準位。ENDF-102 の 3.4.6）の粒子対を、非弾性散乱の準位と同じ出口の
  チャネルにする。PNT は 1 と既定の 0 を受け、−1 はエラーにする。ほかの MT の荷電粒子の粒子対はエラーのまま
- 荷電粒子のチャネルは入射にならないので、断面積は |U_ab|² = 4|X_ab|² で、位相シフト（Coulomb の位相を含む）を使わない
- シフト因子は Coulomb の波動関数の S_l（式 D.84）で、剛体球の式は使えない。4 つの版に SHF=1 の荷電粒子の対は無いので、SHF=1 はエラーにする
  （S は Riccati の積分の S なので、要るときは同じ関数から返せる）
- 共鳴の寄与は粒子対の MT（MT600、MT800）の MF3 に足し、MT103〜107 は部分の MT600〜849 の和で作り直す（今の xs の和の規則と NJOY2016
  RECONR の `anlyzd` と同じ）。5 評価の MF3 の MT103 と MT107 は分離共鳴の範囲で MT600〜649、MT800〜849 の和と同じ値で、共鳴を
  二重に数えない（2026-11-20 に評価のファイルの MF3 を読んで確かめた）。プロセスの側（`chain_reactions`）は作り直した MT103、MT107 を読む
- 共鳴エネルギーで透過因子が 0 にアンダーフローするチャネルに幅があれば、換算幅が定義できないのでエラーにする

NJOY との違い: NJOY RECONR（samm）の荷電粒子の透過因子は、案 A のとおり η ≥ 5 の障壁の下で近似の式になり、Ca-40 の MT600 と MT800 は NJOY と
桁で違う（Cl-35 は 2.9e-6 以下）。CLAUDE.md の「方針」のとおり NJOY を真似ず、共鳴エネルギーで障壁の下にあり η ≥ 5 の荷電粒子の幅があれば、
その数と MT を `tracing::warn!` で出す（`RMatrixLimitedXs::njoy_approximated_widths`）。差の大きさは examples/depletion/verification.md の
「版ごとのライブラリの分離共鳴の R 行列の制限形の荷電粒子の粒子対」。

## 確認

- `penetrability_matches_multiprecision_values`: l = 0〜4 と 19 の (η, ρ) で、mpmath の 80 桁の値と |ln P| × 2e-15（下限 2e-14）の相対で
  合う（参照値は crates/xs/src/resonance/coulomb/tests/cases.py が書き出す）。`penetrability_underflows_to_zero_far_below_the_barrier`、
  `zero_eta_is_the_hard_sphere_penetrability`
- `proton_channel_matches_single_level_breit_wigner_with_coulomb_penetrability`: 陽子のチャネルの単一準位の (n,p) が、η と ρ を重心系の
  エネルギーから別の道筋で求めた Breit-Wigner の式と合う。`endothermic_proton_channel_is_closed_below_its_threshold`、
  `rejects_widths_whose_coulomb_penetrability_underflows`、`rejects_shift_factor_for_charged_pairs`、
  `rejects_charged_pairs_outside_charged_particle_levels`、`rejects_charged_pairs_without_penetrability`、
  `lists_widths_whose_coulomb_penetrability_njoy_approximates`
- 5 評価の NJOY RECONR との比較と `xs nuclide`、既存の LRF=7 の 11 評価の出力が変わらないことは examples/depletion/verification.md の
  「版ごとのライブラリの分離共鳴の R 行列の制限形の荷電粒子の粒子対」

## 帰結

- 5 評価が `xs nuclide` を通る。2A のライブラリの放射化の核種に Cl-35 と Ca-40 を戻すのは、本番の版の決定の後（#0264）
- Ca-40 の 0 K の再構成の時間の 85% が透過因子で（perf、手元、2026-11-20）、再構成は 7.4 秒（NJOY RECONR は全体で 1 秒）。同じエネルギーで
  spin group をまたいで同じチャネル（粒子対、l、半径）の透過因子を何度も求めているので、速くするならエネルギーごとに再利用する
- Ca-40 の (n,p)、(n,α) と、それを通る全断面積と弾性散乱は NJOY と違う。NJOY の処理したライブラリ（OpenMC など）と比べるときの差の原因になる
