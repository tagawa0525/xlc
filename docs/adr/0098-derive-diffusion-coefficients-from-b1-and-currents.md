---
status: accepted
date: 2026-10-09
supersedes: [ADR-0044]
superseded-by:
issues: ["#0021"]
---

# ADR-0098: 少数群の拡散係数は、反射境界の集合体では B1 の基本モードの中性子流から、非反射の格子ではまず outflow の輸送断面積から求める

## 背景

拡散係数は M0 で D = 1/(3Σt) とした（[ADR-0044](0044-collapse-to-two-groups.md)。`crates/lattice/src/fewgroup.rs` の `homogenize`。輸送補正なし）。
D = 1/(3Σt) を確かめる比較はない。その後、ライブラリは散乱行列を P0〜P3 の Legendre のモーメントで持ち（[ADR-0076](0076-store-legendre-scattering-moments.md)）、
MOC は非等方散乱を解き、散乱の次数が 1 以上なら、解く途中で領域ごとに中性子束の高次のモーメント（実の球面調和関数の
符号つきの成分）を持つ（[ADR-0075](0075-solve-moc-with-anisotropic-scattering.md)。いまの MOC の解の出力には含めない）。M3 では、ADF を GET で作り
（[ADR-0094](0094-define-adf-by-one-node-nem.md)）、反射体の定数を作る（[ADR-0095](0095-make-reflector-constants-in-2d-with-normalized-adf.md)）。

- GET の ADF は、定数を作った体系では D の誤差も吸収する。D の良し悪しは、定数を別の体系（単一集合体の定数を
  カラーセットや炉心で使う）で使ったときに出る
- 反射体（水）では、水素の前方散乱で輸送断面積が全断面積よりかなり小さく、D = 1/(3Σt) は D を小さく見積もる
- 外周が鏡面反射の単一集合体では、正味の中性子流が 0 で、MOC の解から D を直接は求められない
- 原理で決まり、手で合わせる係数を持たないこと

## 検討した案

- D = 1/(3Σt) のまま（[ADR-0044](0044-collapse-to-two-groups.md)）: 散乱の非等方性を入れない。水の多い領域と反射体で D を小さく見積もる
- outflow（out-scatter）の輸送補正 Σtr,g = Σt,g − Σ_g' Σs1(g→g')、D = 1/(3Σtr): 広く使われる（MPACT 4.4 Theory
  Manual、式 (5.24)〜(5.26)。手元の写し）。軽い核種で誤差が出る。MPACT は軽い核種に inflow の補正、水素に neutron leakage
  conservation を使う（同 5.3.3 節。B. C. Yee, E. W. Larsen, B. Kochunas, PHYSOR 2016 を引く。読んでいない）
- inflow（in-scatter）の輸送補正 Σtr,g = Σt,g − Σ_g' Σs1(g'→g) φ1_g' / φ1_g（中性子流の1次のモーメント φ1 で重みづけ）:
  より厳密だが、中性子流のスペクトル φ1 が要る。STREAM は減速材と近似のバックリングの 1D の P_N の計算で φ1 を求める
  （D. Lee, S. Choi, K. S. Smith, "On the Diffusion Coefficient of Two-step Method for LWR analysis," Trans. KNS Autumn
  Meeting (2015)、式 (10)）
- 臨界スペクトル（B1 または P1）の計算: 均質化した多群の定数で、臨界のバックリングを探し、B1 の方程式の中性子流 J_g から
  D_g = J_g / (B φ_g) とし、中性子束で少数群に縮約する（同上、式 (1)〜(5)）。B1 の方程式は P1 の散乱行列を含み、in-scatter
  の効果も入る。Polaris は臨界スペクトルの計算を既定で行い（Mode = SEARCH）、P1 と B1 を選べる（SCALE 6.3 の Polaris の
  マニュアル、3.2.7.5 節 option<CRITSPEC>）。CASMO-5 と HELIOS の扱いは一次資料で確かめていない
- MOC の正味の中性子流と中性子束の勾配から D = −J / ∇φ とする: 外周が鏡面反射の集合体では J = 0 で定まらない。勾配が
  小さいところで不安定になる。GET の ADF とは別の自由度を、局所の比で決めることになる

比較の例: STREAM と PARCS の2段階の計算を MCNP6 と比べ、outflow + B1 は k −39 pcm、集合体出力の差の RMS 0.76%、最大
1.79%、inflow + P1 は −58 pcm、0.62%、1.29%、inflow + 基本モード（D = 1/(3Σtr)）は −66 pcm、0.57%、1.28%、outflow +
基本モードは −122 pcm、2.36%、4.70% だった（同上の表 I）。反射体は燃料と並べた 2×1 のカラーセットから、臨界スペクトルの
計算なしに作った（同上 2.2 節）。詳しい比較は S. Choi, K. S. Smith, H. Kim, T. Tak, D. Lee, J. Nucl. Sci. Technol. 54(6),
705–715 (2017)、https://doi.org/10.1080/00223131.2017.1299648 （要旨だけ読んだ）。

## 決定と理由

少数群は、[ADR-0044](0044-collapse-to-two-groups.md) のとおり境界 0.625 eV の2群とし、断面積は領域・群ごとの反応率を保存するように均質化・縮約する
（ADR-0044 から引き継ぐ。比べた案の記録はない）。

拡散係数は、格子の境界によって次のように求める。どちらも、ライブラリの P1 の散乱行列と MOC の解から、手で合わせる係数
なしに決まる。

- 外周が鏡面反射の単一集合体（燃料の定数）: 集合体を細かい群のまま中性子束で均質化し（全断面積、P0 と P1 の散乱行列、
  νΣf、χ）、B1 の方程式で臨界のバックリングを探し、D_g = J_g / (B φ_g) を B1 の中性子束で少数群に縮約する。正味の
  中性子流が 0 の格子では、漏れのある体系での中性子流のスペクトルを、均質な無限媒質の基本モードとして求めるほかないため。
  B1 は P1 の散乱行列（in-scatter を含む）と、P1 より高い次数の効果（関数 γ(B²/Σt²)。B → 0 で 1）を含み、P1 より原理に近い。B² < 0（k∞ < 1）では γ の解析接続を使う
- 外周が鏡面反射でも核分裂のない格子（水だけ、核分裂性でない集合体）では、臨界のバックリングが定まらない。このときは
  B² → 0 の極限の B1 の中性子流の式（P1 の式、γ = 1）3 Σt,g d_g − 3 Σ_g' Σs1(g'→g) d_g' = φ_g を、均質化した MOC の
  中性子束 φ で解き、D_g = d_g / φ_g を MOC の中性子束で縮約する。臨界の場合と同じ式の極限で、in-scatter を含み、
  中性子流のスペクトルを決めるのに使えるのは格子自身の中性子束だけだからである（outflow の補正にすると、同じ鏡面反射の
  格子で、核分裂の有無で in-scatter の扱いが変わる）。いまの MOC は核分裂のない格子を解かない（固有値計算だけ）ので、
  この場合は定数を作る経路ができたときに使われる
- B1 には P1 の散乱行列が要る。ライブラリに P1 のモーメントがない（P0 だけで処理した）ときはエラーにする
- 非反射の格子（反射体、カラーセット）: まず outflow の輸送断面積 Σtr,g = Σt,g − Σ_g' Σs1(g→g') を細かい群で領域ごとに
  作り、Σtr を中性子束と体積で重みづけして少数群に均質化し（真空の領域は Σtr = 0 で、中性子束と体積は重みに入る）、
  D_G = 1/(3⟨Σtr⟩_G) とする。D を直接に重みづけると、ピンの隙間のような真空の領域で D が無限大になる。中性子流のスペクトルを使わずに定まり、
  [ADR-0044](0044-collapse-to-two-groups.md) の 1/(3Σt) より原理に近い。反射体では outflow の補正が漏れを大きく見積もる傾向がある（上の Lee らの論文）ので、
  inflow の補正に進めるかを #0021 で決める
- inflow の補正の候補は、MOC の中性子流の1次のモーメント φ1 を重みにする形である。φ1 は領域ごとの符号つきのベクトルの
  成分で、成分が打ち消す領域や、φ1_g が 0 に近い群では、比 φ1_g' / φ1_g と φ1 の重みが定まらない。反射体の隣では
  高速群と熱群の中性子流の向きが逆になる。定まる形の候補は、ノードの正味の中性子流（面の正味の中性子流から作るベクトル）の
  向きへの射影を重みにすることで、そのときも次を守る: MOC の散乱の次数が 1 以上であること（でなければエラー）。ノードの
  正味の中性子流がほぼ 0（対称な位置）なら、黙って別の定義に替えず、エラーにする。重みの定義は #0021 で、文献（STREAM は
  近似のバックリングでの 1D の P_N の計算で φ1 を求める）と比べて決める

少数群の断面積の縮約の重みは、いまのとおり MOC の中性子束とし、B1 のスペクトルで重みづけし直さない（臨界スペクトルで
断面積を重みづけし直すかは、k∞ の保存と GET の ADF との整合に効くので、カラーセットでの比較のあとに別に決める）。

この ADR は [ADR-0044](0044-collapse-to-two-groups.md) を置き換える（拡散係数の D = 1/(3Σt) をやめ、2群の境界と縮約は上に書き直して引き継ぐ）。

## 確認

- B1（`crates/lattice/src/b1.rs`）: `isotropic_scattering_gives_one_third_total_at_zero_buckling`（Σs1 = 0、B² = 0 で
  1/(3Σt)）、`isotropic_scattering_follows_gamma_at_finite_buckling`（Σs1 = 0 で D = 1/(3γΣt)、B² の正負とも）、
  `diagonal_p1_scattering_gives_the_transport_cross_section`（P1 が対角だけなら B² = 0 で 1/(3(Σt − Σs1))）、
  `one_group_critical_mode_matches_the_transport_dispersion_relation`（1群の等方散乱では B1 は厳密で、臨界の B² と D が
  輸送方程式の平面波の分散関係 c = x/arctan(x)、c = y/artanh(y) に 10⁻⁸ で一致）、`critical_buckling_sign_follows_k_inf`
  （|k − 1| < 10⁻¹⁰、B² の符号は k∞ − 1 の符号）、`modes_beyond_the_pole_count_as_too_negative`、
  `negative_buckling_without_a_bracket_is_an_error`、`nan_is_rejected`、`zero_buckling_current_solves_the_p1_equation`、
  `gamma_series_joins_the_closed_form`、`no_fission_is_an_error`
- 格子（`crates/lattice/src/fewgroup.rs`）: `b1_diffusion_matches_a_hand_computed_two_region_lattice`（χ の違う2つの
  核分裂性の領域と上方散乱のある格子を、実装の均質化を使わない手書きの2群の B1 の式と比べる。1群の D は B1 の中性子束で
  縮約し、MOC の中性子束で重みづけした値と 0.1% 以上違う）、`reflective_lattice_takes_d_from_the_b1_mode`、
  `reflective_lattice_without_fission_uses_the_p1_current_of_the_moc_flux`、
  `open_lattice_without_p1_scattering_uses_one_third_total`、`open_lattice_uses_the_outflow_transport_cross_section`、
  `outflow_counts_a_void_gap_with_its_flux`、`diffusion_needs_p1_moments`（P1 がなければ、どちらの境界でもエラー）。
  縮約が k∞ を保つことは `two_group_constants_reproduce_k_inf`、`one_group_collapse_preserves_k`（ADR-0044 から）
- VERA 1B と 2B（支持格子なし）の2群の D は、1/(3Σt) の 0.662、0.289 cm と 0.655、0.278 cm から、B1 で 1.419、0.469 cm と
  1.446、0.454 cm になった（臨界の B² は 2.82 × 10⁻³、2.85 × 10⁻³ /cm²。examples/vera/verification.md、63cb20d、
  2026-10-09）
- 単一集合体の定数を使う体系で、D の作り方（1/(3Σt)、outflow、B1）ごとに、core の NEM と同じ多群定数の xlc の
  MOC を比べること（4-2D、REF1-2D）は、まだしていない（「結果」）

## 結果

- 非反射の格子の均質化（反射体とカラーセット、[ADR-0095](0095-make-reflector-constants-in-2d-with-normalized-adf.md)）はまだないので、outflow の D は
  `diffusion_coefficients` にあるが、少数群定数の出力（`homogenize`、反射境界だけ）からはまだ使われない
- D の作り方の良し悪しは、4-2D と REF1-2D の xlc の MOC との比較（ADR-0094、ADR-0095 の閾値の案）で確かめる。
  それまで 1/(3Σt) と B1 のどちらが core に良いかは測っていない
- 反射体の D を inflow の補正にするかは #0021
- B1 の中性子流の式の Σt と P1 の散乱行列は、中性子流でなく中性子束で重みづけして均質化した近似である。外周が鏡面反射の
  MOC では正味の中性子流が 0 で、中性子流の重みが得られないため。この近似の大きさは、D の作り方の比較（4-2D）で
  ほかの誤差と一緒に見る
- 単一集合体の定数の ADF（GET、[ADR-0094](0094-define-adf-by-one-node-nem.md)）は平らな均質な解で決まり、D に依らない。D は、別の体系で使ったときの漏れの
  精度に効く。反射体の ADF は D に依るので、D を変えたら反射体のテーブルを作り直す
- 細部は、xlc の格子計算（MOC）の結果との比較を通じて段階的に詰め、新しい ADR か issue に記録する（2026-10-08 の
  ユーザーの方針）
