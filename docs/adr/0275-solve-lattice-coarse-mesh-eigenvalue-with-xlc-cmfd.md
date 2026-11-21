---
type: ADR
status: accepted
date: 2026-11-22
requires: [ADR-0003, ADR-0158, ADR-0232, ADR-0247]
supersedes: []
superseded-by:
issues: ["#0241", "#0328"]
---

# ADR-0275: lattice の CMFD の粗メッシュの固有値問題（2 次元の格子）を xlc-cmfd の Coarse で解き、反復と変化の見積もりを core と共通にし、対角のブロックの解き方は呼ぶ側が選ぶ

## 背景

core の CMFD（[ADR-0232](0232-accelerate-nem-with-node-cmfd.md)）は、粗メッシュの固有値問題のブロック Gauss-Seidel のべき乗法を共通のクレート
`xlc-cmfd` の `Coarse` に置いた。lattice の `cmfd::Coarse`（`solve_gauss_seidel`）は同じ反復の複製で、ADR-0232 の案 S3 は、lattice の実行ファイルが
次に変わる作業で lattice を同じクレートに移すとした（#0241。[ADR-0003](0003-extract-shared-computation-on-second-use.md)）。変化から残りの誤差を
見積もる `Remaining` は、xlc-cmfd、core の NEM の外側の反復、lattice（MOC の外側の反復、一般化随伴の反復、粗メッシュの問題）の 3 か所にあった。

2 つの形は次の 4 点で違った（#0241 の「移すときに合わせるもの」）:

- 変化の見積もり: lattice は変化 0 だけを収束とみなし、xlc-cmfd と core は相対 1e-14 以下の変化を収束とみなす
- 負の解の確認: lattice の Gauss-Seidel は生成率で割った後に、xlc-cmfd は割る前に確かめる
- 対角のブロック: lattice はセルごとの LU 分解（`xlc_numeric::lu::Lu`）で解き、xlc-cmfd は逆行列を平らな配列に持って掛ける
- 組の重み: lattice は対称性の代表のセルの生成率を組のセルの数（`sizes`）だけ数え、xlc-cmfd には重みが無い

lattice だけが使うものに、1 行か 1 列の格子のブロック三重対角の逆反復（`solve_chain`。[ADR-0158](0158-solve-one-row-coarse-mesh-directly.md)）と、
一般化随伴の粗メッシュの補正の固定源の問題（`solve_projected`。[ADR-0247](0247-generalized-perturbation-density-derivatives.md)）がある。
`solve_projected` は、2 次元の格子では、対角のブロックを使うブロック Gauss-Seidel の 1 回を GMRES の作用素にする。

判断の軸: 正確さ > 実行速度 > 保守（CLAUDE.md の「方針」）。lattice の実行ファイルは 2026-11-20 の変更（[ADR-0268](0268-solve-unlisted-fuel-nuclides-by-group-constants.md) など）で
既に変わっていて、ステージが changed になることは見送る理由にしない（[ADR-0258](0258-keep-dvc-and-deps-forms-and-carry-changed-as-rebuild-backlog.md)）。

## 検討した案

変化の見積もり:

- 案 R1（採用）: xlc-cmfd の形（相対 1e-14 以下の変化を収束とみなす）を公開し、lattice の MOC、一般化随伴、粗メッシュの問題と core の NEM が使う
- 案 R2: lattice の形（変化 0 だけ）に寄せる。許容が丸めの水準より小さいと、揺れる変化で縮小率が 1 を超えて止まらない（xlc-cmfd の注釈の理由）

負の解の確認:

- 案 N1（採用）: 生成率で割る前に確かめる（xlc-cmfd の形）
- 案 N2: 割った後に確かめる（lattice の形）。すべてのセルで負の解は、生成率も負なので割ると正に見える。対角が −0.2 の 1 群の 1 つのセルの問題で、
  lattice の形は k = −0.5 を返した（`cmfd::tests::a_coarse_solution_negative_in_every_cell_is_an_error` の RED、main 630459cba のコード、2026-11-21）

対角のブロック:

- 案 B1: セルごとの LU 分解の前進・後退代入（lattice の形。`xlc_numeric::lu::Lu` と同じ演算）を、lattice と core の両方に使う。分解はすべてのセルで
  1 つの平らな配列に持ち、行の入れ替えは右辺を読む並びにする。core の全炉心の BOC の固有値計算が 13.0 秒で、案 B2 の 10.0 秒より 3 割長い（Problem 9 の形の
  台帳とサイクル、ホウ素 1100 ppm、燃焼度 0、Xe なし、平らな初期値、`RAYON_NUM_THREADS=1`、手元の続けて流した目安。入力は examples/depletion/p9-timing/trial-0244/boc-1100.toml の
  コピー、2026-11-22）。外側の反復は 47 回で同じ。逆数を掛ける後退代入（13.1 秒）でも差は縮まらなかった（コミットしていない試行）。core の 2 群のブロックは
  条件数が 8 以下で（確認の「core のブロックの条件数」）、逆行列の誤差は LU と同じ桁なので、core では精度を得ずに時間だけを失う（いったん採り
  （e24a1bd4a）、マージの前の親のセッションの指摘（2026-11-22）で案 B3 に改めた）
- 案 B1′: 案 B1 を、セルごとの `Lu` の列で持つ。core の同じ問題が 14.0 秒で、平らな配列の 13.0 秒より長い
- 案 B2: 逆行列を平らな配列に持って掛ける（xlc-cmfd の形）を、lattice と core の両方に使う。逆行列を掛ける形は後退安定でなく、条件数の悪いブロックで
  誤差が LU の最大 5.4 倍になる（確認の「ブロックの誤差」）。lattice のブロックは、補正で対角が減ると条件数が 10⁶〜10¹⁰ になりうる
- 案 B3（採用）: `xlc_cmfd::Coarse` を作るときに、呼ぶ側が解き方（`BlockSolve::Lu` か `BlockSolve::Inverse`）を選ぶ。既定は置かない。反復の本体
  （Gauss-Seidel、変化の見積もり、負の解の確認、重み）は 1 つのまま

組の重み:

- 案 W1（採用）: xlc-cmfd の `Coarse` に生成率の重み（`Coarse::weighted`、既定は 1）を足す
- 案 W2: lattice が代表の核分裂の項を組の数で拡大して渡す。核分裂源の右辺にも掛かり、方程式が変わる

lattice だけが使う解き方の置き場:

- 案 C1（採用）: 1 行か 1 列の格子の逆反復（`Chain`）と `solve_projected` の GMRES は lattice に残す。`solve_projected` の 2 次元の格子の 1 回は
  xlc-cmfd の `Coarse::sweep` で行う
- 案 C2: `Chain` も xlc-cmfd に移す。core は 1 行か 1 列のノードの炉心を解かないので、ADR-0003 の「2 つ目のステージが使う時点」に当たらない

## 決定と理由

lattice の CMFD の粗メッシュの固有値問題は、2 次元の格子（と代表が 1 つの格子）では `xlc_cmfd::Coarse` を、代表の組のセルの数を重みにして解く
（案 W1）。1 行か 1 列の格子の逆反復は lattice に残す（案 C1）。共通の形は次のとおり:

- **変化の見積もり**: `xlc_cmfd::Remaining`（1 の桁で測った変化（1 の桁の k の変化と中性子束の相対変化）が 1e-14 以下なら収束とみなす）を、lattice の MOC の外側の反復、一般化随伴の反復、
  粗メッシュの問題（逆反復を含む）と、core の NEM の外側の反復が使う（案 R1）
- **負の解**: 1 回の後、生成率で割る前に、中性子束が正の有限であることを確かめる（案 N1）
- **対角のブロック**: `Coarse::new` と `Coarse::colored` の引数 `BlockSolve` で呼ぶ側が選ぶ（案 B3。既定は置かない）。lattice は `Lu`、core は `Inverse`。
  基準は、補正でブロックの条件数が大きくなりうるか: lattice の細かい群のブロックは対角が補正で減りうるので `Lu`。core の 2 群のブロックは、内側の面の項が
  D̃ の選び方（[ADR-0232](0232-accelerate-nem-with-node-cmfd.md)）で正、外周の面の項（外向きの正味の中性子流 / φ）が真空と鏡面反射で負にならず、
  条件数が小さい（確認の「core のブロックの条件数」で最大 9.0）ので `Inverse`。core で群の数を増やすか、ブロックの対角を減らす補正を入れるときは、条件数を再測定して選び直す
- **一般化随伴の 1 回**: `solve_projected` の 2 次元の格子の A⁻¹ の 1 回は、`xlc_cmfd::Coarse::sweep`（固定源を足したブロック Gauss-Seidel の 1 回）で行う

理由:

- 案 R1: 1e-14 は、k の変化（k は 1 の桁なので相対と同じ）と中性子束の相対変化の丸めの水準で、lattice の粗メッシュの許容（k 1e-10、中性子束 1e-8 から
  緩める。[ADR-0110](0110-solve-cmfd-as-finely-as-the-moc-change.md)）と MOC の許容より 4 桁以上小さい
- 案 N1: 案 N2 は、すべてのセルで負の解を見逃す
- 案 B3: 正確さと実行速度が衝突するのは条件数の悪いブロックだけで、lattice には LU の精度を、core には逆行列の速さを残せる（core の結果と時間は
  main と同じ。確認の「core」）。解き方が 2 つになる保守の費用は、どちらも失わないことと釣り合う。逆行列を掛ける形は
  模型の近似でなく、同じ方程式を丸めの桁で解く（CLAUDE.md の「速度のための近似は入力の項目として足し」に当たらない）。収束した NEM の解は粗メッシュの
  問題の解き方で変わらず、core の出力は main とバイトで同じ（確認の「core」）。LU を平らな配列に持つのは、案 B1′ より 7 % 短いため
- 案 W1: 重みは生成率の和にだけ掛かり、解の比の大きさ（生成率 1 の規格化）を決める。重み 1 の積は値を変えないので、重みを付けない core の解は変わらない
- 案 C1: ADR-0003 の基準。`solve_projected` の 1 回を xlc-cmfd に置くのは、対角のブロックと隣との結合を lattice に二重に持たないため

### 帰結

- lattice の粗メッシュの Gauss-Seidel の複製と、core と lattice の `Remaining` の複製がなくなる
- lattice と core の実行ファイルが変わり（core はバイトが変わるが、出力は同じ）、`bin/lattice` と `bin/core` を deps にするステージは changed になる（main 630459cba の時点で、
  `binary-forms/lattice.md5` を deps にする 210 と `binary-forms/core.md5` を deps にする 97 のステージは、どれも既に changed だった）
- core の固有値計算の結果と時間は変わらない（確認の「core」）
- lattice の 2 次元の格子の出力は、変わる前とビットで同じ（確認の「lattice」）。変わるのは一般化随伴の 1 回の丸め（k で割る代わりに 1/k を掛ける）だけ

### 確認

- ブロックの誤差: VERA 2A の細かい群の媒質（69 群。main 630459cba の lattice、`medium` の出力、2026-11-21）から V (diag(Σt + 漏れ) − Sᵀ) の形の
  ブロックを作り、吸収を 1〜10⁻⁶ 倍に、漏れを 0、ふつう（4D/h）、負（補正で対角が減る −2D/h）にした 21 個で、xlc_numeric の LU と xlc-cmfd の
  逆行列の形の演算を Python で同じ順に行った解を、残差を有理数で厳密に求めて 4 回直した解と比べた（コミット 56e4978b9 の trials/0241/blocks.py）。相対誤差の
  最大は、漏れのあるブロック（条件数 11〜12）で LU 5.9×10⁻¹⁶〜1.0×10⁻¹⁵、逆行列 8.1×10⁻¹⁶〜1.3×10⁻¹⁵、漏れが負のブロック（条件数 6.5×10⁶〜
  9.0×10⁹）で LU 3.5×10⁻¹⁴〜2.2×10⁻¹²、逆行列 1.0×10⁻¹³〜2.8×10⁻¹²（最大で 5.4 倍）、漏れのないブロック（条件数 1.1×10²〜8.6×10⁷）は
  どちらも同じ桁
- 下の比較の「前」は main 630459cba、「後」は 2096445e1 のビルド（`cargo build --release`）で、どちらも手元で続けて流した（2026-11-22）。入力は
  コミット 56e4978b9 の trials/0241/ と examples/depletion/p9-timing/trial-0241/、examples/depletion/trial-0241-gen/（例の入力のライブラリと出力の
  パスだけを書き換えたコピー）、比べたのは同じコミットの trials/0241/compare.py
- lattice: VERA 2A（UO2。8 分の 1 の対称性で代表のセル 55、組の重みが 1 でない）、2G（AIC の制御棒）、2B-vacuum（真空の辺。代表のセル 190）、
  3A の軸方向の平板（1 行の格子。逆反復）で、標準出力、ピン出力、少数群定数（3A は 58 のノードの出力）は前と後でビットで同じ。MOC の反復
  （25、24、25、182 回）と粗メッシュの反復の和（7,913、7,643、11,668、1,771 回）も同じ。粗メッシュの区間（trace の `coarse` の和）は 0.92 秒から
  0.86 秒、0.90 秒から 0.86 秒、4.70 秒から 4.34 秒、18.24 秒から 18.38 秒（3A の逆反復のコードは変えていない）
- 一般化随伴: examples/depletion/issue-0279/lattice-3x3-few-group.toml を核種 3 つにし、`direct` と `shielded` を除いた入力（3 × 3 の格子、
  [sensitivity.few-group]）で、感度の 60 個の値の相対差の最大は、前と後（143b3bf3f と 2096445e1 のどちらも）で 8.0×10⁻¹⁴〜8.2×10⁻¹⁴（許容 1e-10）。
  同じビルドを続けて 2 回流しても値は揺れ、前の 2 回の間で 1.3×10⁻¹⁴（揺れの原因は確かめていない。候補は並列の和の順。#0328）
- core: 全炉心の BOC の固有値計算の標準出力と、反復ごとの k と変化のログの行は、前と後で同じ。1 スレッドの時間は前 10.09、10.08 秒、後 10.10、
  10.25 秒。Problem 9 の形の燃焼の時間の内訳の入力（examples/depletion/p9-timing/core-p9.toml のコピー、32 スレッド）は、燃焼の履歴（msgpack）と
  標準出力がバイトで同じで、時間は前 65.0 秒、後 65.1 秒。案 B1 の形（143b3bf3f）では、BOC が 10.2 秒から 13.2 秒、燃焼が 65.2 秒から 70.6 秒で、
  点の k の相対差の最大 2.1×10⁻¹⁴、数密度 3.1×10⁻¹¹ だった
- core のブロックの条件数: lattice の少数群定数（VERA 2A、2G、2B-vacuum の集合体と、3A の軸方向の 58 のノード。反射体を含む）から、core の
  `assemble` の形のブロック V (diag(Σt − Σgg) − 群の間の散乱) を、径方向 21.5、10.75、5.375 cm（`radial-nodes` が 1、2、4）、高さ 4、8、50 cm の
  ノードで作り、内側の面の項（D̂ = 0 の有限差分の 2D/h）を足した形と足さない形（鏡面反射の外周に当たる）で、2 ノルムの条件数の最大は 9.0
  （上部の水の反射体、5.375 cm 角、高さ 4 cm。足さない形は 4.6）（コミット 66f7d6036 の trials/0241/core_blocks.py）。面の項は高速群の対角を
  多く増やすので、細かいノードほど条件数が上がるが、項が大きいほど対角は D の比に近づく。core の例の入力の群の境界は、2 群（[1e7, 0.625, 1e-5]）だけ。
  core の実行の中で組み立てたブロック（D̂ を含む）の条件数は測っていない
- テスト: `cmfd::tests::a_coarse_solution_negative_in_every_cell_is_an_error`、`cmfd::tests::representatives_count_their_orbits_in_the_production`
  （lattice）、`tests::weights_count_each_cell_as_many_times_as_it_stands_for`、`tests::unit_weights_leave_the_solution_bit_for_bit`、
  `tests::weights_must_be_positive_and_one_per_cell`、`tests::a_sweep_solves_each_block_by_the_lu_substitution`、
  `tests::an_inverse_block_solve_multiplies_by_the_inverse_from_the_lu_factors`、`tests::changes_at_rounding_level_count_as_converged`（xlc-cmfd）
