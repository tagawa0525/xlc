---
type: ADR
status: proposed
date: 2026-11-21
requires: [ADR-0003, ADR-0158, ADR-0232, ADR-0247]
supersedes: []
superseded-by:
issues: ["#0241"]
---

# ADR-0275: lattice の CMFD の粗メッシュの固有値問題（2 次元の格子）を xlc-cmfd の Coarse で解き、変化の見積もりと対角のブロックの解き方を lattice、core と共通にする

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

- 案 B1（採用）: セルごとの LU 分解の前進・後退代入（lattice の形）を、xlc-cmfd でも使う
- 案 B2: 逆行列を平らな配列に持って掛ける（xlc-cmfd の形）

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

- **変化の見積もり**: `xlc_cmfd::Remaining`（相対 1e-14 以下の変化を収束とみなす）を、lattice の MOC の外側の反復、一般化随伴の反復、
  粗メッシュの問題（逆反復を含む）と、core の NEM の外側の反復が使う（案 R1）
- **負の解**: 1 回の後、生成率で割る前に、中性子束が正の有限であることを確かめる（案 N1）
- **対角のブロック**: セルごとの LU 分解の前進・後退代入で解く（案 B1）。core の結果も丸めの桁で変わる
- **一般化随伴の 1 回**: `solve_projected` の 2 次元の格子の A⁻¹ の 1 回は、`xlc_cmfd::Coarse::sweep`（固定源を足したブロック Gauss-Seidel の 1 回）で行う

理由:

- 案 R1: 1e-14 は、k の変化（k は 1 の桁なので相対と同じ）と中性子束の相対変化の丸めの水準で、lattice の粗メッシュの許容（k 1e-10、中性子束 1e-8 から
  緩める。[ADR-0110](0110-solve-cmfd-as-finely-as-the-moc-change.md)）と MOC の許容より 4 桁以上小さい
- 案 N1: 案 N2 は、すべてのセルで負の解を見逃す
- 案 B1: 逆行列を掛ける形は後退安定でなく、条件数の悪いブロックで誤差が大きくなりうる（確認の「ブロックの誤差」）
- 案 W1: 重みは生成率の和にだけ掛かり、解の比の大きさ（生成率 1 の規格化）を決める。重み 1 の積は値を変えないので、重みを付けない core の解は変わらない
- 案 C1: ADR-0003 の基準。`solve_projected` の 1 回を xlc-cmfd に置くのは、対角のブロックと隣との結合を lattice に二重に持たないため

### 帰結

- lattice の粗メッシュの Gauss-Seidel の複製と、core と lattice の `Remaining` の複製がなくなる
- lattice と core の実行ファイルが変わり、`bin/lattice` と `bin/core` を deps にするステージは changed になる

### 確認

- ブロックの誤差: VERA 2A の細かい群の媒質（69 群。main 630459cba の lattice、`medium` の出力、2026-11-21）から V (diag(Σt + 漏れ) − Sᵀ) の形の
  ブロックを作り、吸収を 1〜10⁻⁶ 倍に、漏れを 0、ふつう（4D/h）、負（補正で対角が減る −2D/h）にした 21 個で、xlc_numeric の LU と xlc-cmfd の
  逆行列の形の演算を Python で同じ順に行った解を、残差を有理数で厳密に求めて 4 回直した解と比べた（手元の試行 trials/0241/blocks.py。コミットしていない）。相対誤差の
  最大は、漏れのあるブロック（条件数 11〜12）で LU 5.9×10⁻¹⁶〜1.0×10⁻¹⁵、逆行列 8.1×10⁻¹⁶〜1.3×10⁻¹⁵、漏れが負のブロック（条件数 6.5×10⁶〜
  9.0×10⁹）で LU 3.5×10⁻¹⁴〜2.2×10⁻¹²、逆行列 1.0×10⁻¹³〜2.8×10⁻¹²（最大で 5.4 倍）、漏れのないブロック（条件数 1.1×10²〜8.6×10⁷）は
  どちらも同じ桁
