---
status: open
labels: [lattice, core, maintenance]
opened: 2026-11-06
closed:
milestone:
blocked-by: []
resolution:
---

# lattice の CMFD の粗メッシュの固有値問題を xlc-cmfd に移す

## 問題

core の CMFD（ADR-0232）は、粗メッシュの固有値問題のブロック Gauss-Seidel のべき乗法を共通のクレート `xlc-cmfd` の `Coarse` に置いて使う。lattice の
`cmfd::Coarse`（crates/lattice/src/cmfd.rs の `solve_gauss_seidel`、1 行か 1 列の格子のブロック三重対角の逆反復 `solve_chain`。ADR-0158）は同じ
反復の複製を持つ（ADR-0003 は 2 つ目のステージが使う時点で共通のクレートに切り出すと決めた）。今移すと lattice の実行ファイルが変わり、`bin/lattice` を
deps にするステージ（15 の dvc.yaml の dvc.lock に約 380 行。うち examples/vera/p5 が 218）がすべて changed になるので、ADR-0232 は次に lattice の実行
ファイルが変わる作業に回した（`xlc-numeric` に置くと、lattice のコードを変えなくても xs と lattice の実行ファイルのバイトが変わったので、別のクレートにした）。

移すときに合わせるもの:

- 変化の見積もり（`Remaining`）: lattice は丸めの水準の扱いがなく（変化 0 だけを収束とみなす）、`xlc-cmfd` は相対 1e-14 以下を収束とみなす
- 負の解の確認: `xlc-cmfd` は生成率で割る前に確かめる。lattice の Gauss-Seidel の形は割った後に確かめる
- 対角のブロック: `xlc-cmfd` は逆行列を平らな配列に持ち、lattice はセルごとの LU 分解を持つ（結果の丸めの桁が変わりうる）
- 組の重み（`sizes`。対称性の代表のセルの数）: `xlc-cmfd` には無い

関係する ADR: なし

## 完了の条件

lattice の実行ファイルが変わる作業（lattice のコードの変更）と同じブランチで、lattice の粗メッシュの固有値問題を `xlc-cmfd` に移し、lattice の
CMFD のテストと VERA の格子計算の例の結果が許容の中で変わらないことを確かめる。

## 経緯

- 2026-11-06: ADR-0232 の案 S3 で、core だけが `xlc-cmfd` を使う形にした
