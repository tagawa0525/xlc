---
status: superseded
date: 2026-10-14
requires: [ADR-0001, ADR-0002, ADR-0003, ADR-0007, ADR-0025, ADR-0034, ADR-0159, ADR-0160]
supersedes: [ADR-0002]
superseded-by: ADR-0166
issues: []
---

# ADR-0162: 燃焼のソルバを共通のクレート xlc-depletion に置き、段階のクレートが依存してよい共通のクレートに加える

## 背景

[ADR-0025](0025-deplete-with-cram.md) と [ADR-0003](0003-extract-shared-computation-on-second-use.md) は、燃焼のソルバ（CRAM）を共通のクレートに置くとした（lattice の燃焼と core のミクロ燃焼が使い、
冷却計算を輸送計算なしで単独で動かすため）。[ADR-0002](0002-keep-stage-crates-independent.md) は、段階のクレート（xs、lattice、core）が依存してよいのを
`xlc-format` と `xlc-numeric` だけとした。燃焼のソルバは疎な複素 LU に faer（[ADR-0007](0007-use-faer-for-linear-algebra.md)）を使い、faer はまだどのクレートの
依存にもない。置き場として、`xlc-numeric` のモジュールか、新しい共通のクレートかを決める。

- 段階が古くなる範囲: DVC の計算の段階は、呼ぶ実行ファイル（strip したもの）のバイトが変わったときだけ走る
  （[ADR-0159](0159-depend-stages-on-built-binaries.md)、[ADR-0160](0160-strip-stage-binaries.md)）。xs の段階は 60 を超え、走らせ直しに時間がかかる。燃焼のコードの変更で、燃焼を使わない xs の
  実行ファイルが変わらないこと
- 依存の広がり: 燃焼を使わない段階と、冷却計算に、要らない依存（faer とその依存の gemm、pulp など、輸送と NEM のコード）を
  持ち込まないこと
- 規則の数: 段階が依存してよいクレートの一覧（[ADR-0002](0002-keep-stage-crates-independent.md)）を、必要がなければ増やさないこと

## 検討した案

- 案 A: `xlc-numeric` のモジュール `depletion` にする。[ADR-0002](0002-keep-stage-crates-independent.md) の一覧を変えない。xs は `xlc-numeric` に依存するので、燃焼の
  コードを足したり変えたりすると、`xlc-numeric` の codegen unit の分け方とインライン化が変わり、xs が使う関数の機械語の並びが
  変わりうる（release の設定に LTO はなく、codegen unit は既定の 16。変わるかは確かめていない）。変われば xs のすべての段階が
  古くなる。xs のビルドとテストにも faer とその依存が入る
- 案 B: 新しい共通のクレート `xlc-depletion`（crates/depletion）にする。xs の実行ファイルは燃焼のコードを含まないので、燃焼の
  コードを変えても xs の段階は古くならない（`build` は走るが、xs のバイトは変わらない）。冷却計算は `xlc-format`（チェーンと
  テーブルを読む）と `xlc-depletion` に依存し、`xlc-numeric` の輸送と NEM のコードを持たずに済む。[ADR-0002](0002-keep-stage-crates-independent.md) の一覧を変える。
  ただし、依存を共有するクレート（num-complex、num-traits など）の版と feature はワークスペースで統合されるので、faer を
  足したり版を上げたりして Cargo.lock の共有の依存が変わると、xs の実行ファイルも変わりうる。守れるのは燃焼のコードの変更までで、
  依存の変更は守れない（採用）
- 案 C: `xlc-numeric` の任意の機能（Cargo の feature）で faer と燃焼のモジュールを切り替える。ワークスペースをまとめてビルド
  する `build`（`cargo build --release`）では `xlc-numeric` の feature が統合されて xs にも燃焼のモジュールが入るので、燃焼の
  コードの変更については案 A と同じになる
- 案 D: lattice に置き、core が使うときに切り出す（[ADR-0003](0003-extract-shared-computation-on-second-use.md) の「2 つ目の段階が使う時点で」）。[ADR-0003](0003-extract-shared-computation-on-second-use.md) 自身が燃焼のソルバを
  M4 で共通にする予定に挙げ、[ADR-0025](0025-deplete-with-cram.md) が共通のクレートに置くと決めてある。冷却計算は lattice に依存できない

## 決定と理由

[ADR-0002](0002-keep-stage-crates-independent.md) を置き換え、まだ有効な部分を含めて次のとおりにする。

- 段階のクレートは互いに依存しない。後で冷却計算をクレートにするときは、段階のクレートとして扱う。依存してよいのは、データファイルの形式のクレート（`xlc-format`）と、共通の計算のクレート
  （`xlc-numeric` と `xlc-depletion`）だけとする。段階の間の依存をファイル形式だけにする（[ADR-0001](0001-split-into-three-stage-crates.md)）ことを、クレートの依存でも
  守るため
- 燃焼のソルバは `xlc-depletion`（crates/depletion）に置く。いまは faer と数値の型だけに依存し、`xlc-format` と `xlc-numeric` には
  依存しない（受け渡すのは燃焼行列と数密度だけ）。段階の間を渡る燃焼の型（燃焼点の状態など）は、[ADR-0003](0003-extract-shared-computation-on-second-use.md) のとおり分岐計算の
  流し方を決めるときに決め、ファイルで渡す型は `xlc-format` に置く
- 段階の間を渡るデータファイルの型（多群ライブラリ、群構造、少数群定数など）と MessagePack の厳密な読み書きは `xlc-format` に
  置く（M3 の初めに移した）
- ENDF を読むのは xs だけとし、燃焼チェーンも xs が ENDF の崩壊・収率データから作り、データファイルに書く（[ADR-0034](0034-build-full-depletion-chain.md)）

案 B にするのは、燃焼のコードの変更で、使わない xs の実行ファイルと段階が変わらないことを、クレートの境界で保てるからである。
一覧に 1 つ足す手間は、xs の段階を再実行する費用より小さい。

## 確認

各クレートの `Cargo.toml` の依存で確かめる（段階のクレートが別の段階のクレートを依存に持たない）。自動の検査はない。
crates/depletion を足すブランチで、直下の dvc.yaml の `build` を実行して `bin/` の実行ファイルが変わらないことを確かめる
（Cargo.lock の共有の依存が変わらないことを含む）。lattice と core が依存するブランチでも同じく、xs が変わらないことを確かめる。

crates/depletion を足したブランチ（feat/cram-depletion-solver、main の 321efe1c に rebase したもの、2026-10-14）では、Cargo.lock の
共有の依存の版は変わらなかったが、3 つの実行ファイルの md5 がすべて変わった。ワークスペースの全体のビルドで共有の依存の
feature が統合されるため（候補。num-complex に faer-traits の `bytemuck` が足される）で、案 B の但し書きの「依存の変更」に
あたる。扱いは #0114 で決める。

## 帰結

- `xlc-depletion` は CRAM の実装（計画の「進め方」4）で作る。lattice の燃焼（5）で lattice が、core の燃焼（7）で core が依存する
- 依存を足すときに実行ファイルが変わる問題は #0114
- 段階の間で共有する計算の扱いは [ADR-0003](0003-extract-shared-computation-on-second-use.md)、燃焼チェーンの作り方は [ADR-0034](0034-build-full-depletion-chain.md)
- 共通のクレートをさらに足すときは、この ADR を置き換えて一覧を直す
