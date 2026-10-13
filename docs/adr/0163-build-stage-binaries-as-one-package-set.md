---
status: accepted
date: 2026-10-14
requires: [ADR-0159, ADR-0160]
supersedes: []
superseded-by:
issues: ["#0114"]
---

# ADR-0163: build の段階は xs、lattice、core の 3 つのパッケージを 1 回の cargo で作る

## 背景

直下の dvc.yaml の `build`（[ADR-0159](0159-depend-stages-on-built-binaries.md)）は、`scripts/build-binaries.sh` でワークスペースの全体を
`cargo build --release --locked`（`target/dvc`）し、`bin/xs`、`bin/lattice`、`bin/core` を作る。ADR-0159 の案 T1 は、全体を
ビルドする理由を「`-p` で選ぶと機能の統合が変わりうるので、手で `cargo build --release` したものと同じ実行ファイルにするため」と
した。計算の段階は、この実行ファイル（strip したもの。[ADR-0160](0160-strip-stage-binaries.md)）のバイトが変わったときだけ走る。

Cargo は、1 回のビルドで選んだパッケージの全体で、共有する依存の feature を統合する。燃焼のソルバのクレート `xlc-depletion`
（faer に依存し、どの段階もまだ依存しない）を足したブランチ feat/cram-depletion-solver（main の 5c21dfed の上、2026-10-14）で
全体をビルドすると、3 つの実行ファイルの md5 がすべて変わった（xs: e77257b3 → 563bfe80、lattice: 85a3c566 → 9846663c、
core: fe2d8761 → 6cf69103）。crates/depletion を外し、Cargo.toml と Cargo.lock を main に戻すと元の md5 に戻る（#0114）。
原因の候補は feature の統合で、`cargo tree --workspace -e features` では faer-traits が num-complex に `bytemuck` を足し、libm が
入る。lattice と core は num-complex に依存しないので、この 2 つが変わる経路は確かめていない。

- 実行ファイルのバイトは、その実行ファイルの依存の木だけで決まること（ワークスペースに別のクレートを足しても変わらない）

## 検討した案

- 案 A: ワークスペースの全体をビルドしたまま、依存を足すたびに段階を走らせ直す（ADR-0159 の T1）。実行ファイルの依存の外の
  クレートが実行ファイルを変え、そのたびに計算の段階がすべて古くなる
- 案 B: 3 つのパッケージを 1 回の cargo で `-p xlc-xs -p xlc-lattice -p xlc-core` としてビルドする。feature の統合は 3 つとその
  依存の木の中に閉じる。feat/cram-depletion-solver で、crates/depletion を含めたまま 3 つの md5 が元と一致した（採用）
- 案 C: 3 つを別々の cargo（`-p` を 1 つずつ）でビルドする。3 つの間の統合もなくなるが、main（5c21dfed）で lattice
  （85a3c566 → 90e3f50c）と core（fe2d8761 → 7aa410c7）の md5 が変わり（xs は変わらない）、lattice と core の段階を一度
  走らせ直すことになる

## 決定と理由

`scripts/build-binaries.sh` は `cargo build --release --locked -p xlc-xs -p xlc-lattice -p xlc-core --target-dir target/dvc` で
3 つのパッケージだけを 1 回の cargo でビルドする（案 B）。実行ファイルの feature の統合を、3 つの依存の木だけで決め、
ワークスペースに別のクレートを足しても変わらないようにするため。3 つのどれかが新しい依存を持てば（lattice か core が
`xlc-depletion` に依存するなど）、そのとき実行ファイルが変わるのは正当な変化として扱う。案 C は今の実行ファイルのバイトを
変え、段階を走らせ直すことになるので採らない。

ADR-0159 の案 T1 のうち、ワークスペースの全体をビルドして「手で `cargo build --release` したものと同じ実行ファイルにする」
部分を、これで置き換える。`target/dvc` に置き、3 つのパッケージを毎回 `cargo clean -p` してからビルドすること、`bin/` に置く
こと、strip（ADR-0160）は変えない。

## 確認

main（0ceb6fe3、2026-10-14）の上のこのブランチで、変更の前（全体のビルド）と後（3 つの `-p`）の `bin/` の md5 が一致する
（xs e77257b3、lattice 03b6048d、core fe2d8761）。main の dvc.lock の lattice（85a3c566）は 0ceb6fe3 の前のもので、0ceb6fe3 が
lattice を変えたあと `build` を走らせ直していない（そのマージの本文のとおり、lattice の段階はあとで測り直す）。
crates/depletion を含むブランチでも、3 つの `-p` のビルドで md5 が変わらないことを確かめる（feat/cram-depletion-solver）。

## 帰結

- 手で `cargo build --release`（ワークスペースの全体）した `target/release` の実行ファイルは、`bin/` と feature の統合が違いうる。
  perf などで `target/release` を使うとき（ADR-0160）、同じ機械語にしたければ `cargo build --release -p xlc-xs -p xlc-lattice -p xlc-core`
  で作る（docs/environment.md）
- 3 つの間の feature の統合は残るので、3 つのどれかの依存の変更は、ほかの 2 つの実行ファイルも変えうる（案 C の測定）
