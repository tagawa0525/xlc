---
status: closed
labels: [build, dvc]
opened: 2026-10-14
closed: 2026-10-14
milestone: M4
blocked-by: []
resolution: ADR-0163
---

# faer を依存に持つクレートを足すと、build の段階の 3 つの実行ファイルがすべて変わる

## 問題

燃焼のソルバのクレート `xlc-depletion`（ADR-0162。faer に依存し、どの段階のクレートもまだ依存しない）を足した
ブランチ feat/cram-depletion-solver（main の 321efe1c に rebase したもの、2026-10-14）で、直下の dvc.yaml の `build` を実行すると、
3 つの実行ファイルの md5 がすべて変わった（xs: e77257b3 → 563bfe80、lattice: 85a3c566 → 9846663c、core: fe2d8761 → 6cf69103）。
`crates/depletion` を外し、`Cargo.toml` と `Cargo.lock` を main に戻すと、main の dvc.lock と同じ md5 に戻る。このままマージすると、
実行ファイルを deps に持つ計算の段階がすべて古くなる。

原因の候補: `scripts/build-binaries.sh` はワークスペースの全体を 1 回の `cargo build --release` で作り（ADR-0159 の案 T1）、Cargo は
共有する依存の feature をワークスペースの全体で統合する。`cargo tree --workspace -e features` では、faer が加わると num-complex に
`bytemuck`（faer-traits から）が足され、libm が入る。xs は num-complex に依存するが、lattice と core は依存しないので、lattice と
core が変わる経路は確かめていない。

試したこと（ブランチ fix/build-only-stage-binaries、捨てた。2026-10-14）:

- `cargo build --release --locked -p xlc-xs -p xlc-lattice -p xlc-core`（3 つを 1 回の cargo で）: `xlc-depletion` を含めたまま、
  3 つの md5 が main の dvc.lock と一致した。main（321efe1c）でも一致する
- 3 つを別々の cargo（`-p` を 1 つずつ）: main で xs は一致、lattice（85a3c566 → 90e3f50c）と core（fe2d8761 → 7aa410c7）が変わる。
  xs の依存の feature が、今も lattice と core のバイトを決めている

捨てた理由（レビュー）: ADR-0159 の案 T1 は、手で作る `cargo build --release` の `target/release` と同じ実行ファイルにするために
ワークスペースの全体をビルドすると決めており、`-p` にするとその性質（perf などで `target/release` を使う。docs/environment.md、
ADR-0160）が崩れる。ADR-0159 を置き換える決定になる。3 つを 1 回の cargo にしても、lattice と core のどちらかが新しい依存を
持てば残りの 2 つも変わりうる（別々の cargo の測定のとおり結びついている）ので、根本の解決ではない。別々の cargo は結びつきを
切るが、lattice と core の段階を一度再実行することになる。

関係する ADR: [ADR-0163](../adr/0163-build-stage-binaries-as-one-package-set.md)

## 完了の条件

- 次のどれかを ADR で決めてある: (a) ワークスペースの全体のビルドのまま、faer を足したときに段階を再実行する、(b) 3 つを
  1 回の cargo で作る（ADR-0159 を置き換え、`target/release` との違いを docs/environment.md に書く）、(c) 3 つを別々の cargo で作る
  （lattice と core の段階を一度再実行する）
- feat/cram-depletion-solver をマージする前に、その決定で `build` を実行し、`bin/` の md5 と再実行する段階が分かっている

## 経緯

- 2026-10-14: feat/cram-depletion-solver の rebase のあと、`dvc repro -s build` で見つけた。`-p` の 2 つの形を試して捨てた
- 2026-10-14: (b)（3 つを 1 回の cargo）に決めた（ADR-0163）。T1 の「`target/release` と同じ実行ファイル」は置き換え、perf のときの作り方を docs/environment.md に書いた
