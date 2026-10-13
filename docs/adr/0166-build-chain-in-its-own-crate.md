---
status: accepted
date: 2026-10-14
requires: [ADR-0001, ADR-0003, ADR-0011, ADR-0013, ADR-0034, ADR-0159, ADR-0163]
supersedes: [ADR-0162]
superseded-by:
issues: []
---

# ADR-0166: 燃焼鎖を xlc-chain で作り、鎖の型を xlc-depletion に置いて MessagePack で渡す

## 背景

全核種の燃焼鎖（[ADR-0034](0034-build-full-depletion-chain.md)）を作るコード、鎖のデータファイルの形式と置き場、それを作る段階を決める。

[ADR-0162](0162-add-depletion-crate.md) は、ENDF を読むのは xs だけとし、燃焼鎖も xs が作り、段階の間をファイルで渡る型は `xlc-format` に置くとした。
同じ ADR は、燃焼のコードの変更で、燃焼を使わない xs の実行ファイルと段階（60 を超える）が変わらないことを軸にした。
DVC の計算の段階は、呼ぶ実行ファイル（strip したもの）のバイトが変わると走る（[ADR-0159](0159-depend-stages-on-built-binaries.md)）。

測ったこと（このブランチ、main の 60d8e217 の作業ツリー、2026-10-14。コミットしていない試し）: `xlc-format` に、serde の
構造体と読み書きの関数だけの鎖のモジュール（xs、lattice、core のどれも使わない）を足して `scripts/build-binaries.sh` で
ビルドすると、`bin/xs`、`bin/lattice`、`bin/core` の md5 がすべて変わった（xs: e77257b3 → f84383de、lattice: 03b6048d → 260f3d16、
core: fe2d8761 → 0c2b88d6）。モジュールを外すと元に戻る。依存するクレートに使わないモジュールを足すだけで、実行ファイルが
変わる（原因は確かめていない。codegen unit の分け方の候補）。xs のクレート（`xlc-xs` のライブラリ）に鎖のコードを足しても、
同じことが起きうる。

- 段階が古くなる範囲: 鎖のコードの変更で、xs、lattice、core の実行ファイルが変わらないこと
- 依存の向き: lattice と core（「進め方」5、7）は鎖を読み、燃焼行列を作る。ENDF は読まない
- 規則の数: [ADR-0162](0162-add-depletion-crate.md) の依存の一覧を、必要な分だけ変える

## 検討した案

鎖を作るコード:

- 案 B1: `xlc-xs` のライブラリのモジュールと、`xs` のサブコマンド: 鎖のコードの変更で `bin/xs` が変わり、xs の段階が
  すべて古くなる（上の測定と同じ理由。`xlc-xs` は直接は測っていない）
- 案 B2: `xlc-xs` のパッケージの別の実行ファイル: `bin/xs` はライブラリをリンクするので、案 B1 と同じ
- 案 B3: ENDF の読み口を共通のクレート（`xlc-endf`）に移し、鎖のクレートと xs がそれに依存する: コードを移すので `bin/xs` が
  変わり、xs の段階を一度すべて走らせ直す。ENDF の読み口が段階の境界から外に出る
- 案 B4: 新しいクレート `xlc-chain`（crates/chain、実行ファイル `chain`）を xs の段階の一部とし、`xlc-xs` のライブラリの
  ENDF の読み口（`endf::Tape` など、公開の API）に依存する（採用）。`xlc-xs` を変えないので `bin/xs` は変わらない。MF8〜10 と
  nfy の読み口は `xlc-chain` に置く

鎖の型（ファイルの形式）の置き場:

- 案 T1: `xlc-format`（[ADR-0162](0162-add-depletion-crate.md) の規則）: 上の測定のとおり、3 つの実行ファイルがすべて変わる
- 案 T2: `xlc-chain`: lattice と core は段階のクレートに依存できない
- 案 T3: `xlc-depletion`（採用）。鎖を読んで燃焼行列を作る lattice と core は `xlc-depletion` に依存する。燃焼行列の組み立て
  （鎖と反応率から行列）も同じクレートに置ける。いまは xs、lattice、core のどれも `xlc-depletion` に依存しないので、
  実行ファイルは変わらない。`xlc-depletion` は MessagePack の厳密な読み書き（[ADR-0013](0013-read-files-strictly.md)）に `xlc-format` の `msgpack` を使う

実行ファイルのビルド:

- `scripts/build-binaries.sh` の 1 回の cargo に `-p xlc-chain` を足す: `xlc-chain` は `xlc-depletion`（faer）に依存するので、
  3 つの実行ファイルの feature の統合に faer の依存が入り、3 つが変わる（[ADR-0163](0163-build-stage-binaries-as-one-package-set.md) の背景の測定と同じ）
- 3 つの cargo のあとに、`-p xlc-chain` を別の cargo でビルドする（採用）。3 つの feature の統合は変わらない

## 決定と理由

案 B4、T3 と、別の cargo でのビルドにする。[ADR-0162](0162-add-depletion-crate.md) を置き換え、まだ有効な部分を含めて次のとおりにする。

- 段階のクレート（xs、lattice、core、後の冷却計算）は互いに依存しない。依存してよいのは、共通のクレート（`xlc-format`、
  `xlc-numeric`、`xlc-depletion`）だけとする。段階の間の依存をファイル形式だけにする（[ADR-0001](0001-split-into-three-stage-crates.md)）ことを、クレートの依存でも
  守るため
- xs の段階は `xlc-xs` と `xlc-chain` の 2 つのクレートからなる。`xlc-chain` は `xlc-xs` のライブラリ（ENDF の読み口）と
  共通のクレートに依存してよく、ほかのクレートは `xlc-chain` に依存しない。ENDF を読むのは xs の段階のクレートだけとする
  （燃焼鎖も xs の段階が作る。[ADR-0034](0034-build-full-depletion-chain.md)）
- 段階の間を渡るデータファイルの型と MessagePack の厳密な読み書きは `xlc-format` に置く。ただし燃焼鎖の型は `xlc-depletion`
  に置く（`xlc-format` に足すと、鎖を使わない実行ファイルが変わるため）。燃焼のソルバ（[ADR-0161](0161-solve-depletion-with-ipf-cram.md)）と、鎖と反応率から燃焼行列を
  組み立てるコードも `xlc-depletion` に置く。`xlc-depletion` は `xlc-format` に依存してよい
- 鎖のファイル: フィールド名付きの map 形式の MessagePack（[ADR-0011](0011-write-data-files-in-messagepack.md)）。`format` は `"xlc-chain"`。読むときは、核種の添字の
  範囲、名前の形、核種の並び（(Z, A, 状態) の昇順）を確かめる
- 段階: examples/depletion/dvc.yaml の `chain` が `bin/chain chain.toml` で `chain.msgpack`（DVC の出力。git で追跡しない）を作る。
  比べる OpenMC の鎖は、同じ dvc.yaml の `openmc-chain` が `Chain.from_endf` で作る
- ビルド: `scripts/build-binaries.sh` は、3 つの実行ファイルを 1 回の cargo で作ったあと（ADR-0163）、`bin/chain` を別の
  cargo で作る

## 確認

- 鎖の型を足したこのブランチ（feat/depletion-chain-type、main の 60d8e217 の上、2026-10-14）で `dvc repro -s build` を走らせた
  `bin/xs`、`bin/lattice`、`bin/core` の md5（e77257b3、03b6048d、fe2d8761）が、main の dvc.lock と同じ
- `xlc-chain` を足し、`bin/chain` を別の cargo で作っても 3 つの md5 が変わらないことは、`xlc-chain` を足すブランチ
  （feat/depletion-chain）で確かめた（同じ 3 つの md5。2026-10-14）
- 鎖の読み書き: crates/depletion/src/chain/tests.rs（`round_trips_bit_exactly`、`written_as_a_map_with_the_format_name`、
  `rejects_another_format`、`rejects_nested_structs_written_as_arrays`、`rejects_an_omitted_decay` ほか）
- 依存の向き: 各クレートの `Cargo.toml`（`xlc-depletion` は `xlc-format` に依存し、xs、lattice、core は `xlc-depletion` に
  依存しない）。自動の検査はない

## 帰結

- `xlc-chain` と、`bin/chain` を作る build の段階、examples/depletion/dvc.yaml の段階は、全核種の鎖の規則を実装するブランチ
  （feat/depletion-chain）で足す。このブランチは鎖の型だけを足す
- lattice と core が `xlc-depletion` に依存するとき（「進め方」5、7）、その実行ファイルは変わる（正当な変化）
- `xlc-chain` は `xlc-xs` の ENDF の読み口の公開の API に依存するので、その API を変えると `xlc-chain` も直す
- 依存するクレートに使わないモジュールを足すと実行ファイルが変わる原因は確かめていない。`xlc-format` に型を足すたびに
  xs の段階が古くなるので、原因と対策（codegen unit の固定など）は #0117 で調べる
- 共通のクレートを足すとき、段階のクレートの依存を変えるときは、この ADR を置き換えて一覧を直す
