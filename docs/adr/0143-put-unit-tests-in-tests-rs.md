---
status: proposed
date: 2026-10-12
requires: [ADR-0123, ADR-0130]
supersedes: []
superseded-by:
issues: ["#0089"]
---

# ADR-0143: 単体テストを各モジュールの tests.rs に置き、本番のファイルにテストのコードを書かない

## 背景

単体テストの置き方が 2 つ混ざっている。2026-10-12 の main（1842033d）では、`crates/*/src` でテストを持つファイルの
うち 64 個がファイルの末尾に `mod tests { … }`（または `pub(crate) mod tests { … }`）を書き、11 個が
`#[cfg(test)] mod tests;` で `<モジュール>/tests.rs` に分けている（例: crates/lattice/src/branches/tests.rs、crates/format/src/table/tests.rs）。

ファイルの中に置いたテストは、大きなファイルで本体より長い。

| ファイル                        | 行数 | テストの行数（割合） |
| ------------------------------- | ---- | -------------------- |
| crates/lattice/src/fewgroup.rs  | 3668 | 2056（56%）          |
| crates/lattice/src/moc.rs       | 3378 | 1889（56%）          |
| crates/core/src/nem.rs          | 3028 | 2035（67%）          |
| crates/lattice/src/shielding.rs | 3019 | 1710（57%）          |
| crates/lattice/src/input.rs     | 2153 | 1015（47%）          |

行数は 3668 から数百まで切れ目なく並び、分ける・分けないのしきい値を置く自然な位置がない（テストを中に持つ 64 個の
ファイルで、テストのモジュールの `#[cfg(test)]` より前の行数の中央値は 308.5）。

判断の軸:

- 読む量: 本体を読むときにテストを読まずに済むか。LLM のエージェントはファイルを読んで理解するので、混ざったテストは
  そのまま読む量になる
- パスでの区別: grep、diff、`git log --stat` で、テストか本体かがパスだけで分かるか。TDD の RED と GREEN のコミット
  （テストと実装を別にコミットする）も、触ったファイルで確かめられる
- 機械的に守れるか: 人やエージェントの判断に任せず、検査で止められるか

## 検討した案

- 案 A: 今のまま、置き方を決めない。2 つの書き方が残り、大きなファイルは大きいまま
- 案 B: しきい値（例えば 1000 行。超えるのは 13 個）を超えるファイルだけ分ける。変更は少ないが、しきい値に根拠がなく、
  ファイルが育ってしきい値を超えるたびに移す作業が出る。2 つの書き方が残る
- 案 C: テストを持つすべてのモジュールで `<モジュール>/tests.rs` に分け、本番のファイルの中のテストのコードを検査で禁じる。
  `#[cfg(test)] mod tests;` の子モジュールは `use super::*;` で非公開の項目に触れられるので、単体テストとしての扱いは
  ファイルの中に置いたときと同じ。clippy の `allow-unwrap-in-tests` などの設定（clippy.toml）も、`#[cfg(test)]` の付いた
  モジュールとして効く（既存の 11 個の tests.rs で効いている）
- 案 D: 案 C のテストのファイルを `src` の外に置き（例: `#[cfg(test)] #[path = "../../unit/foo.rs"] mod tests;`）、
  DVC の deps（`crates/*/src`）から外す。テストを変えても段階が古くならない（ADR-0123 が「`src` の外に出すのは別の
  決定」として残した点）。代わりに、すべてのモジュールに `#[path]` が要り、Rust のモジュールの規則（ファイルの場所が
  モジュールの階層を表す）から外れる。テストの変更で段階が古くなることは、ADR-0123 が許している（案 H）

## 決定と理由

案 C を採る。

- テストは `#[cfg(test)] mod tests;`（ほかのモジュールのテストが使う補助を持つなら `pub(crate) mod tests;`）で宣言し、
  中身は Rust のモジュールの規則が探す `tests.rs` に置く。`foo.rs` のテストは `foo/tests.rs`、`lib.rs` と `mod.rs` の
  テストは同じディレクトリの `tests.rs`。`main.rs` のテストは `#[path = "main/tests.rs"]` で `src/main/tests.rs` に置く
  （xs、lattice、core は `src` に `lib.rs` と `main.rs` の両方を持ち、どちらも `mod tests;` だと同じ `src/tests.rs` を
  指してぶつかる）
- テストだけが使う項目（試料を作る関数、テストだけが呼ぶメソッド）も `tests.rs` に置く。メソッドは `tests.rs` の中の
  `impl` に書ける（固有の `impl` は同じクレートのどのモジュールにも置ける）。2026-10-12 には、crates/format/src/fewgroup.rs
  の `sample`、`sample_pins`、crates/numeric/src/pins/node.rs の `modes_at`、`flux` がこれにあたる。ほかのモジュールの
  テストが使う項目（`sample` は crates/format/src/table/tests.rs が使う）は、`pub(crate) mod tests;` の中に置き、
  呼び出し側は `crate::fewgroup::tests::sample` のように `tests` を通して引く
- 本番のファイルに残る `#[cfg(test)]` は、`mod tests;` の宣言と、パスが `tests` のモジュールを通る `use`（ほかの
  モジュールのテストの補助を引くもの。2026-10-12 には crates/xs/src/endf/mod.rs の `pub(crate) use mf5::tests as mf5_tests;`、
  `pub(crate) use tape::tests::rec;`、crates/xs/src/resonance/mod.rs と crates/xs/src/multigroup/mod.rs の同じ形の
  `use`）だけにする。クレートの入口（`lib.rs`、`main.rs`）の `#![cfg_attr(test, allow(clippy::disallowed_macros))]`（テストのビルドで `assert!` を許す。Cargo.toml の
  lint の注記）は、テストのコードでなく lint の設定なので残す
- crates/docs の一貫性の検査（#0013）に、`crates/*/src` の `tests.rs` 以外のファイルが、上で残すもの以外のテストの
  コード（`#[test]`、`mod tests { … }`、条件に `test` を含む `cfg`、`cfg_attr`、`cfg!`）を持たないことを確かめるテストを
  足す。照合の方法と、残すものを減らせるか（補助の `use` を `tests.rs` に移すなど）は、実装のブランチ（#0089）で
  検査のテストとともに決める

案 D を採らないのは、古くなった段階の走らせ直しは ADR-0123 で許した費用で、そのためにモジュールの規則から外れた
置き方をすべてのテストに求めるのは釣り合わないからである。

しきい値を置かないのは、根拠のある値がなく、すべてで分ければ検査が「本番のファイルにテストがない」の 1 つで済むから
である。小さなファイルではファイルが 1 つ増えるが、パスで区別できる利点は大きさによらない。

## 確認

- crates/docs の検査が、本番のファイルの中の `mod tests { … }` とテストだけの項目を見つけて失敗すること（検査のテストで、
  そうしたファイルを与えて確かめる）。移したあとの main で検査が通ること
- 移すのはテストのコードだけで、`#[cfg(test)]` の項目は release のビルドに入らない。ファイルの末尾の `mod tests { … }`
  と、後ろに本番の項目が続かないテストだけの項目（crates/format/src/fewgroup.rs の `sample`、`sample_pins`）だけを移す
  コミットでは、本体の行番号（`tracing` の位置の情報に入る）も変わらないので、移す前後で
  `cargo build --release` の xs、lattice、core の実行ファイルがビット単位で同じことを確かめる。実行ファイルには
  worktree の絶対パスが入る（examples/m0/dvc.yaml の注記。ADR-0123）ので、前後の 2 つのコミットは同じ worktree で
  checkout してビルドする。比べる前に、移す前のコミットを、間に `cargo clean` を挟んで 2 回ビルドし、実行ファイルが一致することを確かめる（同じ
  パスでのビルドし直しがビット一致するかは、ADR-0123 の案 I で確かめていない）。一致しなければ、比べ方を実装の
  ブランチで決め直す
- テストだけの項目のうち、本体の途中にあるもの（2026-10-12 には crates/numeric/src/pins/node.rs の `modes_at`、`flux`。
  後ろに本番の `gradient`、`average` が続き、lattice と core の実行ファイルに入る）を移すと、後ろの本体の行番号がずれ、
  実行ファイルの位置の情報だけが変わる。これは別のコミットにし、そのコミットの本番のファイルの差分が `#[cfg(test)]` の
  項目を消す行だけであることを、差分で確かめる。release のビルドに入らない項目を消しただけなので、出力は変わらない。
  この確認のために DVC の段階を走らせることはしない（段階は、下の「帰結」のとおり次の `dvc repro` で走る）
- 移す前後で、`cargo test --workspace` が出す `passed` の件数の和が同じこと

## 帰結

- 実装（64 個のファイルの移動と検査）は、この ADR のあとの別のブランチで行う（#0089）。テストの移動と検査の追加を別のコミットにする
- 2026-10-12 に入口でテストを持つのは、xs の `main.rs`（`src/main/tests.rs` になる）と crates/docs の `lib.rs`
  （`src/tests.rs` になる）である
- ADR と issue には、テストの場所をファイルのパスで書いたものがある（例: ADR-0132 の crates/core/src/nem.rs）。
  採用済みの ADR は直せない（docs/adr/README.md の「運用」で、直せるのは誤字、リンク、`status` と `superseded-by` だけ）ので、決めた時点の
  場所の記録として残る。proposed の ADR と開いた issue のパスは、移すブランチで `tests.rs` に直す
- `examples/*/dvc.yaml` の段階は deps に `crates/*/src` のディレクトリを持つので、テストを移すと、出力が変わらなくても
  deps の md5 が変わり、段階が古くなる。ADR-0123（案 H）のとおり、これを許し、次の `dvc repro` で走らせる。
  `dvc commit` で deps だけを取り込むことはしない（ADR-0130 の案 L と同じく、`dvc.lock` が、記録した deps から
  作っていない出力を指す）
