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

単体テストの置き方が 2 つ混ざっている。2026-10-12 の main（14793ea7）では、`crates/*/src` でテストを持つファイルの
うち 64 個がファイルの末尾に `mod tests { … }`（または `pub(crate) mod tests { … }`）を書き、11 個が
`#[cfg(test)] mod tests;` で `<モジュール>/tests.rs` に分けている（例: crates/lattice/src/branches/tests.rs、crates/format/src/table/tests.rs）。

ファイルの中に置いたテストは、大きなファイルで全体の半分前後を占める。

| ファイル                        | 行数 | テストの行数（割合） |
| ------------------------------- | ---- | -------------------- |
| crates/lattice/src/fewgroup.rs  | 3668 | 2056（56%）          |
| crates/lattice/src/moc.rs       | 3378 | 1889（56%）          |
| crates/core/src/nem.rs          | 3028 | 2035（67%）          |
| crates/lattice/src/shielding.rs | 3019 | 1710（57%）          |
| crates/lattice/src/input.rs     | 2153 | 1015（47%）          |

行数は 3668 から 57 まで切れ目なく並び、分ける・分けないのしきい値を置く自然な位置がない（テストを中に持つ 64 個の
ファイルで、テストのモジュールの `#[cfg(test)]` より前の行数の中央値は 308.5）。

判断の軸:

- 読む量: 本体を読むときにテストを読まずに済むか。LLM のエージェントはファイルを読んで理解するので、混ざったテストは
  そのまま読む量になる
- パスでの区別: grep、diff、`git log --stat` で、テストか本体かがパスだけで分かるか。TDD の RED と GREEN のコミット
  （テストと実装を別にコミットする）も、触ったファイルで確かめられる
- 機械的に守れるか: 人やエージェントの判断に任せず、検査で止められるか

## 検討した案

- 案 A: 今のまま、置き方を決めない。2 つの書き方が残り、大きなファイルは大きいまま
- 案 B: しきい値（例えば 1000 行。テストを中に持つ 64 個のうち超えるのは 13 個）を超えるファイルだけ分ける。変更は少ないが、しきい値に根拠がなく、
  ファイルが育ってしきい値を超えるたびに移す作業が出る。2 つの書き方が残る
- 案 C: テストを持つすべてのモジュールで `<モジュール>/tests.rs` に分け、本番のファイルの中のテストのコードを検査で禁じる。
  `#[cfg(test)] mod tests;` の子モジュールは `use super::*;` で非公開の項目に触れられるので、単体テストとしての扱いは
  ファイルの中に置いたときと同じ。clippy の `allow-unwrap-in-tests` などの設定（clippy.toml）も、`#[cfg(test)]` の付いた
  モジュールとして効く（既存の 11 個の tests.rs で効いている）
- 案 D: 案 C のテストのファイルを `src` の外に置き（例: crates/x/src/foo.rs に `#[cfg(test)] #[path = "../unit/foo.rs"] mod tests;` と書き、crates/x/unit/foo.rs に置く）、
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
  テストが使う項目（`sample`、`sample_pins` は crates/format/src/table/tests.rs が使う）は、`pub(crate) mod tests;` の中に置き、
  呼び出し側は `crate::fewgroup::tests::sample` のように `tests` を通して引く
- 本番のファイル（`crates/*/src` の `tests.rs` 以外のファイル）に残る `#[cfg(test)]` は、`mod tests;` と `pub(crate) mod tests;`
  の宣言（`main.rs` では `#[path]` を付けたもの）と、パスが `tests` のモジュールを通る `use`（ほかの
  モジュールのテストの補助を引くもの。2026-10-12 には crates/xs/src/endf/mod.rs の `pub(crate) use mf5::tests as mf5_tests;`、
  `pub(crate) use tape::tests::rec;`、crates/xs/src/resonance/mod.rs と crates/xs/src/multigroup/mod.rs の同じ形の
  `use`）だけにする。クレートの入口（`lib.rs`、`main.rs`）の `#![cfg_attr(test, allow(clippy::disallowed_macros))]`（テストのビルドで `assert!` を許す。Cargo.toml の
  lint の注記）は、テストのコードでなく lint の設定なので残す
- 本番のファイルが、上で残すもの以外のテストのコード（`#[test]`、`mod tests { … }`、条件に `test` を含む `cfg`、
  `cfg_attr`、`cfg!`）を持たないことを、`run-checks.sh` で走る検査で確かめる。検査の置き場所、照合の方法、残すものを
  減らせるか（補助の `use` を `tests.rs` に移すなど）は、実装のブランチ（#0089）で検査のテストとともに決める

案 D を採らないのは、古くなった段階を再実行する費用（テストを変えたときだけ）のために、モジュールの規則から外れた
`#[path]` をすべてのモジュールに求めるのは釣り合わないからである。案 C で `#[path]` を使うのは、`lib.rs` とぶつかる
`main.rs` の 1 か所だけである。

しきい値を置かないのは、根拠のある値がなく、すべてで分ければ検査が「本番のファイルにテストがない」の 1 つで済むから
である。小さなファイルではファイルが 1 つ増えるが、パスで区別できる利点は大きさによらない。

## 確認

- 検査が、本番のファイルの中のテストのコードを見つけて失敗し、移したあとの main で通ること
- 移す前後で、release の実行ファイルが変わらないこと（変わるなら、`#[cfg(test)]` の項目を消して後ろの本体の行番号が
  ずれたことによる位置の情報だけであること。2026-10-12 には、テストだけの項目の後ろに本体が続き、行がずれるのは
  crates/numeric/src/pins/node.rs だけ）と、実行されるテストの集まりが変わらないこと。確認の方法は実装のブランチ
  （#0089）で決め、結果はこの節に書く

## 帰結

- 実装（64 個のファイルの移動と検査）は、この ADR のあとの別のブランチで行う（#0089）。テストの移動と検査の追加を別のコミットにする
- 2026-10-12 に入口でテストを持つのは、xs の `main.rs`（`src/main/tests.rs` になる）と crates/docs の `lib.rs`
  （`src/tests.rs` になる）である
- ADR と issue には、テストの場所をファイルのパスで書いたものがある（例: ADR-0132 の crates/core/src/nem.rs）。
  採用済みの ADR は直せない（docs/adr/README.md の「運用」で、直せるのは誤字、リンク、`status` と `superseded-by` だけ）ので、決めた時点の
  場所の記録として残る。proposed の ADR と開いた issue のパスは、移すブランチで `tests.rs` に直す
- `examples/*/dvc.yaml` の段階は deps に `crates/*/src` のディレクトリを持つので、テストを移すと deps の md5 が変わり、
  xs、lattice、core の段階がほぼすべて古くなる。移すブランチで `dvc repro -P` を Slurm で 1 回流し、`dvc.lock` を
  コミットする（ADR-0123 の「実行して `dvc.lock` をコミットする」のとおり）。出力の md5 が移す前と変わったら、移動でなく
  並列の和の順序（#0060）などの原因を調べる。流す時期は、CPU の使用の計画を持つ実装の担当が決める
- ADR-0123（案 H）はテストの変更で段階が古くなることを許したが、この移動はすべての段階を一度に古くするので、古いまま
  残さない。古い上流の段階が main に残ると、ほかのブランチが /dvc-relock で数秒の
  下流の段階だけを再実行することができず、上流から実行することになる。`dvc commit` で deps だけを取り込むこともしない。
  `dvc.lock` が、記録した deps のコードで実行していない出力を指すことになる（ADR-0130 の案 L と同じ形）
