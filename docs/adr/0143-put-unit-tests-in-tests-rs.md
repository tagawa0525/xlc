---
status: accepted
date: 2026-10-13
requires: [ADR-0130, ADR-0159, ADR-0160]
supersedes: []
superseded-by:
issues: ["#0089"]
---

# ADR-0143: 単体テストを各モジュールの tests.rs に置き、本番のファイルにテストのコードを書かない

## 背景

単体テストの置き方が 2 つ混ざっている。2026-10-12 の main（1842033d）では、`crates/*/src` でテストを持つファイルの
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

案 D を採らないのは、古くなった段階を走らせ直す費用（テストを変えたときだけ）のために、モジュールの規則から外れた
`#[path]` をすべてのモジュールに求めるのは釣り合わないからである。案 C で `#[path]` を使うのは、`lib.rs` とぶつかる
`main.rs` の 1 か所だけである。

しきい値を置かないのは、根拠のある値がなく、すべてで分ければ検査が「本番のファイルにテストがない」の 1 つで済むから
である。小さなファイルではファイルが 1 つ増えるが、パスで区別できる利点は大きさによらない。

## 確認

実装のブランチ（#0089、refactor/0089-move-tests。分岐点は main の 387743fc）で確認した（2026-10-13。6b39ad13 から rebase し直して、数と strip の比較をやり直した）。

- 検査: scripts/check-no-inline-tests.sh を、`.pre-commit-config.yaml` の hook `no-inline-tests` と、
  crates/docs/tests/consistency.rs の `production_sources_hold_no_test_code`（`cargo test` に含まれるので `run-checks.sh` でも走る。
  `run-checks.sh` は pre-commit の設定のうち design-records の hook しか走らせない）から呼ぶ。分岐点（64 個のファイルにテストがある）では
  失敗し、移したあとでは通る。照合は行単位の awk で、`#[test]`、`#[should_panic]`、`mod tests {`、test を含む `cfg` と `cfg_attr`、
  `cfg!(` を禁じ、決定で残すもの（`#[cfg(test)]` の次の行が `mod tests;`、`pub(crate) mod tests;`、`#[path = "…tests.rs"]` に続く
  `mod tests;`、パスが `tests` を通る `pub(crate) use`。入口の `#![cfg_attr(test, allow(clippy::disallowed_macros))]`）だけを許す。補助の
  `pub(crate) use` を `tests.rs` に移して残すものを減らすことは、しなかった（呼び出し側が `crate::endf::mf5_tests` のように引いているため）
- 実行されるテスト: クレートごとに `cargo test -p xlc-<crate> -- --list` の並べ替えた出力が、移す前後で完全に同じ（core 151、docs 57
  （単体 52 と統合 5）、format 66、lattice 334、numeric 40、xs 315（単体 305 と main.rs の 10））
- 実行ファイル: 分岐点と移したあとの `cargo build --release --workspace` の xs、lattice、core を `strip` して比べると、3 つともバイト単位で
  一致した（`.text` と `.rodata` は strip の前でも一致）。位置の情報の差は strip で消える。ADR-0159 の build の段階（strip なし）の実行ファイルは、
  main の 703e611c と一致しなかった。節ごとに比べると、`.text` などの割り当てられる節は一致し、`.symtab` と `.strtab` の
  シンボル名の `.llvm.<数字>` の添字だけが違った。この差を ADR-0160 が strip で取り除き、903ba01c の build の段階が作る実行ファイルの md5
  （xs `e77257b3…`、lattice `85a3c566…`、core `fe2d8761…`）が、移したあとでも main の `dvc.lock` と一致した（2026-10-14）
- テストだけの項目 4 個（format の `sample`、`sample_pins`、numeric の `modes_at`、`flux`）は、それぞれ fewgroup/tests.rs（`pub(crate) mod tests;`）と
  pins/node/tests.rs の `impl` に移した

## 帰結

- 実装（64 個のファイルの移動と検査）は、この ADR のあとの別のブランチで行う（#0089）。テストの移動と検査の追加を別のコミットにする
- 2026-10-12 に入口でテストを持つのは、xs の `main.rs`（`src/main/tests.rs` になる）と crates/docs の `lib.rs`
  （`src/tests.rs` になる）である
- ADR と issue には、テストの場所をファイルのパスで書いたものがある（例: ADR-0132 の crates/core/src/nem.rs）。
  採用済みの ADR は直せない（docs/adr/README.md の「運用」で、直せるのは誤字、リンク、`status` と `superseded-by` だけ）ので、決めた時点の
  場所の記録として残る。proposed の ADR と開いた issue のパスは、移すブランチで `tests.rs` に直す
- DVC の計算の段階は、deps に `crates/*/src` でなく、build の段階（直下の dvc.yaml）がビルドして `bin/` に置く実行ファイルを持つ
  （ADR-0159）。テストを移しても、strip した実行ファイルがバイト単位で変わらなければ（ADR-0160。「確認」のとおり一致した）、
  計算の段階は古くならない。移すブランチで `dvc repro -P` は流さない。直下の dvc.lock は、build の段階の deps（`crates`）の md5 だけが
  変わり、outs（`bin/` の md5）は main のまま変わらない
- 実行ファイルが一致しなければ、移動でなく、実行ファイルを変えた別の原因（並列の和の順序（#0060）、ビルドの設定など）を調べる。
  古い段階を `dvc commit` で deps だけ取り込むこともしない。`dvc.lock` が、記録した deps のコードで走らせていない出力を指すことに
  なる（ADR-0130 の案 L と同じ形）
