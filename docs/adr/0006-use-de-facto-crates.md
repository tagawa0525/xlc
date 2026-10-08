---
status: accepted
date: 2026-10-06
supersedes: [ADR-0005]
superseded-by:
issues: []
---

# ADR-0006: 事実上の標準があるものはそのクレートを使い、定番のない数値計算は自前で実装する

## 背景

最初は標準ライブラリを優先し、数値計算も自前とした（[ADR-0005](0005-prefer-std-and-own-numerics.md)）。並列処理、複素数、ファイルの読み書き、エラー、
ログを自前で書くと、その保守を抱えることになる。

- 保守のコスト

## 検討した案

- 標準ライブラリを優先し、自前で書く（[ADR-0005](0005-prefer-std-and-own-numerics.md)）: 保守が重い
- 事実上の標準のクレートを使い、定番のないものだけを自前で書く（採用）

## 決定と理由

依存クレートは、事実上の標準があるものはそれを使う。並列処理の rayon、複素数の num-complex、ファイルの読み書きの
serde、エラーの thiserror と anyhow、ログの tracing など。自前で書いたコードの保守の方が、定番のクレートへの依存
より重いためである。

定番のない数値計算（ENDF の処理、共鳴の計算、原子炉の分野の特殊関数、CMFD の演算子と反復の組み立て、CRAM など）は、
問題に特化しているので自前で実装する。中の線形問題は faer で解く（[ADR-0007](0007-use-faer-for-linear-algebra.md)）。

複素数は num-complex とする（faer も num-complex の複素数を扱える）。xs の共鳴の計算にあった自前の複素数の型は
置き換えた。

## 確認

方針で、確かめるテストはない。

## 結果

- 線形代数は [ADR-0007](0007-use-faer-for-linear-algebra.md)、小さな密行列は [ADR-0008](0008-keep-own-small-dense-solvers.md)
- エラーは [ADR-0016](0016-type-errors-with-thiserror-and-anyhow.md)、ログは [ADR-0018](0018-log-with-tracing-to-stderr.md)
