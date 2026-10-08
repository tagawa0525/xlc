---
status: superseded
date: 2026-10-06
supersedes: []
superseded-by: [ADR-0010, ADR-0011]
issues: []
---

# ADR-0009: 段階間の中間ファイルは自前のテキスト形式とする

## 背景

段階間はファイルで受け渡す（[ADR-0001](0001-split-into-three-stage-crates.md)）。M0 の縦切りで、その形式を決める必要があった。

## 検討した案

- 自前のテキスト形式（採用）

## 決定と理由

中間ファイルは当面、自前のテキスト形式とする。人が読めて diff を取れるため、検証の段階で便利である。サイズや速度が
問題になった時点で見直す。M0 では次を書いた。

- 多群ライブラリ（`xs` → `lattice`）: 群境界、核種ごとの σt, νσf, σf, χ, 散乱行列
- 少数群定数（`lattice` → `core`）: 2群の Σt, Σa, νΣf, Σf, χ, 散乱行列, D

## 確認

M0 の例（examples/m0/README.md）で3段階をつないで確かめた。

## 結果

M2 で、データファイルは MessagePack（[ADR-0011](0011-write-data-files-in-messagepack.md)）、入力ファイルは TOML（[ADR-0010](0010-write-inputs-in-toml.md)）に改め、自前のテキスト形式は
M2 の中で移行した。
