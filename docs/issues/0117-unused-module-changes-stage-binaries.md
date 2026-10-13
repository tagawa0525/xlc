---
status: open
labels: [build, dvc]
opened: 2026-10-14
closed:
milestone:
blocked-by: []
resolution:
---

# 依存するクレートに使わないモジュールを足すと段階の実行ファイルが変わる原因を調べる

## 問題

`xlc-format` に、xs、lattice、core のどれも使わないモジュール（serde の構造体と読み書きの関数だけのチェーンの型）を足して
`scripts/build-binaries.sh` でビルドすると、`bin/xs`、`bin/lattice`、`bin/core` の md5 がすべて変わった（main の 1f50ae7a の
作業ツリーでコミットせずに試した。2026-10-14。ADR-0166 の背景）。モジュールを外すと元に戻る。段階の実行ファイルが変わると
計算の段階が古くなり、xs の段階（60 を超える）を再実行することになる。このため ADR-0166 はチェーンの型を `xlc-format` でなく
`xlc-depletion` に置いた。共通のクレートに型や関数を足すたびに同じことが起きうる。

原因の候補（確かめていない）: codegen unit の分け方（release は LTO なし、codegen unit は既定の 16）で、モジュールを足すと
既存の関数の入る codegen unit が変わる。

関係する ADR: なし

## 完了の条件

- 実行ファイルが変わる原因を確かめてある（変わった節と記号を比べる）
- 使わないコードを足しても実行ファイルが変わらない設定（codegen unit の数、LTO など）があるか、その費用（ビルドの時間、
  計算の速さ）と合わせて決め、決めたことを ADR にしてある

## 経緯

- 2026-10-14: 燃焼チェーンの型の置き場を決めるときに見つけた（ADR-0166）
