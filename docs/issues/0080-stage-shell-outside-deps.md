---
status: open
labels: [verification]
opened: 2026-10-12
closed:
milestone:
blocked-by: []
resolution:
---

# 段階を実行する devShell の設定が DVC の deps に表れない

## 問題

DVC のステージは default の devShell で走る（[ADR-0130](../adr/0130-drive-reference-results-with-dvc.md)）。段階の deps は
`rust-toolchain.toml` と `flake.lock` を含むが、`flake.nix` を含まない（[ADR-0123](../adr/0123-drive-stages-with-dvc.md)）。
そのため `flake.nix` で default の devShell の設定（shellHook の環境変数、C コンパイラ）を変えたり、ツールチェーンを
override したりすると、段階のビルドが変わっても段階は changed にならない。止めているのは flake.nix のコメントだけ
（ADR-0135）。

関係する ADR: [ADR-0135](../adr/0135-keep-dev-tools-out-of-rust-toolchain.md)

## 完了の条件

段階のビルドの環境を deps で表す方法を決める。候補は、段階の shell の store パスを記録する段階（`references-shell` と
同じ形）を Rust の段階の deps にする、ビルド用の shell を分ける、flake の check でツールチェーンが
`rust-toolchain.toml` だけから作られることを確かめる、など。

## 経緯

- 2026-10-12: ADR-0135 のレビューで、ADR-0123 の前からある限界として挙がった
