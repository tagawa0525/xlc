---
status: open
labels: [xs]
opened: 2026-10-08
closed:
milestone: M3
blocked-by: []
resolution:
---

# 散乱行列の Legendre の次数が xs の入力にない

## 問題

[ADR-0076](../adr/0076-store-legendre-scattering-moments.md) は、ライブラリの散乱行列の Legendre の次数 L を核データの処理の入力（既定 3）と決めている。実装では、xs の
コードに `SCATTERING_ORDER = 3` と固定してあり（crates/xs/src/main.rs）、xs.toml に項目がない。格子計算の側は
`scattering_order`（省略するとライブラリにある全ての次数）で選べる。

関係する ADR: [ADR-0076](../adr/0076-store-legendre-scattering-moments.md)

## 完了の条件

xs の入力に次数の項目（省略すると 3）を加え、次数を変えたライブラリを作れることをテストで確かめる。

## 経緯

- 2026-10-08: ADR のレビューで、決定と実装の食い違いとして見つかった
