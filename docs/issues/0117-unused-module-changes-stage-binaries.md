---
status: open
labels: [build, dvc]
opened: 2026-10-14
closed:
milestone:
blocked-by: []
resolution:
---

# 依存するクレートに使わないモジュールを足すとステージの実行ファイルが変わる原因を調べる

## 問題

`xlc-format` に、xs、lattice、core のどれも使わないモジュール（serde の構造体と読み書きの関数だけのチェーンの型）を足して
`scripts/build-binaries.sh` でビルドすると、`bin/xs`、`bin/lattice`、`bin/core` の md5 がすべて変わった（main の 70e6ffa8 の
作業ツリーでコミットせずに試した。2026-10-14。ADR-0166 の背景）。モジュールを外すと元に戻る。ステージの実行ファイルが変わると
計算のステージが古くなり、xs のステージ（60 を超える）を再実行することになる。このため ADR-0166 はチェーンの型を `xlc-format` でなく
`xlc-depletion` に置いた。共通のクレートに型や関数を足すたびに同じことが起きうる。

試したモジュール（`crates/format/src/chain.rs` に書き、`lib.rs` に `pub mod chain;` を足して、default の devShell で
`scripts/build-binaries.sh`）:

```rust
use serde::{Deserialize, Serialize};
use crate::msgpack::{self, MsgpackError};
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Product { pub nuclide: u32, pub count: u32 }
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Decay { pub mode: String, pub ratio: f64, pub products: Vec<Product> }
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Nuclide { pub name: String, pub half_life: Option<f64>, pub decays: Vec<Decay> }
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Chain { pub format: String, pub nuclides: Vec<Nuclide> }
impl Chain {
    pub fn to_msgpack(&self) -> Result<Vec<u8>, MsgpackError> { msgpack::to_named(self) }
    pub fn from_msgpack(bytes: &[u8]) -> Result<Self, MsgpackError> {
        let v = msgpack::read_value(bytes)?;
        let c: Chain = msgpack::from_value(v)?;
        for n in &c.nuclides { for d in &n.decays { for p in &d.products {
            if p.nuclide as usize >= c.nuclides.len() { return Err(MsgpackError::NotMap); } } } }
        Ok(c)
    }
}
```

md5 は xs: e77257b3 → f84383de、lattice: 03b6048d → 260f3d16、core: fe2d8761 → 0c2b88d6。

原因の候補（確かめていない）: codegen unit の分け方（release は LTO なし、codegen unit は既定の 16）で、モジュールを足すと
既存の関数の入る codegen unit が変わる。

関係する ADR: ADR-0166

2026-11-18 の main（8369db18。直下の dvc.lock の `build` を付け直した後）では、examples/depletion/p2/dvc.yaml の xs のすべての核種のステージと
`xs`、examples/depletion/dvc.yaml の `chain` が、deps の `bin/xs`（dvc.lock は f98234e2、`build` の出力は 28e28365）と `bin/chain`（4a9d4768 と
d14e59ff）だけで changed で、出力（depletion.msgpack 975be970、chain.msgpack 1f1de180）は dvc.lock のとおり。xs と chain のクレートの変更で
変わったのか、ここの原因で変わったのかは確かめていない。その下流のステージを `dvc repro` すると、ライブラリの作り直しが先に走る（#0176 の
Gd の確認は、格子計算をステージにせず流した）

## 完了の条件

- 実行ファイルが変わる原因を確かめてある（変わった節と記号を比べる）
- 使わないコードを足しても実行ファイルが変わらない設定（codegen unit の数、LTO など）があるか、その費用（ビルドの時間、
  計算の速さ）と合わせて決め、決めたことを ADR にしてある

## 経緯

- 2026-10-14: 燃焼チェーンの型の置き場を決めるときに見つけた（ADR-0166）
