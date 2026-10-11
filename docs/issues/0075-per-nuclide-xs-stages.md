---
status: open
labels: [xs]
opened: 2026-10-12
closed:
milestone:
blocked-by: []
resolution:
---

# xs を核種ごとの段階にし、ライブラリの組み立てを段階にする

## 問題

examples の dvc.yaml（examples/m0、vera、reflector）で、xs は `xs.toml` の全核種を 1 つの段階で処理する。1 つの核種の
入力や扱いを変えただけでも全核種を処理し直し、下流の lattice もすべて走り直す。OpenMC のライブラリの段階は、すでに
核種ごとの段階（`openmc-neutron`、`openmc-thermal`）と組み立て（`openmc-index`）に分けてある（ADR-0130）。

関係する ADR: [ADR-0123](../adr/0123-drive-stages-with-dvc.md)、[ADR-0130](../adr/0130-drive-reference-results-with-dvc.md)

## 完了の条件

xs を核種ごとの出力を作る段階と、それを 1 つのライブラリにまとめる段階に分け、核種 1 つの変更でその核種と組み立てだけが
走ることを `dvc status` で確かめる。分け方（xs の CLI、核種ごとの出力の形式）は ADR にする。

## 経緯

- 2026-10-11: DVC の段階（ADR-0123）の後の作業として決めた（メモリに置いていたものを 2026-10-12 に issue にした）
