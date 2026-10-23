---
status: open
labels: [xs]
opened: 2026-10-12
closed:
milestone:
blocked-by: []
resolution:
---

# m0 と reflector の xs を核種ごとのステージにする

## 問題

vera の xs は範囲、核種ごと、組み立てのステージに分けた（[ADR-0138](../adr/0138-split-xs-into-nuclide-stages.md)）。
examples/m0/dvc.yaml（4 核種）と examples/reflector/dvc.yaml（100 近い核種）の `xs` は、1 つのステージのままで、
1 つの核種の項目を変えても全核種を処理し直す。どちらも ENDF のファイルを 1 つずつ deps に手で並べている。

関係する ADR: [ADR-0138](../adr/0138-split-xs-into-nuclide-stages.md)

## 完了の条件

examples/reflector/dvc.yaml の `xs` を、vera と同じ形（`xs.toml` を `vars` にして核種ごとに `foreach`。`scattering-order` を
`xs.toml` に書く）の 3 ステージに分け、核種 1 つの変更でその核種と組み立てだけが走ることを `dvc status` で確かめる。
m0 は核種が 4 つなので、変更が続くようになってから（分けても得るものが小さい）。

## 経緯

- 2026-10-12: #0075 で vera だけを分け、残りを ADR-0138 の判断で、この issue に分けた
