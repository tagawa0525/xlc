---
status: open
labels: [performance, lattice, core, xs]
opened: 2026-10-12
closed:
milestone:
blocked-by: ["#0089"]
resolution:
---

# 出力を変えない範囲で、測ってから速くする

## 問題

出力をバイト単位で変えない高速化（メモリの確保の使い回し、ホットな経路の `clone` を除くことなど）を、まだ探して
いない。ビルドの設定（LTO、`codegen-units`）はインライン化を変えるので、出力が変わらないかは実行して確かめる。コードのレビュー（2026-10-12）は `.clone()` が多いと指摘した（14793ea7 の `crates` の `*.rs` で、
テストを含めて 483 個）が、どれがホットな経路かは測っていない。

分かっている時間の内訳: VERA 2A の 1 状態は 118 秒で、自己遮蔽 110 秒、MOC 7.3 秒（examples/vera/verification.md、
e4dbc0de、2026-10-12）。自己遮蔽の中は、VERA 2B の 1 状態の perf で `cp::Cylinder::tracks` 45%、`cp::interpolate` 23%（同じ
ファイルの「対称性」、2026-10-10）。自己遮蔽の時間の内訳を測って速くする方法を決めることは #0039 と同じ作業なので、
自己遮蔽の内訳は #0039 で測って記録し、ここからはそれを引く。ここで扱うのは、そのうち出力を変えない高速化である。xs の時間は、examples/branches/verification.md（13 分、32 スレッド）と examples/m3-nuclides/verification.md（29 核種 × 4 温度で
57 秒）にある。core の時間は記録がない。ワークスペースの Cargo.toml に `[profile.release]` の設定はなく
（LTO と `codegen-units` は既定）、CPU の指定は ADR-0109 の `target-cpu=znver5` だけである。

結果を変える高速化（#0039 の、衝突確率を少数の全断面積から内挿する案など）は、ここでは扱わない。

関係する ADR: [ADR-0144](../adr/0144-refactor-and-tune-with-bit-identical-outputs.md)

## 完了の条件

- xs、lattice（自己遮蔽を除く。自己遮蔽は #0039 の記録を引く）、core の例の入力で、perf の時間の内訳を測って記録する
- 内訳の上位の箇所とビルドの設定について、出力を変えない高速化を ADR-0144 の判定で入れるか、効かないと判断する
- 入れたブランチ（ADR-0144 の「確認」に書く最初のブランチを除く）ごとに、ADR-0144 の判定（全桁の msgpack のビット
  一致と、時間の比較）の結果がここにある

## 経緯

- 2026-10-12: コードのレビューを受けて起こした。ADR-0144 で進め方を決めた
