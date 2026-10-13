---
status: open
labels: [performance, lattice, core, xs]
opened: 2026-10-12
closed:
milestone:
blocked-by: []
resolution:
---

# 出力を変えない範囲で、測ってから速くする

## 問題

出力をバイト単位で変えない高速化（メモリの確保の使い回し、ホットな経路の `clone` を除くことなど）を、まだ探して
いない。ビルドの設定（LTO、`codegen-units`）はインライン化を変えるので、出力が変わらないかは走らせて確かめる。コードのレビュー（2026-10-12）は `.clone()` が多いと指摘した（1842033d の `crates` の `*.rs` で、
テストを含めて 483 個）が、どれがホットな経路かは測っていない。

分かっている時間の内訳: VERA 2A の 1 状態は 118 秒で、自己遮蔽 110 秒、MOC 7.3 秒（examples/vera/verification.md、
1f4fcaaf、2026-10-12）。自己遮蔽の中は、VERA 2B の 1 状態の perf で `cp::Cylinder::tracks` 45%、`cp::interpolate` 23%（同じ
ファイルの「対称性」、2026-10-10）。自己遮蔽の時間の内訳を測ることは #0039 と同じ作業なので、
ここで測った記録を #0039 からも引く。ここで扱うのは、そのうち出力を変えない高速化である。xs の時間は、examples/branches/verification.md（13 分、32 スレッド）と examples/m3-nuclides/verification.md（29 核種 × 4 温度で
57 秒）にある。core の時間は記録がない。ワークスペースの Cargo.toml に `[profile.release]` の設定はなく
（LTO と `codegen-units` は既定）、コードの生成の指定は、ADR-0109 の `target-cpu=znver5` と、ADR-0157 のループの先頭の整列（64 バイト）である。

結果を変える高速化（#0039 の、衝突確率を少数の全断面積から内挿する案など）は、ここでは扱わない。

関係する ADR: [ADR-0144](../adr/0144-refactor-and-tune-with-bit-identical-outputs.md)

## 完了の条件

（ソースを変える高速化は、同じファイルに触るテストの移動 #0089 のあとに行う。perf の測定とビルドの設定は待たないので、
`blocked-by` には入れない）

- xs、lattice（自己遮蔽を含む）、core の例の入力で、perf の時間の内訳を測って記録する
- 内訳の上位の箇所とビルドの設定について、出力を変えない高速化を ADR-0144 の判定で入れるか、効かないと判断する
- 入れたブランチ（ADR-0144 の「確認」に書く最初のブランチを除く）ごとに、ADR-0144 の判定（全桁の msgpack のビット
  一致と、時間の比較）の結果がここにある

## 入れた高速化

### J, K 積分の表の構築（xs）

`xs::resonance::line_shape` の J, K の表（β 641 点 × θ 361 点、θ の行ごとに適応 Gauss-Kronrod 積分）を、θ の行について
rayon で並列に作る（ブランチ perf/0088-line-shape-table）。行どうしで和を取らないので、各点の値はスレッド数によらない。
この表の構築が #0091 の遅いテストの時間を決めていた。

比べた候補:

- (a) θ の行を rayon で並列に作る。採った。積分の許容誤差と格子はそのままで、値を変えない。xs にはほかに rayon がないので
  依存を足したが、ワークスペースの依存（lattice と core が使う）で、Cargo.lock に新しい crate は入らない。
  `std::thread::scope` で行を分けても同じことはできるが、行の重さが θ で違うので、仕事の分け方を自分で書くことになる
- (b) 適応積分の許容誤差（1e-12）や表の間隔（`STEP`）を変える。値が変わるので ADR-0144 の範囲の外。採らない
- (c) 表を事前に計算してリポジトリに持つ。ln J, ln K の 2 × 641 × 361 点の f64 で約 3.7 MB になり、読み込みの形式、
  コードから作り直した表と一致するかの検査（それ自体が今の構築の時間を払う）、積分のコードを変えたときの作り直しの手順が要る。
  release の構築は下のとおり 0.27 秒で、本番の時間にはほとんど効かず、得るのはテストの時間だけなので採らない
- (d) 必要な点だけ遅延で作る。表の補間は 4 × 4 点を引くので、θ の行（β の全点を 1 回の積分で求める）単位の遅延になり、
  行ごとの `OnceLock` と排他が要る。`table_matches_direct_integration` は乱数の 400 点で θ の広い範囲を引き、ほかの
  テストと本番でどれだけの行を引くかは測っていないので効きが読めず、(a) で足りるので採らない
- 候補の外: テストの debug ビルドだけ xs を最適化する（`[profile.dev.package.xlc-xs]`）。ワークスペースの設定で、
  デバッグのしやすさも変わるので、ここでは試していない

ADR-0144 の判定（d49fe6a5、分岐点 6b39ad13、2026-10-13、r995 の Slurm の排他のジョブ 268）: 分岐点と変更後の release の
`xs` に examples/m0/xs.toml（H-1, O-16, U-235, U-238、293.6 K、σ0 = 1e10。U-235 と U-238 の非分離共鳴で表を引く。
出力のパスだけ変えた写し）を与え、`RAYON_NUM_THREADS=32` で 2 回ずつ走らせた。変更後は `RAYON_NUM_THREADS=1` でも 1 回走らせた。

| 走行                       | 時間（秒） | 出力の md5                         |
| -------------------------- | ---------- | ---------------------------------- |
| 分岐点 1                   | 55.32      | `1943d0c43e0a5023686b829cc4029967` |
| 分岐点 2                   | 55.00      | 同じ                               |
| 変更後 1                   | 54.32      | 同じ                               |
| 変更後 2                   | 54.92      | 同じ                               |
| 変更後（rayon 1 スレッド） | 54.05      | 同じ                               |

5 つの `wims69` の msgpack は `cmp` でバイト単位で一致した。xs の段階の時間の差（約 0.5 秒）は、表の構築の差（0.24 秒）と
走行のばらつき（同じ実行ファイルで 0.3〜0.6 秒）の程度で、段階の時間は U-238 の処理（約 55 秒）が決めている。
遅くはなっていない。

表の構築: release のテスト `table_matches_direct_integration`（表の構築と 400 点の直接の積分）は、同じジョブで 3 回ずつ
0.26〜0.27 秒から 0.03 秒になった。

debug のテスト（手元、`cargo test -p xlc-xs --lib -- --exact <テスト>`、load average 3〜5、相対の目安）:

| テスト                                                                                | 分岐点（秒） | 変更後（秒） |
| ------------------------------------------------------------------------------------- | ------------ | ------------ |
| `line_shape` の 8 本（`cargo test -p xlc-xs --lib line_shape`）                       | 7.5〜8.1     | 2.5〜3.1     |
| `line_shape::tests::table_matches_direct_integration`                                 | 7.42         | 0.92         |
| `unresolved::tests::capture_is_shielded_more_at_small_background_and_low_temperature` | 7.36         | 0.97         |
| `process::tests::unresolved_shielding_lowers_capture_at_small_background`             | 7.61         | 1.04         |
| `process::tests::processes_each_temperature_and_sigma0`                               | 15.93        | 4.06         |
| `process::tests::scattering_order_sets_the_number_of_moments`                         | 11.99        | 5.58         |

## 経緯

- 2026-10-12: コードのレビューを受けて起こした。ADR-0144 で進め方を決めた
- 2026-10-13: 関数の置き場所で時間が 30% 変わることが分かり（#0105）、ループの先頭を 64 バイトに揃えるビルドの設定を入れた（ADR-0157。測った M0 の xs、lattice、core と VERA 1B、2A の出力はバイト単位で変わらない）
- 2026-10-13: J, K 積分の表の構築を θ の行について並列にした（ブランチ perf/0088-line-shape-table。「入れた高速化」）
