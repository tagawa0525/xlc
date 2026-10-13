---
status: open
labels: [xs, test]
opened: 2026-10-12
closed:
milestone:
blocked-by: []
resolution:
---

# 入力を縮めても速くならない遅いテストの原因と、速くする方法

## 問題

遅いテストの入力を縮める作業（ブランチ test/shrink-slow-tests）で、次のテストは縮めても速くならない、または縮めると
確かめる性質が変わるので縮めなかった。時間は負荷の高い手元（load average 70 以上）で `cargo test -p <crate> <テスト名>`
を 1 本ずつ測った user の CPU 時間で、相対の目安にだけなる。

関係する ADR: なし

### J, K 積分の表の構築が時間を決めるテスト

`xs::resonance::line_shape` の J, K の表（`table()`、`OnceLock`）は、最初に使うときに一度だけ作る。表は
β 641 点 × θ 361 点で、各点で適応 Gauss-Kronrod 積分する。debug ビルドでは作るのに 15 秒ほどかかり、次のテストの時間は
ほぼこの構築で決まっていた（下の時間は 2026-10-12 のもの）。θ の行を並列に作るようにして、1 本ずつ走らせた時間は 0.6〜5.1 秒に
なった（#0088 の「入れた高速化」）。

- `process::tests::processes_each_temperature_and_sigma0`（21 秒）
- `process::tests::scattering_order_sets_the_number_of_moments`（23 秒）
- `process::tests::unresolved_shielding_lowers_capture_at_small_background`（15 秒）
- `unresolved::tests::capture_is_shielded_more_at_small_background_and_low_temperature`
- `line_shape::tests::table_matches_direct_integration`（表の細かさと補間の精度そのものが検証の対象）

試した縮め方: `process` の 2 つの群構造を WIMS 69 群から 6 群（2.5e4, 1.5e4, 1e3, 4, 0.1, 1e-5 eV）に減らしても、時間は
23 秒から 21 秒、21 秒から 19 秒でほとんど変わらなかった（perf で `jk_integrals` が大半）。群数は時間を決めておらず、
入力を縮める余地がない。

`cargo test` は 1 つのプロセスで表を共有するので、表の構築は全体で 1 回しか払わない。nextest はテストごとにプロセスを
分けるので、これらのテストが並列でそれぞれ構築し、各テストの時間が長く見える。

速くする案（どれも本番のコードの変更なので、ADR で決める）: 表の構築を θ について並列にする（入れた。rayon は捨てて
`std::thread::scope` にした。理由は #0088）。表を事前に計算して持つ（#0088 で採らなかった）。`STEP` を粗くする案は、補間誤差（`table_matches_direct_integration` の 1e-6）との関係を測ってから。

### 表の細かさが検証の対象のテスト

`xs::multigroup::free_gas::tests::free_gas_moments_match_the_s_alpha_beta_table`（14 秒）は、細かい自由ガスの表
（`free_gas_table`、α 801 点 × β 752 点）の補間が解析解と 3e-3 で一致することを確かめる。表を粗くすると補間誤差が増え、
確かめる性質が変わる。時間は `ThermalSab::kernel_moments` が決めている。

### S(α, β) の群への積分が時間を決めるテスト

`xs::multigroup::scatter::tests::sab_ignores_mf3_elastic_below_cutoff`（14 秒）は、`sab_matrix` を 2 回作る。確かめる行は
cutoff（1 eV）以下の全群で、群を 1 eV 以下に絞っても時間は変わらなかった（16 秒）。減らせる入力がない。

## 完了の条件

J, K の表の構築を速くする（または nextest の各プロセスで構築しない）方法を ADR で決めて実装し、上のテストの時間が
縮むこと。

## 経緯

- 2026-10-12: ブランチ test/shrink-slow-tests で、入力を縮められるテストは縮め、縮められないものをここに書いた
- 2026-10-13: J, K の表の構築を θ の行について並列にし、上の 5 本の debug の時間が縮んだ（数値は #0088 の「入れた高速化」）
