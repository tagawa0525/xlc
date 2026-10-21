---
status: open
labels: [xs, lattice, verification]
opened: 2026-10-11
closed:
milestone:
blocked-by: []
resolution:
---

# 並列の和の順序を固定して、出力をスレッド数によらずビット一致させる

## 問題

DVC（ADR-0123）は、段階の出力を md5、つまりバイトの一致で記録し、数値の許容幅をもたない。同じ入力とコードでも、
並列の和の順序が実行方法（スレッド数。既定は物理コアの数で、マシンと Slurm の割り当てで変わる。ADR-0067）で変わると、
出力の md5 が変わる。再実行の引き金は deps なので、走らせなければ出力はキャッシュから戻って揺らがないが、次の 3 つで
効く。

- 上流の段階を再実行したとき、結果の意味が同じでも出力の md5 が変わり、下流の段階まで走る
- `dvc.lock` に意味のない差分が出て、別々に実行したブランチの間でぶつかる
- DVC の外で実行した出力が、`dvc status` で「changed outs」になる

2026-10-11 に `crates/*/src` を `par_iter`、`par_chunks`、`thread::scope`、`HashMap`、`reduce` で grep した並列の箇所と、
その和の順序（確かめたのはコードとテストで、出力のバイト一致は確かめていない）:

| 箇所                                                                                            | 中身                                         | 和の順序                                                                                                                                                    |
| ----------------------------------------------------------------------------------------------- | -------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| crates/lattice/src/moc.rs（スイープと `Sums` の集計）                                           | 代表のトラックの組ごとの並列のスイープ       | 組の分け方と、組の順・像の順の和は固定（ADR-0064）。`*result_does_not_depend_on_the_thread_count` の 4 本が、1 スレッドと複数のスレッドで解の一致を確かめる |
| crates/lattice/src/ufg.rs（エネルギー点ごとのトラック）                                         | `par_iter().map().collect()`                 | 点ごとに独立に求め、`collect` が順を保つ。和は逐次のループ                                                                                                  |
| crates/lattice/src/shielding.rs（ピンごとの減速計算）                                           | `par_iter().map().collect()`                 | ピンごとに独立に解き、結果はピンの順に逐次に扱う                                                                                                            |
| crates/xs/src/main.rs（核種ごとの読み込みと処理）                                               | `std::thread::scope` で核種ごとに 1 スレッド | 核種ごとに独立で、結果は入力の順に集める                                                                                                                    |
| crates/lattice/src/shielding.rs、crates/core/src/input.rs、crates/core/src/pins.rs の `HashMap` | 番号の割り当てとキャッシュ                   | 引くだけで、繰り返しの順に依存しない                                                                                                                        |

順序に依存する並列の和は、grep では見つからなかった。足りないのは、段階の出力（msgpack）がスレッド数によらず
バイト一致することの確認である。MOC のテストは解（`solve` の戻り値）の一致を見ていて、自己遮蔽、超詳細群、
2群への縮約、書き出しまでを通した出力は見ていない。xs の段階の出力の一致も確かめていない。

関係する ADR: [ADR-0123](../adr/0123-drive-stages-with-dvc.md)

## 完了の条件

- lattice と xs の段階の出力の msgpack が、スレッド数（1 と複数）を変えてもバイト一致することを確かめるテストがある
- 一致しない箇所があれば、和の順序を固定して直す（修正方法は ADR にする）

## 経緯

- 2026-10-11: ADR-0123 で段階を DVC で駆動し、出力をバイトの一致で記録することにしたので起こした。M0 の lattice は、
  同じスレッド数で 2 回実行して `pincell-fg.msgpack` の md5 が一致した（ADR-0123 の確認）
- 2026-10-14: feat/0060-dvc-examples に m1、mox、c5g7、colorset、split、m3-nuclides の dvc.yaml と OpenMC の
  model/run/extract の書き換えをコミットした
- 2026-10-21: mox と c5g7 を Slurm で 1 度実行して dvc.lock を作った（OpenMC の出力が手元になかったため。mox の OpenMC の k は
  前の測定と 6 桁で同じ）。m1、colorset、split、m3-nuclides も同日に Slurm の排他のジョブ 1 本（1 時間 16 分）で実行し、dvc.lock を作った。
  追跡している出力（比較の JSON など）は前の値と変わらなかった
