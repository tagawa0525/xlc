---
status: open
labels: [vera, openmc]
opened: 2026-10-12
closed:
milestone: M3
blocked-by: []
resolution:
---

# バランス（openmc_balance.py）の世代ごとの値を、参照解の実行から取る

## 問題

k の差の内訳（examples/vera/verification.md の「k の差の内訳」）は、`openmc_balance.py` が参照解とは別に OpenMC を
実行して作る `balance.msgpack` を読む。標準のタリーの組（ADR-0152）は、物質 × 71 のビンの中性子束、吸収、生成、
2 個以上の中性子を出す反応（`standard-1/material groups`）と、物質 × 入射 × 出射の nu-scatter（`material transfer` の P0）を
持つので、`balance.msgpack` の平均（`materials` の flux、absorption、production、transfer と `outside`）は、参照解の
statepoint 1 つから作れる。

作れないのは統計誤差である。`compare_balance.py` の σ は世代ごとの値（`batches`）から求め、生成、吸収、中性子束が同じ
粒子の追跡から出る相関と、生成の和での規格化の相関を含む。OpenMC の statepoint は、ビンごとに世代の値の和と 2 乗の和
（`results` の最後の次元の 2 つ）しか持たず、ビンの間の共分散を持たない。物質と範囲の和（ΣΔA）の σ は、ビンの間の
相関を落とすと、過大にも過小にもなりうる（どちらかは確かめていない）。そこで `openmc_balance.py` は世代ごとに
statepoint を書き（`settings.statepoint`）、続く 2 つの差から世代の値を求める。

`balance-1B` が 20 GB（ADR-0150 の背景）なのは、204 の statepoint がそれぞれ源の粒子のバンク（100 万粒子 × 104 バイト
= 104 MB）を持つためで、バランスのタリーは 1 つの statepoint で 0.23 MB（`references/openmc/balance-1B/statepoint.100.h5`、
2026-10-12 に h5py で見た）。参照解の実行で世代ごとに statepoint を書くと、標準の組（ADR-0152 の帰結。集合体で十数 MB）が
世代ごとに入るので、2000 世代では数十 GB になり、同じ手は使えない。

候補（どれも試していない）:

- 参照解の実行を `openmc.lib` で世代ごとに進め（`openmc.lib.iter_batches()`）、世代の終わりにバランスのタリーの和だけを
  読んで小さなファイルに足していく。statepoint は今と同じに書く。#0073 の再開（`openmc -r`）と合わせて、再開の前の世代の
  値をどこに残すかを決める必要がある。#0073 の `run_model`（`model.run`）を置き換えるので、#0073 がマージされてから行う
- `openmc_balance.py` を残し、世代ごとの statepoint に源の粒子のバンクを書かない（`settings.sourcepoint = {"write": False}`）。
  別の実行は残るが、20 GB は約 50 MB（0.23 MB × 204）になる
- 世代をまとめた平均（batch means）で σ を見積もる: 参照解の実行で数十世代ごとに statepoint を残す。世代の間の相関も
  含められるが、残す statepoint の数 × 標準の組の大きさになる

関係する ADR: ADR-0152（フロントマターの `issues` にこの issue がある）

## 完了の条件

k の差の内訳の σ（相関を含む）を、参照解と同じ実行から求められるようにするか、`openmc_balance.py` の別の実行を残すと
決めて、その出力を小さくする。

## 経緯

- 2026-10-12: 標準のタリーの組（ADR-0152）を決めたときに、`openmc_balance.py` の別の実行をなくせるかを調べた。平均は
  作れるが、世代ごとの値が statepoint に残らないので σ を作れない
