---
status: open
labels: [verification, lattice, xs]
opened: 2026-11-04
closed:
milestone: M3
blocked-by: []
resolution:
---

# VERA の格子と P3 の k が、ADR-0205 と ADR-0209 の後の再実行で KENO-VI に対して +34〜+154 pcm 動いた

## 問題

examples/vera と examples/vera/endfb70 のライブラリと格子を、ブランチ chore/rerun-after-adr-0226（流したコミットは 4cfe4f58、分岐点は
main の b4e9b2cd、2026-11-04、Slurm 678〜681）で今の bin で作り直すと、xlc の k が一様に上がり、KENO-VI との差が大きくなった。
前の実行は 05714cea（2026-10-23。examples/vera/verification-endfb70.md の前の版の「計算」）で、結果は 4ad8d9fbc でコミットした。
その dvc.lock の `bin/lattice` は 457aae7a で、05714cea の dvc.lock の値と違う（どのコミットで作った bin かは dvc.lock からは決まらない）。

xlc − KENO-VI（pcm。比較の JSON の `k.diff_pcm`、P3 は `k.diff_keno_pcm`。前 → 今）:

- ENDF/B-VIII.0（本番。`examples/vera/compare-keno-*.json`）: 1B −3.8 → +48.4、1C −31.2 → +66.5、1D −32.9 → +120.6、2A +36.8 → +90.8、
  2B +26.7 → +78.7、2C +15.2 → +108.8、2D −2.9 → +141.4、2E +97.9 → +142.4、2F +145.6 → +187.1、2G +159.1 → +195.0、2H +93.7 → +128.0、
  2I +45.9 → +94.1、2J +153.8 → +195.2、2K +98.7 → +141.6、2Q +146.3 → +202.9（xlc の k の変化は +34〜+154 pcm）
- ENDF/B-VII.0（`examples/vera/endfb70/compare-*.json`）: 1B +167.7 → +202.3、1C +134.1 → +213.9、1D +127.1 → +262.3、2A +145.7 → +185.2、
  2B +154.4 → +191.4、2C +137.2 → +215.3、2D +114.2 → +242.6、2E +121.8 → +154.4、2F +120.4 → +151.4、2G +80.3 → +109.5、2H +104.9 → +131.4、
  2I +142.2 → +176.7、2J +128.3 → +159.2、2K +126.2 → +158.4、2Q +245.1 → +285.6（xlc の k の変化は +27〜+135 pcm）
- P3（`examples/vera/p3/compare-3?.json`、`examples/vera/endfb70/p3/compare-3?.json`）: VIII.0 の 3A +62.0 → +110.1、3B +164.5 → +213.1、
  VII.0 の 3A +157.4 → +191.8、3B +136.7 → +172.8（xlc の k の変化は +34〜+49 pcm）

前と今の実行の deps（examples/vera/dvc.lock の md5）: `bin/lattice`（lattice-1B）457aae7a → 48387d65、`bin/xs`（xs）d3dc3807 → 953c587b、
本番のライブラリ `vera.msgpack` 05132f1d → 0dfa3aca、VII.0 のライブラリ（examples/vera/endfb70/dvc.lock）1d13e690 → 022595ef。前の
ライブラリは今の `bin/lattice` で読めない（`missing field sab_transfer`。ADR-0205 の後の形式）ので、ライブラリと格子の効きを今の bin で
分けることはできない。

関係する ADR: なし

## 原因（1C と 2A）

VIII.0 の 1C と 2A の xlc の k の変化（+97.7、+54.0 pcm）は、ADR-0209 が決定のときに既定の切り替えの前後の k の動きとして記録した値と、
同じ問題、同じ期間の変化で、0.1 pcm 以内で一致する。ADR-0209 の「確認」は、1C の燃焼で「xlc の k は 0 の点で +97.7（1 リング）、+97.2 pcm
（4 リング）…動き」、2A の集合体の燃焼で「xlc の k は 0 の点で +54.0 pcm…動き」と書く。前の実行は ADR-0209 の前の既定（静止した標的、
無限希釈の形）、今の実行は ADR-0209 の既定（`ufg-scattering = "free-gas"`、`ufg-free-gas` の速い方法と `sab = true`、
`ufg-transfer = "ring-flux"`）で解いている。

xlc の 0 の点の k（vera は `compare-keno-<問題>.json` の `k.xlc`、depletion は 1C が `compare-1C-rings1.json` の `k[0].xlc`、2A が
`p2/compare-2A.json` の `k[0].xlc`。表の括弧は値を入れた JSON のコミット。数値を出した xlc のコミットと日付は、vera の前が 05714cea
（2026-10-23、記録の「計算」）、今が 4cfe4f58（2026-11-04）、depletion の 1C の前が adacf8b4（13826a2b9 の本文の da2224ff と同じパッチ。rebase 前の hash で main から辿れない。分岐点 6822de4c、Slurm 457、
2026-10-23）、2A の前が 2c2b8c60
（Slurm 520、50f94ebfa の本文、2026-10-24）、ADR-0209 のときの 1C と 2A が xlc 73d03d37c（2026-10-28、ADR-0209 の「確認」。値は ae6fdf7fc の
JSON）、今が ADR-0226 の bin で、1C は Slurm 660（流したコミットと分岐点は main の 54ab434ff。06d6bd719 の本文、2026-11-03）、2A は Slurm 664
（流したコミット 24a26e32、分岐点 main の a0fd48c0。8488f5dfa の本文、2026-11-03。sacct の Start））:

| 問題 | 実行                                 | 前（ADR-0209 の前の既定） | 今                     | 変化（pcm） |
| ---- | ------------------------------------ | ------------------------- | ---------------------- | ----------- |
| 1C   | examples/vera の keno-pincell        | 1.171410（4ad8d9fbc）     | 1.172387（339858c46）  | +97.7       |
| 1C   | examples/depletion の 1 リング       | 1.1714054（13826a2b9）    | 1.1723822（06d6bd719） | +97.68      |
| 2A   | examples/vera の keno-assembly       | 1.182543（4ad8d9fbc）     | 1.183083（339858c46）  | +54.0       |
| 2A   | examples/depletion/p2 の集合体の燃焼 | 1.1825485（50f94ebfa）    | 1.1830889（8488f5dfa） | +54.03      |

depletion の今の値は、vera の今の実行と同じ `bin/lattice` 48387d65 で解き直したもの（main の dvc.lock）で、ADR-0209 のときの値（ae6fdf7fc の JSON）と
1e−12 以内で同じ。ADR-0209 の後から ADR-0226 の bin までの変更は、この 2 つの k を動かさない。vera と depletion の差は、1C で
前 +0.46、今 +0.48 pcm、2A で前 −0.55、今 −0.59 pcm で、ほぼ一定（ライブラリが違い、depletion の 2A は `ufg-tolerance = 1e-4` の間引いた格子。前は `bin/lattice` の md5 も違い、vera は
457aae7a、depletion の 1C の 13826a2b9 の dvc.lock は 7cf252f4）。
vera の出力は小数点以下 6 桁なので、差の差（0.02、0.04 pcm）は丸めの中。

- 同じ問題: examples/vera/p1/pincell-1C.toml と examples/depletion/lattice-1C-rings1.toml、examples/vera/p2/assembly-2A.toml と
  examples/depletion/p2/lattice-2A.toml は、形状、組成、温度、離散化が同じで、どれも UFG の散乱核と移行の項目を書かない（既定で解く）。
  違いはライブラリ（vera.msgpack と depletion.msgpack。どちらも ENDF/B-VIII.0、WIMS 69 群、同じ σ0）、depletion の燃焼の節、depletion の 2A の
  `ufg-tolerance`、vera の 2A のピン出力の出力（`pin-powers`）、2A の減速材の核種の順（examples/vera/p2/assembly-2A.toml は H-1 が先、
  examples/depletion/p2/lattice-2A.toml は O-16 が先）だけ
- 前の実行は、UFG の散乱核と移行の項目より前の bin: vera の前の結果のコミット 4ad8d9fbc の作成の日時は 2026-10-23 22:56、depletion の 2A の
  前の値のコミット 50f94ebfa は 2026-10-24 02:53、1C の 13826a2b9 は 2026-10-23 02:51。`ufg-scattering` を lattice に入れた最初のコミット
  （ADR-0203 の RED の 316052c1a）は 2026-10-24 13:35、`ring-flux` は 2026-10-28 00:55（0a377fa14）、ADR-0209 のマージは eba3a023b
  （2026-10-28 06:30）。どの日時も作成の日時（author date）
- 確かめていないこと: 上の一致は、vera と depletion が同じ期間に同じだけ動いたことを示し、ADR-0209 の値もその期間の差である。期間の中の
  既定の切り替えでない変更（2026-10-23 から 2026-10-28 のほかのマージ、ADR-0205 の形式を含むライブラリの作り直し）の効きは、ADR-0209 も分けていない。既定の
  結果を変えないことを、ADR-0203、ADR-0204、ADR-0207 はそれぞれの変更について、ADR-0209 は前の既定を書いた入力について確かめた

ADR-0209 はこの動きを OpenMC に対して測って受け入れた（k の差 xlc − OpenMC は、1C の 1 リングで −61〜+73 → +29〜+136 pcm、2A で
+5〜+203 → +56〜+164 pcm）。1C と 2A の KENO-VI に対する差の変化は、この既知の動きの再掲で、新しい未解明の差ではない。

ADR-0205 は別の候補から外す。ADR-0205 は `ufg-free-gas` に `sab`（省略すると false）を足しただけで、既定では使わない。ADR-0205 の
S(α,β) の移行の効き（ADR-0205 の決定の「k を自由ガスから 1C の 0 MWd/kg で +12.5、5 MWd/kg で +10.0、2A-40 で +7.5 pcm 動かし」）は、
ADR-0209 が `sab = true` を既定にしたことで入り、上の +97.7 の中に含まれる。ライブラリの形式の変更（`sab_transfer`）については、ADR-0205 の
帰結が「既定（静止した標的）と `sab` を書かない `"free-gas"` の結果は変わらない」と書き、確かめたのは 1C の 0 と 5 MWd/kg の static、
heavy-cap、heavy-fast3 の `swap_absorption` の出力が ADR-0204 の記録のログと文字まで同じこと。

## 候補（確かめていない）

1C と 2A 以外の問題は、ADR-0209 が測っていない。同じ既定の切り替えを候補とするが、効きは分けていない。

- 既定の切り替えの効きが問題（温度、組成、可燃性吸収体、制御棒、集合体の形状）で違う大きさ: VIII.0 の変化は +34〜+154 pcm で、1D（燃料 1200 K）の +153.5 pcm が大きい理由は分けていない。examples/depletion の p2 の 2F など 7 本は、同じ切り替えの前後で xlc の k∞ が
  −39.1〜+58.2 pcm 動いた（examples/depletion/verification.md の「集合体の燃焼」の「計算」、ブランチ chore/rerun-after-adr-0226）
- VII.0 の変化: 1C は +79.8 pcm で、VIII.0 の +97.7 と違う。ADR-0209 は VII.0 のライブラリで測っていない。評価（S(α,β) を含む）の違いで
  同じ切り替えの効きが変わることが候補
- 05714cea から 4cfe4f58 の間の、既定の切り替えでない変更: 上の「確かめていないこと」のとおり、1C と 2A でも分けていない

既に挙がっている見方（open の issue と ADR の「候補」「確かめていない」の行を「散乱核」「自己遮蔽」「KENO」で検索した）:

- #0131: VII.0 の xlc − KENO-VI が VIII.0 の xlc − OpenMC より約 +130 pcm 大きい。候補は VII.0 のファイルの扱い、AMPX と NJOY の処理の差、
  KENO-VI のモデルの解釈、565 K の S(α,β) の温度の補正。今回の動きは VIII.0 と VII.0 の両方で同じ向きで、VIII.0 の 1C と 2A は ADR-0209 が記録した動きなので、
  #0131 の差とは別（VII.0 は上の候補）
- #0040: 565 K の集合体の k∞ が KENO の 2A より約 170 pcm 大きい（前の版の値）
- #0005: k の差のうち、散乱の移行（共鳴の範囲の自己遮蔽を含む）が未解明で最も大きな候補
- #0001: 散乱核を熱中性子の範囲まで扱えば、#0152 の候補の 1 つ（H の S(α,β)）も扱う

## 完了の条件

1C と 2A 以外の問題（VIII.0 の 1B、1D、2B〜2K と 2Q、VII.0 の全問題、P3）で、ADR-0209 の既定の切り替えとそのほかの変更の効きを、今の bin で ADR-0209 の
前の既定（`ufg-scattering = "static"`、`ufg-transfer = "infinite-dilution"`）を書いた入力を解いて分け、VIII.0 の xlc − KENO-VI と xlc − OpenMC
（examples/vera/verification.md の「OpenMC との比較」）の今の値で、どちらに近づいた（離れた）かを記録に書く。

この計算は後回しにする。1C と 2A の変化は ADR-0209 が OpenMC に対して測って受け入れた差の再掲で、精度の作業より機能の拡充を進めるユーザーの
方針に従う（2026-11-04 17:24。ユーザーは、精度は重要だが、KENO に対する差が OpenMC に対しても同じ傾向なら、
精度を詰めるより機能の拡充を先に進めるよう求めた）。
OpenMC でも同じ傾向かを確かめたのは 1C と 2A だけで、ほかの問題（1D の +153.5 pcm を含む）は確かめていない。ほかの問題の差を詰める
作業が来たら、その問題で上の入力を解く。

## 経緯

- 2026-11-04: ADR-0226 の後の再実行（ブランチ chore/rerun-after-adr-0226）で見つけた
- 2026-11-04: 1C と 2A の変化が ADR-0209 の記録した既定の切り替えの前後の動きと一致することを確かめ、ADR-0205 を別の候補から外した（計算なし）。ほかの問題の切り分けは後回し
