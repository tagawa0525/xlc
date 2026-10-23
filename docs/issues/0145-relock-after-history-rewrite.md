---
status: open
labels: [verification]
opened: 2026-10-22
closed:
milestone: M4
blocked-by: []
resolution:
---

# 履歴の書き換えと opt-level の変更で全例の dvc.lock が deps と食い違う

## 問題

履歴の書き換え（#0141。辞書で語を直す。「ステージ」→「ステージ」など）が、DVC のステージの deps になっている toml と py の
コメントを変えたので、書き換えの前に作った `dvc.lock` が、どの例でも deps と合わない。例: `examples/depletion/p2/xs.toml` の md5 は
lock の `a657fbde` に対して今は `63cb0331`（書き換えの前の版 `~/.local/share/history-rewrite/xlc-2026-10-21/original.git` の HEAD は
`a657fbde`。差はコメントだけ）。出力（depletion.msgpack など）は lock と一致している。

加えて、開発のプロファイルの opt-level を 1 にした変更（bcb1ae2dd。Cargo.toml）で、リポジトリの直下の `build` の deps が変わる。
`build` を作り直すと（2026-10-22、main fcb2289c）`bin/` の 4 つの実行ファイルは lock と同じバイトで、変わるのは直下の
`dvc.lock` の Cargo.toml の md5 だけ（opt-level の変更はリリースのプロファイルに及ばない）。

このため、下流のステージを `dvc repro` で流すと、上流（NJOY、xs、チェーン）も再実行される。2026-10-22 の `dvc status` で deps が
changed のステージは c5g7、core、depletion、depletion/p2、m0、mox、reflector、vera、vera/endfb70、vera/endfb70/openmc、vera/endfb70/p3、
vera/p3、vera/p5 にある。ただし、そのうち `bin/` の差には、書き換えより前に計算を実行せずにマージしたコードの変更（ADR-0130 の
「帰結」のとおり）が混ざっており、すべてが書き換えによるものではない。lock と再実行の扱いは ADR-0130 と ADR-0159 に従い、この issue はそれらを変えない。

### 参照解の deps の .py を文書だけ直しても changed になる（2026-10-23）

ADR-0183 に合わせたコミット（3baba57c0。main に入っている）が examples/depletion/openmc_depletion.py の docstring とコメントだけを直した（112.5 → 112.47 eV と
引く ADR の番号）。このファイルは `openmc-1C` と `openmc-1C-extract` の deps（examples/depletion/dvc.yaml）なので、OpenMC の燃焼（3517e4666 は 1.5 時間 × 2、acfe57b01 の本文は 5 時間と書く。どちらも確かめていない）と
その extract が changed になった。χ の担当が compare-1C を目標に `dvc repro` を流すと OpenMC が始まり、止めた（Slurm 454。acfe57b01 の本文による。sacct では 11 分 55 秒で止めている）。
`dvc commit` で付け直した（main の acfe57b01（f804cba1c でマージ）。lock の古い md5 は 3baba57c0 の直前のファイルで、出力は同じコードから
できている。3517e4666 が上の書き換えのあとに行ったのと同じ手段）。

- /dvc-relock、/git-merge、CLAUDE.md に、`dvc commit` で付け直す手順はない。/dvc-relock が拠る ADR-0130 は、走らせなかったステージの
  changed を「そのままでよい」とする
- 3517e4666 と acfe57b01 の `dvc commit` は、ADR-0159 の案 B（ADR-0143 の帰結と ADR-0130 の案 L も同じ理由で退けた形。出力を変えないかの
  判断が人に残り、lock が記録した deps のコードで実行していない出力を指す）を、その場の判断で使っている。ADR-0159 の移行は同じ手段を
  1 回だけ、変更の直前に changed でなかったステージに限って使った。この 2 回は、その限り方を ADR で決めずに使った
- Rust は ADR-0159 の `build`（ADR-0160 で strip）が実行ファイルを deps にするので、実行ファイルのバイトが変わらない変更（tests.rs だけの
  変更など）は計算のステージを動かさない。ただしコードの行をずらすコメントは、実行ファイルを変えうる（ADR-0159 が測った）。
  deps がソースそのものである .py には、バイトの同一で止める仕組みがなく、文書だけの変更が出力の再計算の候補になる

関係する ADR: なし

## 例ごとの扱い

ユーザーの決定（2026-10-22）: 再実行が 2 時間程度で済む例は再実行して lock を作り直す。それより重い例は再実行せず、ADR-0130 の
とおり結果が要るときに再実行する。

- リポジトリの直下（`build`）: 作り直して relock する（数十秒）
- examples/depletion/p2: ライブラリ（xs の 90 核種、NJOY の openmc-neutron/thermal/index）と `chain`（examples/depletion/dvc.yaml）を
  再実行し、集合体の燃焼の lattice@2A と compare@2A を流して relock する（このブランチ）。2F 以降のジョブは、別の作業として流す（2026-10-22 の作業の分け方。このブランチでは流さない）
- examples/depletion の `openmc-1C`（1.5 時間 × 2）、examples/vera の OpenMC: 再実行しない。結果が要るときに再実行する
- ほかの例: 再実行しない。結果が要るときに再実行する

## 再発を防ぐ仕組みの案（決めていない）

- 履歴の書き換えの手順（nixfiles の history-rewrite）に、書き換えのあとで各例の `dvc.lock` の deps の md5 を付け直す段を入れる。
  出力を作らずに、deps の md5 だけを新しいファイルの md5 に置き換える（中身の差がコメントだけであることを、書き換えの前後の
  ファイルの比較で確かめてから）。手段は `dvc commit <ステージ>`（deps と outs の今の md5 を、再実行せずに lock に書く）。
  ADR-0159 の移行が同じ手段を 1 回だけ使い、対象を「変更の直前の main で `dvc status` が changed を示さなかったステージ」に
  限った。同じ限り方（書き換えの直前の main、控えの original.git の HEAD で changed でなかったステージ）なら、書き換えより前から
  changed だったステージ（`bin/` の差など）を取り込まない。p2 のライブラリは、再実行（1:44:15）の前後で出力の md5 が同じだった
  ので、この形で取り込んでも同じ lock になった（書き換えの直前に changed でなかったかは確かめていない）。ADR-0130 の案 L と
  ADR-0143 が `dvc commit` を退けた理由（記録した deps のコードで実行していない出力を lock が指す）に当たるかは、コメントだけの
  差を「同じコード」とみなせるかで決まる
- 辞書の書き換えの対象から、DVC の deps になっているファイル（toml、py）を外す。外すと語が残り、辞書の hook が次の変更で止める
- deps のコメントの変更を DVC に無視させる（deps をコメントを除いた正規形の hash にする）。DVC の仕組みにはなく、ステージの
  deps を正規形のファイルに替える前処理が要る
- 何もしない（ADR-0130 のとおり、結果が要るときに再実行する）。書き換えは稀で、再実行の費用はそのときに払う

### 参照解の deps の .py の文書だけの変更への案（決めていない。決めるときは ADR）

- (a) 参照解の deps の .py を文書だけ直すマージの前に、`dvc commit` で付け直す手順を /dvc-relock に足す。ADR-0130 の案 L、ADR-0143 の帰結、ADR-0159 の案 B が
  退けた形との整合（コメントと docstring だけの差を「同じコード」とみなす条件と、その確かめ方）を ADR で決める要がある
- (b) ADR-0159 と ADR-0160 の延長で、.py の deps を、docstring とコメントを除いた AST の正規形を出す段の出力に替える。構造で防ぐ案で、
  正規形を作る段（数秒）と、deps を替える dvc.yaml の変更が要る
- (c) 当面は CLAUDE.md の「`dvc repro` を投げる前に走るステージを数える」規則（docs/0150-estimate-dvc-repro-stages。
  main に入っている）で防ぐ。走るステージに参照解が入っていれば、投げずに原因を調べる。事故の費用は防ぐが、changed のままの
  ステージは残る

## 完了の条件

- p2 のライブラリと lattice@2A、compare@2A が relock されている（2F 以降のジョブは別の作業で、この issue の条件に含めない）
- 再発を防ぐ仕組みを決める（採るなら ADR に書き、history-rewrite に入れる）か、採らないと決める

## 経緯

- 2026-10-21: #0141 の書き換えのセッションは、history-rewrite の手順の 6（「xlc なら `/dvc-relock` に従って lock を直す」）に
  沿って `dvc status` で changed を見たが、/dvc-relock の後半の「走らせなかったステージの changed はそのままでよい（ADR-0130）」に
  従って何もせず、ユーザーへの報告の「やらなかったこと」に書いた。issue は開かず、下流を流すときの再実行の費用も見積もらなかった
- 2026-10-22: p2 の集合体の燃焼を流す前に `dvc status` で、下流を流すと上流（NJOY、xs）も再実行されると分かり、この issue を開いた
- 2026-10-22: p2 のライブラリ（xs、NJOY の openmc-neutron/thermal/index）と `chain` を再実行した（Slurm のジョブ 417、排他、1:44:15）。
  p2 の 185 個の出力と `chain` の出力の md5 はすべて以前の lock と同じで、コメントの変更と作り直した実行ファイルはライブラリを変えない。
  続けて #0135 のマージ（main 7790c871。crates/lattice を変える）に rebase して `build` を作り直し（変わった実行ファイルは
  `bin/lattice` だけ。xs、chain、core は同じバイトで、ライブラリの lock はそのまま有効）、lattice@2A（ジョブ 421、排他、1:40:36）と
  compare@2A を流した。直下の `dvc.lock` の crates と `bin/lattice` の md5 の変化は、この rebase によるもので、書き換えと opt-level の
  変更によるものではない
- 2026-10-22: examples/depletion の xs（xs-range、xs-nuclide、xs）、NJOY（openmc-neutron、openmc-thermal、openmc-index）、openmc-chain、
  openmc-1C、openmc-1C-extract は、再実行せずに `dvc commit` で付け直した（ブランチ feat/yield-cutoff-default、3517e4666）。changed の deps は
  どれも、書き換えの前の控え（original.git の HEAD）の md5 が lock と同じで、差はコメントの語だけだった。`bin/xs` は p2 と同じ
  52102d9e から a316fea4 への変化で、p2 の再実行で出力が変わらなかったものと同じ。付け直しで出力の md5 は変わらない（変わったのは deps の
  md5 と size だけ）。lattice-1C と compare-1C は同じブランチで再実行した。lattice-2A は付け直さず、changed のまま残す
- 2026-10-23: 3baba57c0（ADR-0183 の docstring の変更）で openmc-1C とその extract が changed になり、compare-1C の `dvc repro` が OpenMC を始めた
  （Slurm 454）。acfe57b01 で `dvc commit` により付け直した。「問題」の節と「案」の節に追記した
- 2026-10-24: 語の規則を直したコミット（ed95d3d44。#0142）が `examples/depletion/p2/xs.toml` と `examples/depletion/chain.toml` のコメントだけを
  直し、p2 の `xs-range`、`xs` と examples/depletion の `chain` が再び changed になった（同じ問題の 2 度目）。p2 の `xs-nuclide` の 90 本も
  `bin/xs` の md5 の変化（d3dc3807 → c6bdf397。ADR-0188 のビットを保つベクトル化とスレッドの変更）で changed で、`compare@<job>` を目標に
  `dvc repro` を流すとライブラリの作り直しが走る。高速化の後の p2 の再実行（ブランチ chore/p2-rerun-after-speedups）で、投げる前の
  `dvc status --with-deps` で見つけ、ライブラリと `chain` を再実行して出力の md5 を lock と比べた（結果は examples/depletion/verification.md の
  「集合体の燃焼」の「計算」）。同じ時点で、p2 の NJOY のステージ（openmc-neutron、openmc-thermal、openmc-index）と、
  examples/depletion の xs-range、xs、NJOY のステージ、openmc-chain も changed のまま残る（このブランチでは再実行していない。p2 の集合体の
  燃焼の目標の上流には入らない）
