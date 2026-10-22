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

履歴の書き換え（#0141。辞書で語を直す。「段階」→「ステージ」など）が、DVC のステージの deps になっている toml と py の
コメントを変えたので、書き換えの前に作った `dvc.lock` が、どの例でも deps と合わない。例: `examples/depletion/p2/xs.toml` の md5 は
lock の `a657fbde` に対して今は `63cb0331`（書き換えの前の版 `~/.local/share/history-rewrite/xlc-2026-10-21/original.git` の HEAD は
`a657fbde`。差はコメントだけ）。出力（depletion.msgpack など）は lock と一致している。

加えて、開発のプロファイルの opt-level を 1 にした変更（bcb1ae2dd。Cargo.toml）で、リポジトリの直下の `build` の deps が変わる。
`build` を作り直すと（2026-10-22、main fcb2289c）`bin/` の 4 つの実行ファイルは lock と同じバイトで、変わるのは直下の
`dvc.lock` の Cargo.toml の md5 だけ（opt-level の変更はリリースのプロファイルに及ばない）。

このため、下流のステージを `dvc repro` で流すと、上流（NJOY、xs、チェーン）も再実行される。2026-10-22 の `dvc status` で deps が
changed のステージは c5g7、core、depletion、depletion/p2、m0、mox、reflector、vera、vera/endfb70、vera/endfb70/openmc、vera/endfb70/p3、
vera/p3、vera/p5 にある。ただし、そのうち `bin/` の差には、書き換えより前に計算を実行せずにマージしたコードの変更（ADR-0130 の
「帰結」のとおり）が混ざっており、すべてが書き換えによるものではない。

関係する ADR: ADR-0130、ADR-0159

## 例ごとの扱い

ユーザーの決定（2026-10-22）: 再実行が 2 時間程度で済む例は再実行して lock を作り直す。それより重い例は再実行せず、ADR-0130 の
とおり結果が要るときに再実行する。

- リポジトリの直下（`build`）: 作り直して relock する（数十秒）
- examples/depletion/p2: ライブラリ（xs の 90 核種、NJOY の openmc-neutron/thermal/index）と `chain`（examples/depletion/dvc.yaml）を
  再実行し、集合体の燃焼（lattice@2A、2F、2N）と、OpenMC の参照解がそろったものの比較を流して relock する（このブランチ）
- examples/depletion の `openmc-1C`（1.5 時間 × 2）、examples/vera の OpenMC: 再実行しない。結果が要るときに再実行する
- ほかの例: 再実行しない。結果が要るときに再実行する

## 再発を防ぐ仕組みの案（決めていない）

- 履歴の書き換えの手順（nixfiles の history-rewrite）に、書き換えのあとで各例の `dvc.lock` の deps の md5 を付け直す段を入れる。
  出力を作らずに、deps の md5 だけを新しいファイルの md5 に置き換える（中身の差がコメントだけであることを、書き換えの前後の
  ファイルの比較で確かめてから）
- 辞書の書き換えの対象から、DVC の deps になっているファイル（toml、py）を外す。外すと語が残り、辞書の hook が次の変更で止める
- deps のコメントの変更を DVC に無視させる（deps をコメントを除いた正規形の hash にする）。DVC の仕組みにはなく、ステージの
  deps を正規形のファイルに替える前処理が要る
- 何もしない（ADR-0130 のとおり、結果が要るときに再実行する）。書き換えは稀で、再実行の費用はそのときに払う

## 完了の条件

- p2 のライブラリと集合体の燃焼のステージが relock され、比較が OpenMC の参照解のそろったものについて流れている
- 再発を防ぐ仕組みを決める（採るなら ADR に書き、history-rewrite に入れる）か、採らないと決める

## 経緯

- 2026-10-22: #0141 の書き換えのあと、p2 の集合体の燃焼を流す前に `dvc status` で分かった
- 2026-10-22: p2 のライブラリ（xs、NJOY の openmc-neutron/thermal/index）と `chain` を再実行した（Slurm のジョブ 417、排他、1:44:15）。
  p2 の 185 個の出力と `chain` の出力の md5 はすべて以前の lock と同じで、コメントの変更と作り直した実行ファイルはライブラリを変えない。
  続けて #0135 のマージ（main 7790c871）に rebase し、lattice@2A（ジョブ 421、排他、1:40:36）と compare@2A を流した。2F 以降は
  #0135 の後の lattice で流す
