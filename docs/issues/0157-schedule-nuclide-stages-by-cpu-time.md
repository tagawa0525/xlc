---
status: open
labels: [xs, performance]
opened: 2026-10-23
closed:
milestone:
blocked-by: []
resolution:
---

# 核種ごとのステージを、前回の CPU 時間の長い順にスレッドの数を決めて同時に実行する

## 問題

`dvc repro` はステージを 1 つずつ実行するので、核種ごとの xs のステージ（ADR-0138）は 1 核種ずつ走る。核種の中の並列でコアを埋めても
（ADR-0189 の案 D）、examples/depletion の 65 核種で核種のプロセスの和は 1,154 秒で、CPU 時間の和の 1/16（710 秒）の 1.6 倍かかる
（Slurm 479、xlc a7c8546c、2026-10-23）。changed の核種のステージを、前回の CPU 時間の長い順に、スレッドの数（`xs --threads`）を
ceil(CPU 時間 / (CPU 時間の和 / 16)) に決めて合計 16 スレッドの中で同時に実行すると、同じジョブの試行で 720 秒だった（ADR-0189 の案 C。
ユーザーの案、2026-10-23）。

足すものは、DVC の外でステージの cmd を実行する runner、核種ごとの CPU 時間の記録（ENDF の大きさ × 温度の数とはほとんど相関しない。ADR-0189）、終わったステージを
lock に取り込む `dvc commit` の運用。`dvc commit` は ADR-0159 の案 B と #0145 が扱う、lock が実行していない deps を指すおそれに関わる（runner は
今の deps で実行するので当たらないが、実行中に deps が変わったときの扱いが要る）。

関係する ADR: ADR-0189

## 完了の条件

runner を足すか足さないかを ADR で決める。足すなら、出力の md5 が `dvc repro` と同じこと、`dvc commit` の運用の条件、DVC を通した
経過時間を確かめる。

## 経緯

- 2026-10-23: ADR-0189 で案 D を採り、案 C をここに分けた
