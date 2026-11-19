---
type: Issue
status: open
labels: [dvc, depletion]
opened: 2026-11-07
closed:
milestone:
blocked-by: []
resolution:
---

# DVC が核種ごとの表の見出しの多い xs.toml を見出しの数の 2 乗の時間で読み、ステージの 1 操作に数十秒かかる

## 問題

- DVC 3.67.1 は、`dvc.yaml` の `vars` とステージの `params` の TOML を tomlkit 0.15.0 で読む。tomlkit は、`[nuclide.X]` の形の表の見出しを
  足すたびに親の表を deepcopy し、読む時間が見出しの数の 2 乗で延びる。examples/depletion/p2/xs.toml を 1 回読む時間は、見出し 90 で
  0.21 秒、555 で 9.11 秒、同じ 555 核種を `[nuclide]` のインラインの表にすると 0.07 秒（2026-11-07 に手元で tomlkit.parse を測った）
- DVC の 1 つの操作は同じ TOML を何度も読む。555 の見出しの xs.toml で、1 つのステージの `dvc status` は 34.5 秒、cProfile では 108 秒のうち
  tomlkit の読みが 11 回で 103.9 秒だった（2026-11-07、手元）
- Slurm 724（ブランチ feat/0153-kinetics-table-data、2026-11-07）は、`stage_runner.py` の上流の `dvc status` に 1,116 秒、`build` と `xs-range` の
  `dvc repro` に 7,097 秒、xs-nuclide の commit に 3 ステージで 124 秒、16 ステージで 660 秒と 682 秒かかった（1 コア）。xs そのものは 16 本で
  124 秒以内。この速さでは 555 ステージの commit だけで 6 時間を超えるので、見込みの外れたジョブとして止めた
- インラインの表にした後も、xs-nuclide の commit が 1 コアで batch の時間のほとんどを占める。Slurm 726（2026-11-07）の 16 本の batch は、
  xs が 1 本 1 秒ほどで、batch の間隔（runner は commit の間ほかの行を出さないので、commit を含む）は 2 本目から 41.3 秒、42.6、43.6、… 54.0 秒と
  batch ごとに約 1 秒ずつ延びた（runner のログの `committed` の時刻）。親が ps で見た 1 つの batch は 52 秒のうち `dvc commit -f` の 1 プロセス
  （99 %）が約 47 秒で、1 ステージ約 3 秒。間隔が commit したステージの数とともに延びるので、DVC が呼び出しごとに 555 ステージの dvc.yaml と
  伸びていく dvc.lock を読み書きしていることが候補（確かめていない）。直す案の候補は、commit を最後に 1 回にする、batch を大きくする
- Slurm 887（p2 の 545 核種の `stage_runner.py`、排他、xlc d6cc90d96、2026-11-15）の 2 時間 3 分 9 秒のうち、16 コアの計算は約 27 分だけで、
  上流の確認 547 秒、`xs-range` と核種の投入の前の処理 1,317 秒（中身は確かめていない）、`dvc commit` が 1 回約 100 秒（まとめた核種の数に
  よらず、64 回で約 1 時間 40 分。ログは `~/github/slurm-logs/xlc-p2-library-defaults-887.out`）。後に続く p2 の
  `dvc repro examples/depletion/p2/dvc.yaml:compare@<履歴>`（Slurm 888〜895）も、計算の前の同じ確認で 10〜20 分かかる見込み（#0288 が上の 547 秒の確認から見込んだ値で、測っていない）。
  CLAUDE.md の方針の問題の大きさに釣り合わない時間（#0288 から寄せた）
- 2A のライブラリの放射化の核種（ADR-0238）は、インラインの表にした（xs.toml の説明）。同じ形の大きな入力（ほかの例の核種を足すとき）で
  また起きるのを止める検査はない

関係する ADR: ADR-0184、ADR-0193

## 完了の条件

- 大きな TOML の vars を持つ例で、DVC の 1 操作の時間が見出しの数で 2 乗に延びない形を決め、必要なら検査を足す（tomlkit を直す、vars を
  JSON か YAML にする、表の見出しを数える検査など）
- runner の固定費（上流の確認、核種の投入の前の処理、`dvc commit` の回数と 1 回の時間）を、小さい入力でステージと deps の数を変えて測り、
  計算より短くする形（commit をまとめるか投入と回収から切り離す。ADR-0193 の「commit の前後で deps の指紋を比べる」の制約と合わせる）を決める

## 経緯

- 2026-11-07: ADR-0238 で 2A のライブラリに 465 の放射化の核種を足したとき、Slurm 724 が 1 コアの DVC の処理で見積もりを外れて分かった
- 2026-11-10: ADR-0239 で p2 のライブラリの温度の点を 2 から 13 にした作り直し（Slurm 727、区分 record。割り込まれうるジョブ）では、
  xs.toml の表の見出しは 90 で Slurm 726 と同じだが、16 本をまとめた commit の直前の間隔（`committed` の時刻の差）は 30 回で
  99.9〜106.3 秒（726 は 20 回で 41.3〜73.3 秒）だった。16 本より少ない commit を含むすべての間隔は、727 が 66 回で 7.7〜336.6 秒、
  726 が 35 回で 6.5〜104.7 秒。`xs-range` を始めるまで 550.8 秒（726 は 71.4 秒）。核種の並列の実行（ログの `budget=16` の行から
  `runner elapsed=` の行まで）の 4,916.6 秒と、CPU の和を 16 で割った 1,826 秒の差の約 3,090 秒は、xs の外の時間（commit など）と、
  重い核種が終わるまで 16 の枠が埋まらない偏りの分を含みうる。runner は commit の間に終わった xs を commit の後に `done` として
  記録するので、核種ごとの xs の終わりの時刻はログから読めず、内訳は確かめていない（ADR-0239 の帰結）。xs.toml は 48,579 から
  80,286 バイトになった。延びた原因（TOML の大きさ、dvc.lock の大きさ、出力のハッシュ）は確かめていない
- 2026-11-10: 上の項のレビュー（2 周目）の読み。726 と 727 はどちらも一時停止されていない（`sacct -j 726,727 -X -o Suspended` が
  00:00:00）。runner（examples/stage_runner.py の `execute`）は、空いたスレッドに次を投入した後に `finish(done)` の `dvc commit` を
  同期で呼ぶので、commit の間に終わった xs は回収されず、次も投入されない。レビューの集計（親は数え直していない）では、727 の
  runner が commit の中にいた時間（commit の直前の行から `committed` の行まで）は 66 回の和で 3,657.0 秒、`done` の
  elapsed × threads の和は 16 × 4,916.6 秒の 97.9% で、枠が空いていた偏りの分は経過時間で約 104 秒以下。軽い核種の xs は数秒で
  終わるので、後半は 16 の枠の大半が commit を待つ。直す候補: commit を投入と回収から切り離す（確かめていない）
- 2026-11-19: #0291 のブランチ（chore/0291-terms-in-example-inputs）で、全体の `scripts/count-changed-stages.py`（中は `dvc status --json`）が
  手元で 15 分の timeout（`timeout 900`）に収まらず、数を取れなかった。ステージ 1 本ずつの `dvc status <yaml>:<ステージ>` は返った（時間は測っていない）
- 2026-11-19: DVC の tomlkit だけを 0.15.2 にした（8d5b23da3。flake.nix の `dvcTomlkit`）。2 乗の原因は tomlkit 0.15.0 の既知の誤り
  （python-poetry/tomlkit の issue 479）で、0.15.1 が線形にした。DVC 3.67.1 は PyPI の最新で、tomlkit の指定は `>=0.11.1` だけなので、
  古いのは固定した nixpkgs（2026-10-06）の tomlkit だった。p2（xs.toml は 80,273 バイト、見出し 91、656 ステージ）で、1 コア、手元で測った値:
  1 ステージ（`xs-nuclide@U-234` など）の `dvc status` は、直す前（main 0768dee0f の worktree、出力あり）が 4.3 秒、直した後（8d5b23da3 と
  同じ flake.nix のこのブランチの worktree、出力は戻していない）が 1.3〜1.8 秒。p2 の `dvc.yaml` の全体の `dvc status --json` は、
  このブランチの worktree（出力は戻していない）で、直す前の dvc が 590 秒の打ち切りで終わらず、直した後が 171 秒。直す前の 1 ステージの
  cProfile（main 0768dee0f の worktree、出力あり）は 12.1 秒のうち tomlkit の読みが 3 回で 10.2 秒、`dvc.lock`（837,390 バイト）の読みが
  1.2 秒、ファイルのハッシュが 0.007 秒（出力が変わっていないので、DVC が覚えている md5 を引いた場合の値）
- 2026-11-19: 直した後の p2 の全体の `dvc status` の cProfile（504.5 秒。8d5b23da3 と同じ flake.nix のこのブランチの worktree、出力は
  戻していない）では、`dvc/dependency/param.py` の `read_param_file` が
  636 回で 496.3 秒を占める。DVC は params を持つステージごとに同じ xs.toml を読み直し、読んだ結果を再利用しない（1 回の tomlkit の
  読みが cProfile で 0.66 秒、うち `_validate_out_of_order_table` の deepcopy が全体で約 203 秒。`dpath` でキーを引くのが 66.4 秒）。
  `dvc.lock` の読みは上位に出ない。この測定では出力が無く、DVC はハッシュを取っていないので、ハッシュの時間は分からない
  （出力を書いた直後と、新しい worktree に出力を戻した直後に、ファイルを読み直す時間は測っていない）。測った範囲で残る固定費は、
  ステージの数 × params のファイルを 1 回読む時間である。
  直す候補（確かめていない）: params のファイルを DVC が tomlkit を通さない形式（JSON、YAML）にする、核種ごとの params をやめて
  ファイルの全体を deps にする、DVC に読んだ params のファイルを再利用させる（上流への提案か手元のパッチ）。同じ TOML の標準の読み手
  （tomllib）での読みは 0.02 秒
- 2026-11-19: main 0768dee0f の worktree で、直す前の dvc の `dvc --wait-for-lock status`（p2 の 1 ステージずつ）を 4 本同時に流すと、4.2、8.3、12.2、16.3 秒に 1 本ずつ終わった
  （`--wait-for-lock` が無いと 3 本が「Unable to acquire lock」で終了コード 255）。DVC は worktree ごとの排他を操作の間ずっと持ち、
  ステージの cmd を実行する間だけ外す（`dvc/stage/run.py` の `unlocked_repo`）。DVC の中でステージを並列に実行する機能は上流の
  issue（iterative/dvc の 755。2018-06 から開いたまま）にある。ファイルのハッシュは、ディレクトリの出力の中の 2 MiB を超える
  ファイルが 2 つ以上のときだけ並列で（dvc-data 3.18.3 の `hashfile/build.py`）、ファイル 1 つずつの deps と outs は逐次に取る
- 2026-11-19: `dvc commit` の 1 回約 100 秒（Slurm 887）が tomlkit の版で縮むかは測っていない。`dvc commit -f` は共有のキャッシュの
  run-cache に項目を足すので、キャッシュを分けた使い捨ての複製で測る
