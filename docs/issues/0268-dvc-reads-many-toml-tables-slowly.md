---
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
- 2A のライブラリの放射化の核種（ADR-0238）は、インラインの表にした（xs.toml の説明）。同じ形の大きな入力（ほかの例の核種を足すとき）で
  また起きるのを止める検査はない

関係する ADR: なし

## 完了の条件

- 大きな TOML の vars を持つ例で、DVC の 1 操作の時間が見出しの数で 2 乗に延びない形を決め、必要なら検査を足す（tomlkit を直す、vars を
  JSON か YAML にする、表の見出しを数える検査など）

## 経緯

- 2026-11-07: ADR-0238 で 2A のライブラリに 465 の放射化の核種を足したとき、Slurm 724 が 1 コアの DVC の処理で見積もりを外れて分かった
