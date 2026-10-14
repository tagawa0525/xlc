---
name: openmc-reference
description: VERA の OpenMC の参照解を作る手順。ライブラリ、入力の TOML、run と extract、バランス（balance）、スレッド数、出力の置き場、Slurm の名前と見積もり、worktree での DVC の出力の戻し方。
---

# VERA の OpenMC の参照解

xs と同じ ENDF/B-VIII.0 から、手元の NJOY で作ったデータで OpenMC を動かす（ADR-0124）。手順の全体とツールの版は
docs/environment.md の「OpenMC」「VERA の参照解」。計算の条件（核種と温度、粒子数、世代、非アクティブバッチ、乱数のシード、
スレッド数）は入力の TOML に書き、スクリプトに既定値はない。足りないキーと知らないキーで止まる（ADR-0130）。

コマンドはリポジトリ（worktree）の直下で実行する。`references/openmc` に `cd` して `../..` で指すと、worktree でなく
main のスクリプトが使われる（`references` は main への symlink）。

## worktree の準備

`references` が main への symlink になっていること（`ls -l references`。なければ docs/environment.md の `ln -s`）。
DVC の出力（`examples/vera/openmc-library/` など）は worktree では空なので、`dvc.lock` が指す版を戻す:

```sh
dvc checkout examples/vera/dvc.yaml
```

キャッシュ（`references/dvc-cache`）に出力がなければ戻せない。そのときは、下の `dvc repro` で作る。

## ライブラリ

`examples/vera/openmc-library/cross_sections.xml` と、`neutron/`、`thermal/` の HDF5 は、`examples/vera/dvc.yaml` の
`openmc-neutron@<核種>`、`openmc-thermal@<散乱則>`、`openmc-index` の段階が作る。温度は引数でなく
`examples/vera/openmc-library.toml` の `temperatures` に書く。核種や温度を足したら、その項目だけ `dvc repro` が作る。

- 565 K（2A、2Q）は、燃料、冷却材、被覆管の核種を 565 K でも作り、H in H2O の S(α,β) は 550 K と 573.6 K の内挿にする
  （ADR-0080）。OpenMC の温度の内挿のモードは、温度を 1 つしか持たない核種が 1 つでもあると全体を最も近い温度の
  モードに戻すので、内挿で解く問題の核種には 2 つ以上の温度を書く
- `references/openmc/endfb80`（600 K だけ）と `endfb80-565` は、段階にする前に手で作ったもので、段階が作る HDF5 と
  バイト単位で一致した（ADR-0130）。新しい計算は `examples/vera/openmc-library/` を使う

## 段階にしてある問題（1B）

段階は `openmc-1B-model`（スクリプトと入力 → model.xml。数秒）と `openmc-1B`（model.xml とライブラリ → statepoint）に分かれる
（ADR-0136）。スクリプトを直してモデルが変わらなければ model.xml の md5 が同じで、`openmc-1B` は再実行しない。model.xml が変わる
直し（タリーを足すなど）は、`openmc-1B` を再実行する。`model.xml` の生成だけを確かめるときは、
`dvc repro -s examples/vera/dvc.yaml:openmc-1B-model`（上流は `references-shell` だけで、数秒。手元でよい）。

`dvc repro` は `default` の devShell で実行する（`references` の devShell の中では cargo と dvc がない）。ライブラリの
63 段階を含み、全体で 40〜60 分かかるので、全体を Slurm で流す:

```sh
~/.claude/scripts/slurm-run.sh -t 90 vera-reference-stages 'dvc repro examples/vera/dvc.yaml'
```

終わったら `dvc.lock` をコミットする。

## 段階にしていない問題（1C〜2K、2Q）

`examples/vera/p1/openmc-1B.toml` と同じ形の入力を作る（置き場は `p2/openmc-<問題>.toml`）。キーは `problem`、
`cross-sections`（入力のファイルからの相対）、`particles`、`batches`、`inactive`、`seed`、`threads`。既存の参照解の条件:

| 問題       | particles | batches | inactive | seed | threads |
| ---------- | --------- | ------- | -------- | ---- | ------- |
| 1B〜1D     | 1000000   | 250     | 50       | 1    | 16      |
| 2A〜2K、2Q | 1000000   | 400     | 100      | 1    | 16      |

```sh
nix develop .#references -c python3 examples/vera/openmc_model.py model examples/vera/p2/openmc-2A.toml references/openmc/vera-2A-model
nix develop .#references -c python3 examples/vera/openmc_model.py run examples/vera/p2/openmc-2A.toml \
  references/openmc/vera-2A-model references/openmc/vera-2A
nix develop .#references -c python3 examples/vera/openmc_model.py extract examples/vera/p2/openmc-2A.toml \
  references/openmc/vera-2A references/openmc/vera-2A-results
cp references/openmc/vera-2A-results/pins.msgpack examples/vera/p2/openmc-2A.msgpack
```

- `model` は `model.xml` を書く（ライブラリのパスとスレッド数を含まない。ADR-0136）。`run` はそれを解いて作業ディレクトリに
  statepoint を残し、`extract` が k を `k.json` に、2 で始まる問題のピン出力を `pins.msgpack` に書く。抽出だけを変えるときは
  `run` を再実行しない（ADR-0159）
- `p2/openmc-<問題>.msgpack` は git で追跡しない（ADR-0121）。worktree では main のものをコピーして使う
- 2E-empty は 2E の Pyrex の棒を抜いたもの。問題の名前は `openmc_model.py` の docstring にある一覧
- k は `k.json` の `k` と標準偏差 `k-std`（xlc との差は `compare_k.py <lattice の出力> <k.json>`）。記録の書き方は、スキル verification-record

## 釣り合い、多群、均質な媒質

位置引数のスクリプトは、スレッド数を入力でなく環境変数で受ける（`OMP_NUM_THREADS=16`）。省略すると全コアを使う。

```sh
nix develop .#references -c sh -c 'OMP_NUM_THREADS=16 python3 examples/vera/openmc_balance.py \
  examples/vera/openmc-library/cross_sections.xml 1B references/openmc/balance-1B examples/vera/vera.msgpack 1000000 250 50 1'
```

引数は `<cross_sections.xml> <問題> <作業ディレクトリ> <xlc のライブラリ> [粒子数/世代] [世代数] [非アクティブバッチ数] [乱数のシード]`。
`openmc_homogeneous.py`（`<cross_sections.xml> <媒質 (TOML)> <作業ディレクトリ> <xlc のライブラリ> …`）も同じ形。
`openmc_mg.py` は `<多群定数 (msgpack)> <問題> <作業ディレクトリ> …` で、手順は docs/environment.md の末尾。
`openmc_balance.py` の statepoint は世代ごとに残り（1B で 200 本）大きい。使い終わった作業ディレクトリは、
`balance.msgpack` を残して statepoint を消してよいか、ユーザーに確かめる。

## AWS の自己完結な実行（Problem 3、5-2D）

ADR-0155、0168、0171。手順は docs/environment.md の「自己完結な実行」。要点:

- `aws_run.sh --type auto`: 候補と速度は `examples/vera/aws_types.toml`。スポットの vCPU のクォータの残りに収める
- 起動し直しは `--job <ジョブ>`。`--statepoint-interval` は 3A、3B で 50、5B、5C-2D で 25
- `--retry-failed` は、S3 の `runs/<name>/<ジョブ>/log` で原因を見て、直ってから使う
- S3 には最新の statepoint が 120 秒後に写る（それより前の中断では失う）
- `aws login` の期限は約 12 時間。切れると手元の監視と起動し直しが止まる（VM は動き続ける）。切れたらユーザーに `aws login` を頼む
- 監視と起動し直しは `examples/vera/aws_watch.sh <name>`（docs/environment.md に使い方）

## Slurm の名前と見積もり

スキル slurm-run に従い、`-t` は分の数、名前は `<ツール>-<題材>`。過去の実績（16 スレッド。`slurm-history.sh --by-name`、
verification.md の「参照解」）:

| 計算                                | 実績              | 名前の例                |
| ----------------------------------- | ----------------- | ----------------------- |
| 1B の `dvc repro` 全体              | 38〜57 分（3 回） | `vera-reference-stages` |
| 1B の OpenMC（statepoint の total） | 17.8 分           | `openmc-1B`             |
| 2A〜2K、2Q の OpenMC                | 31〜50 分         | `openmc-2A`             |
| `openmc_balance.py` 1B（250 世代）  | 59 分             | `openmc-balance-1B`     |
| `openmc_balance.py` 2K（400 世代）  | 97 分             | `openmc-balance-2K`     |
| 均質な無限媒質                      | 5〜21 分          | `openmc-homog-h-u5`     |

実績の 1.3 倍ほどを `-t` に書く。長い計算で時間を記録しないもの（参照解を出し直す、方向を決める試行）は、`-c 16` で
排他にせず流してよい。そのときの時間は記録に書かない（スキル slurm-run）。
