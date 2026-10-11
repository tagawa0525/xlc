---
name: verification-record
description: examples/*/verification.md に測定の結果を書く・直す手順。数値を出したコミットの取り方、rebase 後の付け替え、k の差の内訳の道具、記録に書く試行の扱い、Slurm の時間の注意。
---

# 検証の記録を書く

結果を `examples/<問題>/verification.md` に書く。書き方の規則（今の姿だけを書く、結果ごとにコミットと日付を添える、
測り直すのは計算を変えたときだけ）は CLAUDE.md と ADR-0117、ADR-0122 にある。ここには、その手順の細部を書く。

## 数値を出す前

- 記録に書く計算は、流す前に main に rebase する（スキル topic-start の「記録に書く計算の前に main に rebase する」）。
  計算のあとの rebase は hash の付け替えと relock を増やす
- 数値を出すコードと入力は、コミットしてから走らせる。`git status --short` が空であることを確かめる（追跡しない
  ファイルは無視される）。作業ツリーが汚れていると、記録のコミットを checkout しても同じ結果が出ない
- 記録に書くコミットは `git rev-parse --short=8 HEAD`（8 桁。既存の記録と同じ）。日付は走らせた日
- 段階（`dvc.yaml`）で走らせたら `dvc.lock` をコミットする（ADR-0123）。記録には、走らせたコミットと、ブランチの上なら
  分岐点（`git merge-base HEAD main`）を書く。rebase で deps が変わっても走らせ直さず、分岐点から復元する（ADR-0130 の
  「帰結」。復元は分岐点にそこまでのブランチのコミットを当て直し、`dvc status` の md5 で確かめる）
- 結果を作る側が複数あるときは、それぞれのコミットを書く（例: 「OpenMC は 62b800c0、比べたのは ee3cd6ec、2026-10-11」）

## 計算は Slurm で、時間は排他のジョブのものだけ

数値を記録に書く計算は、`~/.claude/scripts/slurm-run.sh` で流す（スキル slurm-run）。CLAUDE.md が示す
`-t <見積もり> <名前> '<コマンド>'` の形で、既定の排他で流す。CLAUDE.md は直接の `sbatch` を使わないと定め、hook が止めるとも書く（hook の所在は確かめていない）。

- 記録に書く時間は、排他のジョブの `TIME:` 行の `elapsed`、または OpenMC が statepoint に書く計算時間（total。既存の
  「参照解」はこちら）
- `-c` で排他にせず流したジョブの時間は書かない。ログの `TIME:` 行に `shared (-c N): 時間は当てにならない` が付いている。
  数値（k など）は書いてよい

## 書き方

- 書き換える節は、古い測定、誤りの注記、置き換えた設定での結果を消して、今の結果だけにする。ほかの節の数値を
  打ち消す注記は書かない。経緯は git に任せる（`git show <コミット>:<ファイル>` で辿れる）
- 再現のコマンドは、リポジトリの直下から走らせる形でコードブロックに書く（「k の差の内訳」の節の形）。python3 は
  `nix develop .#references -c` の中、cargo と dvc は `default` の devShell の中
- 確かめていない原因は「候補」「確かめていない」と書く（CLAUDE.md）。保存則（k と 1/A など）が閉じていることを、
  比べる前に書いておく

## rebase のあと

記録のコミットの hash は、`remap.commitRefs` が true のリポジトリでは、rebase のあとに post-rewrite hook が（rebase の途中の amend では呼ばず、rebase の終わりに一度だけ呼ぶ）
`commit-refs remap` で付け替える。コミットはしないので、`git diff` で付け替えた記録を確かめてコミットする（件名の例:
`docs: point the records at the rebased commits`）。測り直さない（ADR-0122）。

- hook が失敗したときは、出力の最後に「直し方」が出る。`$(git rev-parse --absolute-git-dir)/commit-refs-remap.input` に
  対応表が残るので、`commit-refs remap <そのファイル>` で呼び直す。`ERROR:` は、何も書き換えずに止まったときに出る。出力の理由を読んでから直す。
  理由を読んで、対応表か文書のどちらが誤りかを確かめてから直し、呼び直す
- `~/.claude/scripts/remap-commit-refs.sh` は使わない（`commit-refs` に置き換わった古い写しで、今も残っている）。`commit-refs` を使う
- filter-repo のような履歴の書き換えは `commit-refs rewrite-history <対応表>`

## k の差の内訳（OpenMC との差を物質とエネルギーの範囲に分ける）

外周が鏡面反射の格子で、k の差を正味の吸収の差に分ける。確かめた問題は 1B と 2K で、2J は止まる。リポジトリの直下で、
この順に走らせる（結果の読み方と式は各スクリプトの docstring と verification.md の「k の差の内訳」）:

```sh
# xlc 側。対称性なしで MOC を解く（k は lattice と 6 桁まで同じ）
cargo run --release -p xlc-lattice --example balance -- examples/vera/p2/assembly-2K.toml xlc-2K.msgpack
# OpenMC 側。16 スレッドで 1B は約 1 時間、2K は約 1.6 時間（Slurm で、スキル openmc-reference）
python3 examples/vera/openmc_balance.py examples/vera/openmc-library/cross_sections.xml 2K references/openmc/balance-2K \
  examples/vera/vera.msgpack 1000000 400 100 1
# 比べる。末尾は範囲の境（eV）。省略すると 15030、4、0.625。2K は 27 群を分けるので 9.877 を足す
python3 examples/vera/compare_balance.py xlc-2K.msgpack references/openmc/balance-2K/balance.msgpack \
  examples/vera/vera.msgpack 15030 9.877 4 0.625
```

- 群ごと・核種ごとに分けるのは `compare_groups.py`（`openmc_balance.py` の世代ごとの `out_scatter` と `within_scatter`
  が要る）。均質な無限媒質は `homogeneous_medium.py` と `openmc_homogeneous.py`
- 移行の差の k への影響は `swap_transfer` の例で測る。水の 14〜27 群の散乱行列の行を OpenMC の率 / 中性子束に置き換える:

  ```sh
  cargo run --release -p xlc-lattice --example swap_transfer -- examples/vera/p2/assembly-2K.toml \
    references/openmc/balance-2K/balance.msgpack water 14 27
  ```

  群は 1 始まりで高いエネルギーから数え、`<最初の群> <最後の群>` を含む
- 釣り合いが閉じること（OpenMC のすべてのエネルギーの 1/A が k の推定値と統計誤差の中で合う）を、差の分け方を書く前に確かめる
- 内訳の「残り」は原因ごとの分け方ではない（2 次の項、範囲の外、間接の影響を含む）。原因とは書かない

## 記録に書く試行

試行の入力は `trials/`（無視される）に作る。記録に数値を書く試行だけ、次の順で扱う（ADR-0122）:

1. `git add -f trials/<入力>` でその入力だけをコミットし、記録にそのコミットの hash と入力のパスを文字で書く
   （パスは `trials/` から消えるので、リンクにしない）
2. マージの前に、記録に書かない試行を `git rm -r` で消すコミットを入れる（マージは `--no-ff` なので、消す前のコミットは
   main から辿れる）。残ったままだと、マージのコミットの前の hook が走らせる `cargo test`（`main_reports_trials_left_by_a_merge`）が止める
3. 繰り返し比べる道具と、採用した入力は `git mv trials/<入力> examples/<問題>/` で移し、記録のパスと hash を移した先と
   移したコミットに書き換える。計算（コード、入力、環境）を変えていなければ測り直さない
