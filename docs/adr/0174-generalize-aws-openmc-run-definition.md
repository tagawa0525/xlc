---
status: proposed
date: 2026-10-14
requires: [ADR-0155, ADR-0168, ADR-0171]
supersedes: []
superseded-by:
issues: []
---

# ADR-0174: Generalize the AWS OpenMC run with a run definition file

## 背景

2026-10-14 の方針で、数時間かかる OpenMC の解析は AWS（スポット、`--type auto`。ADR-0171）で流し、xlc 自身の計算は手元の Slurm で流す。
AWS の道具（`aws_run.sh`、`aws_user_data.sh`、`aws_common.sh`、`aws_fetch.sh`、`aws_status.sh`、`aws_watch.sh`。ADR-0155、0168）は
`examples/vera` の Problem 3 と 5-2D（`openmc-3?`、`openmc-5?-2D`）の参照解だけを前提に書かれている。燃焼の参照解
（`examples/depletion/openmc_depletion.py`。M4-5。ADR-0173 の枝で作業中）は同じ OpenMC の閉包を使うが、道具に載らない。

この ADR は、道具が vera に固定している箇所と燃焼の参照解の要件を調べ（「調査の結果」）、どう一般化するかを決める案を比べる。
実装はこの ADR を承認してから別のブランチで行う。

判断の軸:

- 既存の挙動（VERA の m3 の取り出し、`--job` による起動し直し、statepoint の同期。ADR-0168）を壊さない
- 後方互換を作らない（利用者がいない）。形式を変えたら古い形式を残さない
- 変更の大きさと、同期や再開の論理を 2 か所で直し続けないこと

## 調査の結果

### 1. 道具が `examples/vera` に固定している箇所

| 箇所                                                                                                    | 固定の中身                                                                                                                                                                                                                                                  |
| ------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `aws_run.sh` の `tool_files`                                                                            | `src/` に置く道具が `examples/openmc_*.py` 4 個、`examples/vera/openmc_{tallies,model,model_3d,core_2d}.py`、`run_remote.sh` の 9 個の固定の列                                                                                                              |
| `aws_run.sh` の入力の確認                                                                               | `[[ -d examples/vera ]]`、`examples/vera/openmc-library/cross_sections.xml` の存在とその中の全 `path=` のファイルの存在、入力の `cross-sections` が `examples/vera/openmc-library/cross_sections.xml` を指すこと                                            |
| `aws_run.sh` の S3 への置き方                                                                           | `examples/vera/openmc-library` を `store/openmc-library/<鍵>/` に sync。置くデータはこのディレクトリ 1 つだけ（燃焼の chain の XML を置く場所がない）                                                                                                       |
| `aws_common.sh` の `job_script`、`aws_run.sh` の確認                                                    | ジョブ名が `openmc-3?` か `openmc-5?-2D` でなければ止まる。名前から `openmc_model_3d.py` か `openmc_core_2d.py` を選ぶ（`run_remote.sh` の case と二重）                                                                                                    |
| `aws_user_data.sh`                                                                                      | VM でライブラリを `examples/vera/openmc-library` に取る。流すのは `bash examples/vera/run_remote.sh [--restart] [--statepoint-interval N] <入力>` の固定                                                                                                    |
| `run_remote.sh`                                                                                         | `results/<ジョブ>/run` に run、続けて extract（`$PYTHON <script> extract <入力> <run> <out>`）、`done` に終了コード。スレッド数は物理コア数 ÷ 入力の数                                                                                                      |
| `aws_fetch.sh`                                                                                          | `[[ -d examples/vera ]]`。`--verify` は `job_script` で extract の道具を選び、`$dest/src/<script> extract <入力> <run> <check>` を VM の結果とバイト比較                                                                                                    |
| `aws_watch.sh`                                                                                          | 再起動の `--statepoint-interval` は `openmc-3?` で 50、それ以外は 25 の固定                                                                                                                                                                                 |
| `aws_user_data.sh` の `sync_out`、`aws_common.sh` の `job_latest_statepoint`、`aws_run.sh` の起動し直し | 進み具合を `run/statepoint.<世代>.h5` の存在で測る。最新は書き途中でありうるので mtime が 120 秒より古いものだけ写す（ADR-0168）。再開の有無は S3 の statepoint の有無、「完了なのに done がない」の判定は入力の `batches` と最後の世代の statepoint の比較 |

固定されていないもの: 閉包（`.#openmc`、`.#openmc-python`）の S3 への置き方と鍵、機種の選択（ADR-0171）、起動の記録（`launches/`）、
`done` と `SHA256SUMS` の流れ、`aws_status.sh`（`INPUTS` から ジョブ名と最新の statepoint を出すだけ）、`aws_terminate.sh`（タグ XlcRun）。
`examples/openmc_{common,input,run,resume}.py` は problem 1〜5 と c5g7、mox が共有するので、vera の外の statepoint を使う走行
（`openmc_run.py` と `openmc_resume.py`）は、入力と道具の置き場の一般化だけで載る見込み（c5g7 と mox は確かめていない）。

結果の使われ方: vera の m3 は `aws_fetch.sh` で `references/openmc/aws/m3/<ジョブ>/` に取り出し、DVC の段階の出力にはせず、比較の段階が
そのパスを deps にする（`examples/vera/verification.md` の 3A）。

### 2. 燃焼の参照解に固有の要件

`openmc_depletion.py run <toml> <outdir>` は、`openmc.deplete.CECMIntegrator`（CE/CM、燃焼度の点 6 個）を `final_step=True` で 1 回の
`integrate` で通す。ステップごとに輸送を 2 回（先頭と中点）、最後に 1 回で、6 点なら 13 回（`transport-seconds.json`）。

- **statepoint の形がちがう**: OpenMC の輸送の途中の statepoint は設定しておらず（`statepoint-interval` なし）、書かれるのはステップの先頭の
  `openmc_simulation_n<k>.h5`（24 KB。`write_bos_data`）と、各輸送の最後の世代の `statepoint.<batches>.h5`（20 MB。輸送のたびに同じ名前で上書き）だけ。どちらも再開には使えない
  （OpenMC の depletion は輸送の途中から再開できない）。進み具合と再開の状態は `depletion_results.h5`（`StepResult.save` がステップごとに
  追記する。1C の 1 リングで 3.9 MB）にある。現在の `statepoint.*.h5` を数える同期と、`batches` との比較は燃焼に使えない。
- **OpenMC の再開は途中からの再開に向かない**（OpenMC 0.15.2 のソース `deplete/abc.py` の読み。**実機で確かめていない**）:
  `prev_results` を渡すと `integrate` は `Results[-1]` の `data[0]`（先頭の組成）と反応率、`time[-1]` から始まり、`_i_res = len(prev) - 1`。
  完了した走行（最後の項目が `final_step` の輸送の結果。時間幅 0）を延長する使い方は成り立つが、ステップが保存された直後の中断
  （最後の項目が時間幅のあるステップの先頭の組成だけで、終わりの組成が保存されていない状態）では、先頭の組成を終わりの時刻から
  始めてしまうように読める。また `depletion_results.h5` は `mode="a"` で追記するので、書いている最中の中断で壊れうる。
- **失うのは多くて 1 ステップ**: 1C の 1 リング（粒子 20 万・130 世代、16 スレッド）は全体で 5473 秒、1 ステップ（輸送 2 回）は約 14 分
  （`openmc_simulation_n<k>.h5` の mtime の間隔）。ステップの途中の輸送を再開する必要はなく、ステップ単位で十分。
- **chain とライブラリ**: `openmc-chain.xml` は 27 MB（DVC の出力。`.gitignore` 済み）、`openmc-library/` は 431 MB（vera は 631 MB）で、
  どちらも DVC の出力。ライブラリは NJOY で作る段階（`openmc-neutron@核種` など 90 余り）で、VM には NJOY と ENDF（数 GB）がない。
  vera と同じく手元で作って S3 に置く（`store/` に置いたままで月 $0.02 以下）。VM の閉包（`.#openmc`、`.#openmc-python`）は
  `msgpack`、`h5py`、`numpy` を含み、`openmc_depletion.py` の import は満たす見込み（`referencesPython` の依存。実行して確かめていない）。
- **時間の測定が混ざる**: 燃焼の参照解は輸送と CRAM の時間（`transport-seconds`、`cram-seconds`）を記録し、xlc の時間との比較に使う。
  AWS の機種はまちまちで、中断で別の機種に変わるので、この時間は比べる値として使えない（CLAUDE.md の「時間を記録する測定」は
  Slurm）。AWS で流した結果は時間を記録に使わないと決める必要がある。
- **`openmc_depletion.py` を変えると重い段階が古くなる**: `openmc-1C` の deps に入っているので、ここに再開を足すと `dvc.lock` が古くなり、
  再び 1.5 時間以上の計算が要る。vera が `openmc_resume.py` を別に置いた（`openmc_run.py` を変えず problem 1、2 を古くしないため）のと同じく、
  再開は別のファイルに置き、`openmc_depletion.py` の `build`、`TimedOperator`、`read` を import する。
- **DVC との接続は vera と同じく切る**: `openmc-1C@rings1` の `outs` は AWS で流した結果の置き場にならない。結果は
  `references/openmc/aws/<名前>/<ジョブ>/` に取り出し、比較の段階がそのパスを deps にする。

## 検討した案

- **案 A: 実行の定義を入力のファイルで渡す**: `examples/<問題>/aws.toml` に、アップロードする道具、データ（ライブラリ、chain）、VM の入口
  のスクリプト、進み具合の種類（`statepoint` か `steps`）を書く。`aws_run.sh` などは vera の固定をやめ、定義を読む。vera の定義は
  今の固定の値をそのまま書く。長所: 同期、再開、機種選択（ADR-0171）、監視、起動の記録を 1 か所で保てる。次の題材
  （c5g7、mox、将来の炉心の燃焼）は定義 1 つで載る。短所: 定義の読み書き（TOML の部分集合を bash で読む。`aws_types.toml` と同じ流儀）と、
  進み具合の種類ごとの分岐が増える。
- **案 B: example ごとに薄い aws の入力（スクリプト）を置く**: `examples/depletion/aws_run.sh` が `aws_common.sh` を source し、
  `tool_files` やライブラリのパスだけを替えて、同じ関数を呼ぶ。長所: 定義の形式を作らない。短所: `aws_user_data.sh` の同期の論理
  （ADR-0168 の mtime、`--delete` の順序、SIGTERM の扱い）が statepoint に書き込まれているので、燃焼の `steps` の進み具合には書き直しが要り、
  2 つの user-data が分岐する。`aws_status/fetch/watch/terminate` も複製か、vera のパスを引数に取る形になる。ADR-0168 のような修正が
  2 か所になる。2 つ目の複製は許容範囲（CLAUDE.md の「重複は 3 回まで」）だが、同期の論理は競合状態が絡み、複製が最も壊れやすい部分。
- **案 C: 燃焼だけ別の道具**: `examples/depletion/aws_*.sh` を新しく書く。vera の道具を触らないので m3 を壊す恐れがない。短所: 案 B の短所に加え、
  S3 の置き場、機種選択、スポットのクォータ、起動の記録、監視を丸ごと複製し、保守コストが 2 倍になる。ステップ単位の再開の論理だけが
  燃焼に固有なのに、固有でない約 1000 行を複製する。

## 決定と理由

**案 A**。同期と再開の論理は競合状態を含み、ADR-0168 のような修正が続いた部分で、複製すると直し漏れが起きる。燃焼に固有なのは
(1) 入口のスクリプトとデータの置き場、(2) 進み具合の種類（ステップ単位の完了の印）で、定義の項目になる。複製を許容できる範囲内の案 B でも
(2) のために user-data の分岐が要るので、最初から定義にしたほうが総量は少ない。

### 定義の形（案）

`examples/<問題>/aws.toml`。`aws_run.sh <名前> <入力>...` の入力が `examples/<問題>/` の下にあれば、その定義を使う（`--def` は作らない。
1 つの実行の入力はすべて同じ <問題> の下にあること）。定義は `src/` に置き、`launch.env` に `DEFINITION=` で残す。

```toml
tools = ["examples/openmc_common.py", ...]          # src/ に置く（入力は含めない）
entry = "examples/depletion/run_remote.sh"          # VM: bash <entry> [--restart] [--threads N] [--statepoint-interval N] <入力>
extract = "examples/depletion/openmc_depletion_resume.py"  # fetch --verify: $PYTHON <extract> extract <入力> <run> <out> が out にファイルを書く
progress = "steps"                                  # "statepoint"（今の挙動）か "steps"
library = "examples/depletion/openmc-library"       # cross_sections.xml と全 path= の存在を確かめて S3 に置く
data = ["examples/depletion/openmc-chain.xml"]      # そのほかの入力ファイル（ハッシュで store/data/ に置く）
jobs = '^openmc-1C-[a-z0-9]+$'                      # 受け付けるジョブ名（今の job_script の case）
```

契約: 入口は `results/<ジョブ>/run` に走行の状態を、`results/<ジョブ>/done` に終了コードを書く（今の `run_remote.sh` と同じ）。
`extract` の出力は out の下のファイルで、`aws_fetch.sh --verify` が VM の結果とバイト比較する。燃焼の `openmc_depletion.py extract` は
出力を 1 ファイル（msgpack）で受けるので、再開用の別ファイルに、`<out>/<ジョブ>.msgpack` を書く extract を置く。

### 燃焼の進み具合（`progress = "steps"`）の案

OpenMC の再開は完了した走行の延長に向くので（調査の結果 2）、**1 ステップを 1 区間として、区間の終わりごとに完了の印を置く**。
`run/step-<k>/depletion_results.h5` と、`integrate` が戻ったあとの空の印 `run/step-<k>/complete`（k は終わったステップの数）。

- 区間 k+1 は、`step-<k>/depletion_results.h5` を `step-<k+1>/` に写してから `Results` で読み、`prev_results` に渡して、
  `continue_timesteps=True` のとき全ステップの列のうち残りを 1 ステップだけ `integrate(final_step=True)` する。`final_step=True` の最後の
  輸送は次のステップの先頭の輸送そのもの（組成も出力密度も同じ）なので、輸送の回数は増えない（13 回のまま）。
- 中断された区間は、`complete` のない `step-<k+1>/` を捨て、最後の `complete` のある区間から始め直す。壊れうるのは書いていた最中のファイルだけで、
  元の `step-<k>/` は触らない。失うのは多くて 1 ステップ（1C で約 14 分）。
- S3 への同期は `complete` のある `step-<k>/` だけ（`depletion_results.h5` は 4〜数十 MB、`openmc_simulation_n<k>.h5` は 24 KB）。
  statepoint の mtime の規則（ADR-0168）は使わない。最新だけを S3 に残す（古い区間は完了のとき捨てる）。
- 完了の判定は、最後の区間の `complete`（k = 燃焼度の点の数）。`aws_run.sh` の起動し直しの「最後の世代の statepoint があって done がない」は、
  `progress` ごとに最後の印を見る関数に分ける。`--statepoint-interval` は `steps` では受け付けず 2 で止まる。`aws_watch.sh` は `steps` で渡さない。
- 区間に分けた走行が分けない走行（今の `openmc_depletion.py run`）と同じ k と数密度を出すことは、実装の最初の PR で小さい入力
  （粒子 1000、燃焼度 2 点）で確かめる。確かめるまで AWS の結果を参照解として使わない。

### 変更の大きさの見積もり

| 変更                                                                                                | ファイル数 | 行数（追加/削除）の目安   |
| --------------------------------------------------------------------------------------------------- | ---------- | ------------------------- |
| `aws_run.sh`（定義を読む、入力の確認の一般化、`data` の置き方）                                     | 1          | +70 / -45                 |
| `aws_user_data.sh` と `test_aws_user_data.sh`（ライブラリと data の取得、入口、`steps` の同期）     | 2          | +90 / -25                 |
| `aws_common.sh`（定義の読み込み、`job_script` の削除、進み具合の関数の分岐）                        | 1          | +60 / -15                 |
| `aws_fetch.sh`、`aws_status.sh`、`aws_watch.sh`                                                     | 3          | +30 / -15                 |
| `examples/vera/aws.toml`（今の固定の値）、`run_remote.sh` の入口の契約への調整                      | 2          | +25 / -5                  |
| `examples/depletion/aws.toml`、入口 `run_remote.sh`、再開 `openmc_depletion_resume.py` とそのテスト | 4          | +300（テスト 120 を含む） |
| `docs/environment.md`、`.claude/skills/openmc-reference/SKILL.md`、`flake.nix` の check             | 3          | +60 / -20                 |

合計 16 ファイル前後、+630 / -125 行ほど。道具を `examples/vera/` から `examples/aws/` に移すかは別の判断で、移すと参照するドキュメント、
スキル、`flake.nix` の check の約 30 か所のパス変更が加わる（ロジックの変更はない）。燃焼の道具が `examples/vera/aws_run.sh` を呼ぶのは
紛らわしいので、この実装と別のコミットで移すことを勧める（案 A の実装の前に移す）。

既存の挙動を守る手順: (1) 先に `aws.toml` の読み込みだけを入れ、vera の定義が今の固定の値と一致することを `--dry-run` の出力
（user-data を含む）の差がないことで確かめる、(2) `steps` を足す。後方互換は作らず、`tool_files` の固定と `job_script` は消す。

### 帰結

- 今の `aws_run.sh m3 ...` の呼び出しは変わらない（入力から定義が決まる）。`launch.env` に `DEFINITION=` が増えるので、この ADR より前の実行
  （m3）の起動し直しは、`DEFINITION` がなければ止まる。m3 は完了後に取り出して消す前提で、読み込みの互換は作らない。
- 時間の測定を含む燃焼の記録（`transport-seconds`、`cram-seconds`）は AWS の結果から使わない。比べる時間は手元の Slurm の値。
- 未決（実装のときに決める）: AWS の結果を `examples/depletion` の比較の段階の入力にするときの dvc.yaml の形（vera の 3A と同じく
  `references/openmc/aws/` のパスを deps にする案）。区間の分割で ステップごとに `Operator` を作り直す（`openmc.lib` の初期化が毎回要る）
  費用は測っていない。

### 確認

実装の前に確かめること（この ADR では確かめていない）:

1. 小さい入力（粒子 1000、燃焼度 2 点）で、OpenMC 0.15.2 の `prev_results` + `continue_timesteps` が、(a) 完了した走行の延長で
   分けない走行と同じ k と数密度を出すか、(b) 区間の保存直後に殺した走行の `depletion_results.h5` から続けると何が起きるか
   （調査の結果 2 の読みの確認）。
2. `.#openmc-python` で `openmc_depletion.py` の import（`msgpack`、`openmc.deplete`）が通ること。
3. 1C の 1 リングを c7a.4xlarge で流し、手元の Slurm の結果と k が統計誤差の範囲で合うこと。
