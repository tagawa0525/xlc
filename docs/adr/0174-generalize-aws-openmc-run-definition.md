---
status: accepted
date: 2026-10-15
requires: [ADR-0155, ADR-0168, ADR-0171]
supersedes: []
superseded-by:
issues: ["#0132"]
---

# ADR-0174: Generalize the AWS OpenMC run with a run definition file

## 背景

2026-10-14 の方針で、数時間かかる OpenMC の解析は AWS（スポット、`--type auto`。ADR-0171）で流し、xlc 自身の計算は手元の Slurm で流す。
AWS のツール（`aws_run.sh`、`aws_user_data.sh`、`aws_common.sh`、`aws_fetch.sh`、`aws_status.sh`、`aws_watch.sh`。ADR-0155、0168）は
`examples/vera` の Problem 3 と 5-2D（`openmc-3?`、`openmc-5?-2D`）の参照解だけを前提に書かれている。燃焼の参照解
（`examples/depletion/openmc_depletion.py`。M4-5。feat/lattice-depletion-rates の枝で作業中で、main にはまだない）は同じ OpenMC のクロージャを使うが、ツールに載らない。

この ADR は、ツールが vera に固定している箇所と燃焼の参照解の要件を調べ（「調査の結果」）、どう一般化するかを決める案を比べる。
案 A を選び、実装は docs/aws-generic-openmc の枝で行った（「実装の結果」）。

判断の軸:

- 既存の挙動（VERA の m3 の取り出し、`--job` による起動し直し、statepoint の同期。ADR-0168）を壊さない
- 後方互換を作らない（利用者がいない）。形式を変えたら古い形式を残さない
- 変更の大きさと、同期や再開の論理を 2 か所で直し続けないこと

## 調査の結果

### 1. ツールが `examples/vera` に固定している箇所

| 箇所                                                                                                    | 固定の中身                                                                                                                                                                                                                                                        |
| ------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `aws_run.sh` の `tool_files`                                                                            | `src/` に置くツールが `examples/openmc_*.py` 4 個、`examples/vera/openmc_{tallies,model,model_3d,core_2d}.py`、`run_remote.sh` の 9 個の固定の列                                                                                                                  |
| `aws_run.sh` の入力の確認                                                                               | `[[ -d examples/vera ]]`、`examples/vera/openmc-library/cross_sections.xml` の存在とその中の全 `path=` のファイルの存在、入力の `cross-sections` が `examples/vera/openmc-library/cross_sections.xml` を指すこと                                                  |
| `aws_run.sh` の S3 への置き方                                                                           | `examples/vera/openmc-library` を `store/openmc-library/<鍵>/` に sync。置くデータはこのディレクトリ 1 つだけ（燃焼の chain の XML を置く場所がない）                                                                                                             |
| `aws_common.sh` の `job_script`、`aws_run.sh` の確認                                                    | ジョブ名が `openmc-3?` か `openmc-5?-2D` でなければ止まる。名前から `openmc_model_3d.py` か `openmc_core_2d.py` を選ぶ（`run_remote.sh` の case と二重）                                                                                                          |
| `aws_user_data.sh`                                                                                      | VM でライブラリを `examples/vera/openmc-library` に取る。流すのは `bash examples/vera/run_remote.sh [--restart] [--statepoint-interval N] <入力>` の固定                                                                                                          |
| `run_remote.sh`                                                                                         | `results/<ジョブ>/run` に run、続けて extract（`$PYTHON <script> extract <入力> <run> <out>`）、`done` に終了コード。スレッド数は物理コア数 ÷ 入力の数                                                                                                            |
| `aws_fetch.sh`                                                                                          | `[[ -d examples/vera ]]`。`--verify` は `job_script` で extract のツールを選び、`$dest/src/<script> extract <入力> <run> <check>` を VM の結果とバイト比較                                                                                                        |
| `aws_watch.sh`                                                                                          | 再起動の `--statepoint-interval` は `openmc-3?` で 50、それ以外は 25 の固定                                                                                                                                                                                       |
| `aws_user_data.sh` の `sync_out`、`aws_common.sh` の `job_latest_statepoint`、`aws_run.sh` の起動し直し | 進み具合を `run/statepoint.<世代>.h5` の存在で測る。最新は書き途中でありうるので mtime が 120 秒より古いものだけコピーする（ADR-0168）。再開の有無は S3 の statepoint の有無、「完了なのに done がない」の判定は入力の `batches` と最後の世代の statepoint の比較 |

固定されていないもの: クロージャ（`.#openmc`、`.#openmc-python`）の S3 への置き方と鍵、機種の選択（ADR-0171）、起動の記録（`launches/`）、
`done` と `SHA256SUMS` の流れ、`aws_status.sh`（`INPUTS` から ジョブ名と最新の statepoint を出すだけ）、`aws_terminate.sh`（タグ XlcRun）。
`examples/openmc_{common,input,run,resume}.py` は problem 1〜5 と c5g7、mox が共有するので、vera の外の statepoint を使う実行
（`openmc_run.py` と `openmc_resume.py`）は、入力とツールの置き場の一般化だけで載る見込み（c5g7 と mox は確かめていない）。

結果の使われ方: vera の m3 は `aws_fetch.sh` で `references/openmc/aws/m3/<ジョブ>/` に取り出し、DVC のステージの出力にはせず、比較の段階が
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
  完了した実行（最後の項目が `final_step` の輸送の結果。時間幅 0）を延長する使い方は成り立つが、ステップが保存された直後の中断
  （最後の項目が時間幅のあるステップの先頭の組成だけで、終わりの組成が保存されていない状態）では、先頭の組成を終わりの時刻から
  始めてしまうように読める。また `depletion_results.h5` は `mode="a"` で追記するので、書いている最中の中断で壊れうる。
- **失うのは多くて 1 ステップ**: 1C の 1 リング（粒子 20 万・130 世代、16 スレッド）は全体で 5473 秒、1 ステップ（輸送 2 回）は約 14 分
  （`openmc_simulation_n<k>.h5` の mtime の間隔）。ステップの途中の輸送を再開する必要はなく、ステップ単位で十分。
- **chain とライブラリ**: `openmc-chain.xml` は 27 MB（DVC の出力。`.gitignore` 済み）、`openmc-library/` は 431 MB（vera は 631 MB）で、
  どちらも DVC の出力。ライブラリは NJOY で作る段階（`openmc-neutron@核種` など 90 余り）で、VM には NJOY と ENDF（数 GB）がない。
  vera と同じく手元で作って S3 に置く（`store/` に置いたままで月 $0.02 以下）。VM のクロージャ（`.#openmc`、`.#openmc-python`）は
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

- **案 A: 実行の定義を入力のファイルで渡す**: `examples/<問題>/aws.toml` に、アップロードするツール、データ（ライブラリ、chain）、VM の入口
  のスクリプト、進み具合の種類（`statepoint` か `steps`）を書く。`aws_run.sh` などは vera の固定をやめ、定義を読む。vera の定義は
  今の固定の値をそのまま書く。長所: 同期、再開、機種選択（ADR-0171）、監視、起動の記録を 1 か所で保てる。次の題材
  （c5g7、mox、将来の炉心の燃焼）は定義 1 つで載る。短所: 定義の読み書き（TOML の部分集合を bash で読む。`aws_types.toml` と同じ流儀）と、
  進み具合の種類ごとの分岐が増える。
- **案 B: example ごとに薄い aws の入力（スクリプト）を置く**: `examples/depletion/aws_run.sh` が `aws_common.sh` を source し、
  `tool_files` やライブラリのパスだけを替えて、同じ関数を呼ぶ。長所: 定義の形式を作らない。短所: `aws_user_data.sh` の同期の論理
  （ADR-0168 の mtime、`--delete` の順序、SIGTERM の扱い）が statepoint に書き込まれているので、燃焼の `steps` の進み具合には書き直しが要り、
  2 つの user-data が分岐する。`aws_status/fetch/watch/terminate` も複製か、vera のパスを引数に取る形になる。ADR-0168 のような修正が
  2 か所になる。2 つ目の複製は許容範囲（CLAUDE.md の「重複は 3 回まで」）だが、同期の論理は競合状態が絡み、複製が最も壊れやすい部分。
- **案 C: 燃焼だけ別のツール**: `examples/depletion/aws_*.sh` を新しく書く。vera のツールを触らないので m3 を壊す恐れがない。短所: 案 B の短所に加え、
  S3 の置き場、機種選択、スポットのクォータ、起動の記録、監視を丸ごと複製し、保守コストが 2 倍になる。ステップ単位の再開の論理だけが
  燃焼に固有なのに、固有でない約 1000 行を複製する。

## 決定と理由

**案 A**。同期と再開の論理は競合状態を含み、ADR-0168 のような修正が続いた部分で、複製すると直し漏れが起きる。燃焼に固有なのは
(1) 入口のスクリプトとデータの置き場、(2) 進み具合の種類（ステップ単位の完了の印）で、定義の項目になる。複製を許容できる範囲内の案 B でも
(2) のために user-data の分岐が要るので、最初から定義にしたほうが総量は少ない。

### 定義の形

`examples/<問題>/aws.toml`。`aws_run.sh <名前> <入力>...` は、各入力の親のディレクトリを上へたどって最初に見つかる `aws.toml` を使う
（`--def` は作らない。1 つの実行の入力はすべて同じ定義であること）。定義は `src/` に置き、`launch.env` に `DEFINITION=` で残す。
起動し直しは S3 の `src/` の定義を使う。

```toml
tools = ["examples/openmc_common.py", ...]          # src/ に置く（入口と定義自身と入力は含めない）
entry = "examples/depletion/run_remote.sh"          # VM: bash <entry> [--restart] [--statepoint-interval N] <入力>
progress = "steps"                                  # "statepoint"（今の挙動）か "steps"
library = "examples/depletion/openmc-library"       # cross_sections.xml と全 path= の存在を確かめて S3 に置く
data = ["examples/depletion/openmc-chain.xml"]      # そのほかの入力ファイル（ハッシュで store/data/ に置く。省略できる）
jobs = '^openmc-1C-[a-z0-9]+$'                      # 受け付けるジョブ名（今の job_script の case）
```

契約: 入口は `results/<ジョブ>/run` に実行の状態を、`results/<ジョブ>/done` に終了コードを書く（今の `run_remote.sh` と同じ）。
入口は `--extract-script <入力>` で、その入力の run と extract を行うスクリプトを出し、`aws_fetch.sh --verify` はそのスクリプトの
`extract <入力> <run> <out>` が out に書くファイルを VM の結果とバイト比較する。燃焼の extract は `<out>/<ジョブ>.msgpack` を書く。

案の段階では定義に `extract` のキーを置いていたが、実装で外した。VERA の 3? と 5?-2D は extract のスクリプトが違い、定義に
1 つのパスしか書けないと `job_script` の case（`run_remote.sh` の case との二重）が残るため、入口に問い合わせる形にして、
ジョブ名とスクリプトの対応を `run_remote.sh` の 1 か所にした。

### 燃焼の進み具合（`progress = "steps"`）

OpenMC の再開は完了した実行の延長に向くので（調査の結果 2）、**1 ステップを 1 区間として、区間の終わりごとに完了の印を置く**。
`run/step-<k>/depletion_results.h5` と、区間の最後に置く印 `run/step-<k>/complete`（内容は `<k>/<全ステップ数>`。k は終わった
ステップの数）。

- 区間 k は別のプロセスで、`step-<k-1>/depletion_results.h5` を `step-<k>/` にコピーして `prev_results` に渡し、残りの 1 ステップだけを
  `integrate` する。最後の区間以外は `final_step=False`。次の区間は、保存した反応率から始めず、先頭の輸送をやり直す
  （OpenMC の `_get_bos_data_from_restart` を差し替える）。輸送の回数は分けない実行と同じ 13 回（6 点）。
- 中断された区間は、`complete` のない `step-<k>/` を捨て、最後の `complete` のある区間から始め直す。壊れうるのは書いていた最中の
  ファイルだけで、元の `step-<k-1>/` は触らない。失うのは多くて 1 ステップ（1C で約 14 分）。
- 最後の区間の `step-<n>/depletion_results.h5` が結果で、区間ごとの `transport-seconds.json` は前の区間までの時間を含む。extract は
  最後の区間だけを読む。
- S3 への同期は、`complete` のある `step-<k>/` だけ。区間のファイルを写し終えてから、`complete` を別の呼び出しでコピーする（印だけが先に
  着いて、欠けた区間を再開が読まないように）。VM の再開は、S3 から `complete` のある最後の区間だけを戻す。statepoint の mtime の
  規則（ADR-0168）は使わない。
- 最後の区間まで終わって `done` がないジョブ（extract の前の中断）は、`steps` では再開できる（区間を流さずに extract を流す）。
  `aws_run.sh` の起動し直しの「最後の世代の statepoint があって done がない」（statepoint は世代が残っておらず再開できない）は、
  `progress` ごとの関数（`job_latest_progress`、`job_progress_finished`）に分けた。`--statepoint-interval` は `steps` では
  受け付けず 2 で止まる。`aws_watch.sh` は `steps` で渡さない。
- 古い区間は消さない。完了のときに消す案は、同期が一覧を取ったあとに区間が消えると、含めたはずの区間の同期が失敗する（statepoint の
  ADR-0168 と同じ競合）ので採らなかった。ただし再開した VM は最後の区間だけを戻すので、成功のときの `--delete` の同期で、前の VM が
  S3 に置いた古い区間は S3 から消える。

### 検討して捨てた方式（実装の途中）

- **保存した反応率から始める再開**（OpenMC の既定。`prev_results` を渡すだけ）: 最初の半ステップの核分裂収率が、直前の輸送の
  スペクトルで決まる `chain.fission_yields`（`CoupledOperator.__call__` が設定し、結果に保存されない）でなく、新しい `Operator` の
  既定の収率で解かれる（OpenMC 0.15.2 のソース `deplete/openmc_operator.py`、`abc.py` の読み）。粒子 1000・世代 20・捨てる 10 の
  小さい入力（burnup 0.1、0.5 MWd/kg）で試すと、2 点目は分けない実行と同じ k（1.141633）だが、3 点目が 1.129814 に対して
  1.121824 だった。差 0.008 は統計誤差（0.010）の内で、系統差とは言えない。収率の違いを測って確かめてはいない。原理的な理由
  （収率が実行の履歴に依存し、保存されない）で、先頭の輸送をやり直す形にした。
- **完了のときに古い区間を消す**: 上のとおり、同期との競合のため。

### 変更の大きさの見積もり

実装前の見積もり（実際は下の「実装の結果」）:

| 変更                                                                                                | ファイル数 | 行数（追加/削除）の目安   |
| --------------------------------------------------------------------------------------------------- | ---------- | ------------------------- |
| `aws_run.sh`（定義を読む、入力の確認の一般化、`data` の置き方）                                     | 1          | +70 / -45                 |
| `aws_user_data.sh` と `test_aws_user_data.sh`（ライブラリと data の取得、入口、`steps` の同期）     | 2          | +90 / -25                 |
| `aws_common.sh`（定義の読み込み、`job_script` の削除、進み具合の関数の分岐）                        | 1          | +60 / -15                 |
| `aws_fetch.sh`、`aws_status.sh`、`aws_watch.sh`                                                     | 3          | +30 / -15                 |
| `examples/vera/aws.toml`（今の固定の値）、`run_remote.sh` の入口の契約への調整                      | 2          | +25 / -5                  |
| `examples/depletion/aws.toml`、入口 `run_remote.sh`、再開 `openmc_depletion_resume.py` とそのテスト | 4          | +300（テスト 120 を含む） |
| `docs/environment.md`、`.claude/skills/openmc-reference/SKILL.md`、`flake.nix` の check             | 3          | +60 / -20                 |

ツールは `examples/vera/` から `examples/aws/` に、定義の実装と別のコミットで先に移した（ロジックの変更なし。参照の約 30 か所のパスを直した）。

既存の挙動を守る手順: (1) 先に `aws.toml` の読み込みだけを入れ、vera の定義が今の固定の値と一致することを `--dry-run` の出力
（user-data を含む）の差がないことで確かめる、(2) `steps` を足す。後方互換は作らず、`tool_files` の固定と `job_script` は消した。

### 実装の結果

- `examples/aws/`（`aws_run.sh`、`aws_common.sh`、`aws_user_data.sh`、`aws_fetch.sh`、`aws_status.sh`、`aws_watch.sh`、
  `test_aws_user_data.sh`、新しい `test_aws_definition.sh`）。定義の読み込みは `aws_common.sh` の `def_parse`（TOML の部分集合。標準の
  bash と grep、sed だけ）と `def_find`。flake の check は `aws-user-data` と `aws-definition`（リポジトリの 3 つの定義を読む）。
- `examples/vera/aws.toml`、`examples/aws/aws-check/aws.toml`（確認用の入力は `examples/aws/` に移ったので、親のディレクトリに置く定義の
  写し）、`examples/depletion/aws.toml`、`examples/depletion/run_remote.sh`、`examples/depletion/openmc_depletion_resume.py`。
- `data` の S3 の置き場は `store/data/<鍵>/<ファイル名>`。launch.env に `DEFINITION=` と、`data` があれば `DATA=`（`パス=鍵`）が増える。
- 新しい実行の `src/` には、ツールと入口のほかに定義自身を置く。起動し直しはツールと入口だけを置き直す。
- `examples/depletion/openmc_depletion.py` は main にまだない（feat/lattice-depletion-rates）ので、`aws.toml` の `tools` にあるが、
  このブランチのテストは OpenMC を使わない部分（区間の置き場の規則）だけ。main に入る順は、feat/lattice-depletion-rates のあと。

### 帰結

- 今の `aws_run.sh m3 ...` の呼び出しは変わらない（入力から定義が決まる）。`launch.env` に `DEFINITION=` が増えるので、この ADR より前の実行
  （m3）の起動し直しは、`DEFINITION` がなければ止まる。m3 は完了後に取り出して S3 から消した。読み込みの互換は作らない。
- 時間の測定を含む燃焼の記録（`transport-seconds`、`cram-seconds`）は AWS の結果から使わない。比べる時間は手元の Slurm の値。
- 未決: AWS の結果を `examples/depletion` の比較の段階の入力にするときの dvc.yaml の形（vera の 3A と同じく
  `references/openmc/aws/` のパスを deps にする案）。区間の分割でステップごとに `Operator` を作り直す費用は、小さい入力
  （1 コア、粒子 1000）で全体 25〜40 秒のうち初期化が数秒の見込みだが、1C の大きさでは測っていない。再開は 1 区間が約 14 分かかるので、
  中断が続くと進まないことがある（`aws_watch.sh` が起動し直すが、`steps` に区間内の保存はない）。

### 確認

実装の前に確かめるとした項目を、次のとおり確かめた（2026-10-15。数値を出した xlc のコミットは、この枝の `feat(depletion)` のコミット。
手元、`nice 19`。重い計算はしていない。粒子 2 万の比較はコミット前の `openmc_depletion_resume.py` で、コミットとの違いは時間の記録、
extract、スレッドの環境変数だけ）:

1. 小さい入力で、区間に分けた実行が分けない実行と同じ k と数密度を出すこと。
   - 粒子 2 万・世代 40・捨てる 20、burnup 0.1、0.5 MWd/kg、1 リング（`examples/depletion/aws-check/medium.toml`、比べるツールは同じディレクトリの `compare_split.py`、乱数のシード 1。OpenMC が
     `OMP_NUM_THREADS` を読まずに約 10 スレッドで走った）で、分けない実行
     （`openmc_depletion.py run`）と、区間に分けた実行（`openmc_depletion_resume.py run`）を比べた。k は 0 点目（燃焼度 0）が一致
     （1.17024 ± 0.00177、同じ輸送）。1 点目は 1.13752 ± 0.00140 に対して 1.13341 ± 0.00181（差 −0.00411、1.8σ）、2 点目は
     1.12419 ± 0.00134 に対して 1.12523 ± 0.00143（差 +0.00104、0.5σ）。参考に、乱数のシードだけ 2 に変えた分けない実行は、種 1 の
     分けない実行に対して 1 点目が −0.00253（1.1σ）、2 点目が −0.00213（1.2σ）で、区間に分けた差と同じ大きさ。数密度（最後の点、
     数密度が 1e-9 /(b cm) を超える 151 核種）の相対差は、中央値が 1.3e-4（種を変えた分けない実行は 2.5e-5）、最大が 5.7e-3
     （種を変えた分けない実行は 1.1e-2）で、U-235 は +5.2e-6（−4.5e-6）、Pu-239 は +2.4e-3（−2.8e-6）、U-238 は −1.1e-6（+1.2e-7）。
     3 点の k と数密度は、種を変えた実行との違いと同じ範囲で、統計誤差で区別できない。ただし粒子が少なく（統計誤差 約 150 pcm）、
     3 点では 1.8σ の 1 点を棄却できる検出力はない（系統差が小さければ見えない）。輸送の回数は両方 5 回（2 点で 2N + 1）。
   - 再開: 粒子 1000 の入力で区間に分けた実行の `step-2/` を消して `--restart` すると、k と全核種の数密度が元の実行と一致した
     （`openmc_depletion_resume.py run` で完全一致。VM の入口 `run_remote.sh --threads 2` では 1e-13 の丸めの違い）。VM と同じく最後の区間だけを残して
     再開しても、extract が書けた（最後の区間に時間が累積）。
   - 区間の保存直後に殺した実行（(b) として元の確認項目にあった）は、区間が `complete` を最後に置き、なければ捨てるので、
     殺された実行の `depletion_results.h5` から続けることはない。この前提を殺して確かめてはいない。
2. `.#openmc-python` で `openmc_depletion.py` と `openmc_depletion_resume.py` の import が通る（`msgpack` 1.2.1、`h5py` 3.16.0、
   OpenMC 0.15.2）。`aws.toml` の `tools` と入口だけを置いた作業ディレクトリで、`run_remote.sh` が小さい入力を最後まで流した（done が 0、
   k は手元と一致）。
3. 1C の 1 リングを c7a.4xlarge で流し、手元の Slurm の結果と k が統計誤差の範囲で合うこと: **未実施**（実機。AWS を起動していない）。
   `aws_run.sh` は `--dry-run` だけを確かめた: VERA の m3、m5 の入力で、user-data がバイト単位で実装前と同じ（`launch.env` の
   `DEFINITION=` の 1 行だけが増える）。確認用の入力は移動で入力のパスが変わるため user-data の入力のパスだけが違う。燃焼の入力
   の `--dry-run` は通り、user-data は steps 用になる。S3 への同期の絞り込み（`--exclude`、`--include` の順序）は、実際の S3 では
   試していない。
