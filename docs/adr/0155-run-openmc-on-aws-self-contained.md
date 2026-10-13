---
status: accepted
date: 2026-10-13
requires: [ADR-0151]
supersedes: []
superseded-by:
issues: ["#0106", "#0107", "#0109"]
---

# ADR-0155: Run OpenMC on AWS self-contained with user-data, S3 and self-termination

## 背景

ADR-0151 は OpenMC の参照解を AWS EC2 の c7a のスポットで流すと決めた。試走（docs/environment.md の「試走」、2026-10-13）と
feat/aws-tools のツールは、手元から ssh でクロージャと入力を送り、VM の上で実行し、ssh で結果を戻し、手元から終了する。この形では、
手元の PC が落ちるか AWS に接続できないと、(1) 結果が VM の EBS にしかなく、(2) 終わった VM の課金が止まらない。スポットの
中断（one-time のスポットは回収で終了し、EBS が消える）では、statepoint が残らず 1 本（本番で約 6 時間）をやり直す。

要求（2026-10-13、ユーザー）:

- 手元の PC が落ちても、AWS に接続できなくても、計算が完了し、結果が残り、課金が止まること
- 終わったら結果を S3 に移し、EC2 と EBS は消す。手元は定期的に S3 から取り出し、検証してから S3 の結果を消す

判断の軸:

- VM が手元に頼らずに始めから終わりまで進み、どの失敗でも最後は終了すること
- 中断と失敗で失うものが、最後の同期からの分に収まること
- 権限と常駐する仕組み（料金のかかる常設の資源、別のサービス）が少ないこと

## 検討した案

### 起動のあとの進め方

- **A. user-data（cloud-init）が全部行う**: 起動のツールが user-data を作って `run-instances` に渡す。VM が Nix を入れ、S3 からクロージャ、
  核データ、入力を取り、流し、S3 に同期し、自分で電源を切る。手元は起動のときだけ要る
- **B. 手元から ssh で進める（今のツール）**: 手元が落ちると止まる。要求を満たさない
- **C. AWS Batch、ParallelCluster**: 待ち行列とスポットの再投入を持つが、計算環境、ジョブ定義、コンテナのイメージ（Nix のクロージャを
  コンテナにする手順）が増える。1 回に数本のジョブには重い

### クロージャの受け渡し

- **A. 手元で `nix copy --to file://<一時ディレクトリ>` の binary cache を作り、`aws s3 sync` で `store/nix/` に置く。VM は
  `aws s3 sync` で取り、`nix copy --no-check-sigs --from file://...`（root で走るので信頼される）**: VM は結果の同期に AWS CLI を
  どのみち使うので、転送のツールが 1 つで済む。クロージャ（`.#openmc`、`.#openmc-python`、62 パス、約 800 MB）は zstd で 212 MB になり、
  `nix copy --to file://` は手元で約 2 秒（2026-10-13、xlc のコミット df72096a）。置き場は store のパスの組のハッシュで
  `store/nix/<鍵>/` に分け、同じクロージャは 1 度だけ置く
- **B. `nix copy --to 's3://<bucket>/store?region=us-east-1'`（Nix の S3 binary cache）**: 退けた理由は 2 つ。(1) nix 2.34 の
  `nix help-stores` は URL の形を `s3://bucket-name` とだけ書き、バケットの中の接頭辞（`store/`）を書いていない（確かめていない）。
  (2) 認証: 手元の認証は `aws login` のセッション（`~/.aws/config` の `login_session`）で、Nix の S3 の資格情報の連鎖がこれを読むか、
  VM の Nix がインスタンスプロファイルを読むかを確かめていない。A なら S3 に触るのは AWS CLI だけになる
- **C. VM でビルドする**: ADR-0151 で退けた（再ビルドしない）

AWS CLI は、flake の nixpkgs の `awscli2`（`nix build --inputs-from . nixpkgs#awscli2`。クロージャ 422 MB）を、VM が cache.nixos.org から
`nix-store --realise` で取る（2026-10-13 に、flake.lock の nixpkgs の `awscli2-2.35.11` が cache.nixos.org にあることを確かめた）。
自前の cache は AWS CLI がないと取れないので、AWS CLI だけは公開の cache から、クロージャのハッシュで固定して取る。

### 自分を終了する方法

- **A. `--instance-initiated-shutdown-behavior terminate` で起動し、VM の中で `shutdown -h now`**: 電源を切れば終了し、EBS は
  DeleteOnTermination で消える。IAM の権限が要らない
- **B. IMDSv2 でインスタンス ID を取り、`aws ec2 terminate-instances`**: `ec2:TerminateInstances` をロールに足す（`ec2:ResourceTag/Project = xlc`
  の条件つき）。VM の中の権限が増えるのに、A と比べて止められる場面は増えない（どちらも VM の OS が動いていることが前提）。
  A は `run-instances --dry-run`（スポット one-time、`--instance-initiated-shutdown-behavior terminate`）が通ることを確かめた（2026-10-13）

安全装置として、user-data の最初に `shutdown -h +<最長の分>`（既定 840 分 = 14 時間）を置く。本番の 1 本は r995 の 16 スレッドで
約 5.5 時間（ADR-0151）で、c7a.4xlarge（16 物理コア、Zen4）での速さは確かめていないので、約 2.5 倍の余裕をとる。時間切れでも
statepoint は S3 にあるので、再開で続けられる。OS が止まって shutdown も走らない場合は、外から止めるもの（CloudWatch の
アラームの terminate など）がないと止まらない。この場面は確かめていない稀なものとし、作らない（`aws_status.sh` で見る）。

### 中断のあとの再開

- **A. 手元から起動し直す**: 起動のツールが S3 の `runs/<名前>/<ジョブ>/` を見て、`done` がなく statepoint があるジョブを
  `--restart` で起動する。中断の自動の再起動は作らない
- **B. Spot Fleet（maintain）や Auto Scaling group で自動の再起動**: 起動テンプレート、フリートの後始末、再起動のたびに
  `--restart` を判断する仕組みが VM 側に要る。中断は数本の計算で数回の見込みで、手元から 1 コマンドで起動し直せば足りる

### 保存先

S3 のバケット 1 つ（us-east-1、パブリックアクセスなし、バージョニングなし）。このアカウントは組織の SCP が `s3:CreateBucket` を拒むので、
バケットは管理者が作る（2026-10-13）:

- `store/nix/`（クロージャの binary cache）、`store/openmc-library/<内容のハッシュ>/`（核データ）: 残す。内容のハッシュで分け、
  ライブラリを作り直しても古い版を上書きしない
- `runs/<名前>/`（起動の記録 `launch.env`、入力とスクリプト `src/`、ジョブごとの結果 `<ジョブ>/`）: ライフサイクルで 30 日で失効。
  取り出して検証したら手元のツールで消す。失効は消し忘れの安全装置
- `runs/` の完了しなかったマルチパートのアップロードは 2 日で消す

EBS のボリュームを DeleteOnTermination=false で残す案は、取り出すのに VM が要り、要求（手元から定期的に取り出す）に合わない。

### VM とジョブの対応

- **A. 1 つのジョブ（入力）に 1 台**: 中断と失敗の影響が 1 本に閉じる（ADR-0151）。再開の判断もジョブごとに閉じる
- **B. 1 台で複数のジョブ（`run_remote.sh` の並列）**: 台数は減るが、1 回の中断で全部が止まり、再開は全部をそろえて起動し直す

### 実装の途中で試して捨てた方式

確認（c7a.large、examples/vera/aws-check/openmc-3A.toml）とレビューで、次を試して捨てた（2026-10-13）。

- **書かれて 60 秒たたない最新の statepoint を同期から外す**（xlc のコミット 87d1227d）: 完了の直後の同期でも最後の世代の statepoint が
  外れ、`--delete` で古いものも消えて、S3 に statepoint が 1 つもなくなった（確認の check6。`aws_fetch.sh` の SHA256SUMS が
  `run/statepoint.200.h5` で失敗）。また長い同期の途中に書き始めたものは外せない。走っている間は、同期の一覧の時点で書き終えた
  statepoint だけを含める方式にした
- **走っている間も `--delete` で同期する**: 一覧のあとに VM で statepoint が消されると、含めたはずのものを S3 から消す（レビューの指摘。
  実行中の間引きが #0106 で効いていないので、今は起きない）。`--delete` は成功のあとだけにし、コピーしてから消す
- **終了コード 128 以上を中断とみなす**（xlc のコミット 3665a524）: OOM（137）などでジョブ自身が落ちたのも中断になり、起動し直すたびに
  同じ失敗を繰り返す。user-data が SIGTERM を受けたかで分けた
- **準備（Nix、クロージャ、核データ、入力の取得）の失敗でも `done` に終了コードを置く**: 一時的な失敗でジョブが恒久的に失敗の扱いになる。
  準備では S3 のジョブに触れないので、`done` を置かずに起動し直せるようにした
- **クロージャの置き場は `nix-cache-info` があれば送らない**: sync が途中で切れても「ある」とみなす。毎回 `--size-only` で sync する
- **起動し直しでは S3 の statepoint をすべて戻す**: 走っている間は S3 の statepoint が消えないので、本番では約 25 個を EBS に戻す。
  最新の 2 個だけ戻す

## 決定と理由

OpenMC の本番は、user-data で自己完結する VM で流す。手元のツールは起動、状態の一覧、取り出しだけを行い、VM の計算の途中に手元は
要らない。

- 1 つのジョブ（入力）に 1 台。中断と失敗の影響が 1 本に閉じる（ADR-0151）
- VM は起動すると、最長時間の shutdown を予約し、Nix を入れ、AWS CLI を cache.nixos.org から、クロージャと核データと入力を S3 から取る。
  `run_remote.sh` で流し、決まった間隔（既定 5 分）ごとに作業ディレクトリ（書き終えた statepoint、summary.h5、ログ）を
  `aws s3 sync` で S3 にコピーする（`--delete` は完了のあとだけ）。完了したら、最後の世代のほかの statepoint を消し、sha256 の manifest を書いて同期し、
  電源を切る（terminate）。終わるときは起動ごとのログを S3 に置く。ジョブを始めたあとの失敗では終了コードを `done` に置く。
  準備（Nix、クロージャ、核データ、入力の取得）での失敗と、SIGTERM を受けたとき（中断）では `done` を置かず、起動し直せるようにする
- ジョブが走っている間は、最新の statepoint（書いている途中でありうる）と同期の一覧の後に書き始めたものを同期しない。写ったものが壊れていても、再開の
  `latest_readable`（examples/openmc_resume.py）が壊れたものを消して 1 つ前から続ける。再開の VM は、S3 から作業ディレクトリに
  戻し終えるまで同期しない（`--delete` で S3 の statepoint を消さないため）。走っている間は statepoint を `--delete` の対象から外すので、S3 の statepoint は完了まで増える（実行中の間引きも効いていない。#0106）
- インスタンスプロファイルのロール `xlc-openmc-vm` は、そのバケットの一覧（`store/`、`runs/`）、`store/` の読み取り、`runs/` の
  読み書きと削除だけ。EC2 の権限は持たない。バケットはツールの変数 `XLC_BUCKET`（既定 xlc-openmc-118645514640）で選ぶ
- 中断のあとは手元から同じ名前で起動し直し、`done` のないジョブを `--restart` で続ける。起動し直すときは、最初の起動で S3 に置いた
  `launch.env` と `src/` を使い、手元の作業ツリーの入力とスクリプトを読まない（途中で変えたものが混ざらない）。user-data を作るツールと AMI は
  起動し直したときのものなので、起動ごとに `runs/<名前>/launches/<時刻>.env` にツールのコミットと AMI を記録する
- 手元は `aws_fetch.sh` で取り出し、manifest で検証し、`--verify` で手元の extract が VM の結果とバイト一致することを確かめ、
  `--delete` で S3 の `runs/<名前>/` を消す
- 入ってくる通信は要らないので、キーペアを付けず、VPC の既定のセキュリティグループ（同じグループのインスタンスからの通信だけを許す）で起動し、IMDSv2 を必須にする。
  様子は S3 のログと `aws ec2 get-console-output` で見る

## 帰結

- 手元が落ちても、計算は最後まで進み、結果は S3 に残り、VM は終了する。中断と時間切れで失うのは、S3 にコピーした最後の statepoint
  からの分。走っている間は最新の statepoint をコピーしないので、最大で入力の `statepoint-interval` の約 2 倍と同期の間隔
- 費用（2026-10-13 の公開の単価。docs/environment.md）: S3 の保存は $0.023/GB・月で、`store/`（zstd のクロージャ 212 MB と核データ 631 MB、
  約 0.85 GB。2026-10-13 の `du`、xlc のコミット df72096a）は月約 $0.02。取り出しは月 100 GB まで無料。EC2 と EBS は VM が動いている間だけ
- VM が S3 にいつ何を書いたかは、起動ごとのログ（`runs/<名前>/<ジョブ>/boot-<起動の時刻>.log`）に残る。中断（SIGTERM）では最後の同期とログが置かれないことがあり（#0107）、
  残るのは決まった間隔の同期の分。AWS CLI の取得は cache.nixos.org に頼る
- 手元から ssh で進めるツール（remote_send.sh、run_remote.sh、remote_fetch.sh、aws_launch.sh など）は残す。run_remote.sh は
  VM の中でも使う

- ツールの整理（テストの書き方、重複、状態の一覧の呼び出しの数など）の残りは #0109

## 確認

2026-10-13、xlc のコミット 6360a566 のスクリプト（3 の起動し直しは、起動のツールだけ 93c36f8c。VM は S3 の `src/` の 6360a566 のスクリプト）、c7a.large のスポット（us-east-1a、
$0.0365/時間）、examples/vera/aws-check/openmc-3A.toml（粒子 1 万・200 世代・捨てる 50・間隔 20）、`--sync-minutes 1`。手元から
ssh はしていない（VM にキーペアがない）。

1. 自己終了: 起動から約 9 分 20 秒で、VM が自分で電源を切り terminated になった（`describe-instances` の StateReason が
   `Client.InstanceInitiatedShutdown`）。user-data の開始から run の開始まで約 1 分 30 秒（Nix、AWS CLI、クロージャ、核データ、入力の取得）。
   インスタンスプロファイル `xlc-openmc-vm` で S3 の `store/` の読み取りと `runs/` の書き込みが通った
2. 取り出し: `aws_fetch.sh --verify --delete` で SHA256SUMS が一致し、手元の extract の 4 ファイルが VM の結果とバイト一致し、
   S3 の `runs/check1/` が消えた
3. 再開: 別の実行で statepoint 100 が S3 に写ったあと手元から終了し、`aws_run.sh <名前>` で起動し直すと、`--restart` で
   statepoint 100 から続けて完了し、自己終了した。k は 1 の実行と 1e-15 の相対差で一致した（1.175018 ± 0.000747）。この確認で、
   起動し直しが statepoint を見落として最初から起動する誤り（`aws s3 ls | grep -q` が pipefail で失敗）を見つけて直した
4. レビューの修正（`done` の扱い、走っているジョブを飛ばす起動し直しなど）のあと、xlc のコミット 7a172cb0 で 3 をもう一度通した。
   statepoint 80 のあと手元から終了すると S3 に `done` は置かれず（SIGTERM）、VM が走っている間の `aws_run.sh <名前>` は
   「VM が走っているので飛ばす」で何も起動しなかった。終了のあとの起動し直しは statepoint 80 から再開して自己終了し、
   `aws_fetch.sh --verify --delete` が通り、k は 1 と同じだった。手元から終了したときは `boot.log` が S3 に置かれなかった
   （シャットダウンの途中で finish の同期が間に合わないか、走らない。どちらかは確かめていない）
5. 費用: 6 台（誤って起動して 1 分で終了した 1 台を含む）で合わせて約 31 分、約 $0.02。S3 は 1 セント未満。user-data が
   途中で失敗したときに `done` に終了コードを書く道は、実際には通していない
