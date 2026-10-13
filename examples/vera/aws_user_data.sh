# shellcheck shell=bash
# aws_run.sh が source する。1 つのジョブを自己完結に流す VM の user-data（cloud-init のシェルスクリプト）を標準出力に書く
# （ADR-0155）。AWS を呼ばないので、単体の検査（test_aws_user_data.sh）でも source できる。
#
# aws_user_data <region> <bucket> <run> <job> <入力> <再開 0|1> <最長の分> <同期の秒> <ライブラリの鍵> <クロージャの鍵> <openmc> <python> <awscli> <間隔>
#   <job>      results/ の下の名前（入力の basename から .toml を除いたもの。run_remote.sh と同じ）
#   <入力>     作業ディレクトリ（S3 の runs/<run>/src/ をコピーしたもの）からの相対
#   <ライブラリの鍵>、<クロージャの鍵>  S3 の store/openmc-library/<鍵>/、store/nix/<鍵>/
#   <間隔>     入力の statepoint-interval を上書きする世代数（正の整数。空なら上書きしない。ADR-0168）
aws_user_data() {
  (($# == 14)) || { echo "aws_user_data: 引数は 14 個（$# 個）" >&2; return 2; }
  local region=$1 bucket=$2 run=$3 job=$4 input=$5 restart=$6 minutes=$7 seconds=$8 library=$9
  local closure=${10} openmc=${11} python=${12} awscli=${13} interval=${14} v
  [[ $restart == 0 || $restart == 1 ]] || { echo "aws_user_data: 再開は 0 か 1（$restart）" >&2; return 2; }
  [[ $minutes =~ ^[1-9][0-9]*$ && $seconds =~ ^[1-9][0-9]*$ ]] ||
    { echo "aws_user_data: 最長の分と同期の秒は正の整数（$minutes、$seconds）" >&2; return 2; }
  [[ $interval =~ ^([1-9][0-9]*)?$ ]] || { echo "aws_user_data: 間隔は正の整数か空（$interval）" >&2; return 2; }
  for v in "$openmc" "$python" "$awscli"; do
    [[ $v =~ ^/nix/store/[a-z0-9]{32}-[^/[:space:]]+$ ]] || { echo "aws_user_data: store のパスではない: $v" >&2; return 2; }
  done
  for v in "$region" "$bucket" "$run" "$job" "$input" "$library" "$closure"; do
    [[ $v =~ ^[A-Za-z0-9._/-]+$ ]] || { echo "aws_user_data: 使えない文字: $v" >&2; return 2; }
  done
  cat <<USERDATA
#!/bin/bash
# aws_run.sh が作った user-data（ADR-0155）。cloud-init が root で 1 度だけ実行する
set -euo pipefail
# 安全装置: 何が起きても最長 $minutes 分で電源を切る（InstanceInitiatedShutdownBehavior=terminate なので終了する）
shutdown -h +$minutes
exec > >(tee -a /var/log/xlc-job.log) 2>&1
TEE=\$!
echo "xlc: start \$(date -u +%FT%TZ)"

export HOME=/root AWS_DEFAULT_REGION=$region
BUCKET=$bucket
PREFIX=s3://\$BUCKET/runs/$run/$job
AWS=$awscli/bin/aws
NIX=/nix/var/nix/profiles/default/bin
WORK=/var/lib/xlc
OUT=results/$job

# 起動ごとのログの名前（起動し直しで前の VM のログを上書きしない）
BOOTLOG=boot-\$(date -u +%Y%m%dT%H%M%SZ).log
# 同期: done とログは finish が最後に置く（S3 に done があれば、ほかのファイルはすべて置いてある）。
# READY（作業ディレクトリが S3 の内容を含むようになった）までは同期しない。再開で S3 から戻す前に --delete でコピーすると、
# S3 の statepoint を消すため
READY=0
# <ディレクトリ> の statepoint のファイル名を世代の順に
statepoints() {
  ls "\$1" 2> /dev/null | grep -E '^statepoint\.[0-9]+\.h5\$' | sort -t. -k2 -n || true
}
sync_out() {
  ((READY)) || return 0
  local filters=() delete=(--delete) f all mtime
  # ジョブが走っている間は、一覧の後に書き始めたものと、書いている途中でありうるものをコピーしない。コピーするのは、一覧の時点で最新で
  # ないもの（書き終えた）と、最新でも mtime が 120 秒より古いもの（OpenMC は数秒で書き終え、書く間は mtime が進む。ADR-0168）。
  # ジョブが成功した（done が 0）ときだけ、時刻によらずすべてコピーする（失敗や中断の直後は、最新の statepoint が書いている途中で
  # 止まったかもしれず、mtime が新しいので除かれる）。そのあいだは --delete もしない（一覧のあとに statepoint が
  # 消されると、含めたはずのものを S3 から消すため）
  if ((INTERRUPTED)) || [[ \$(cat "\$WORK/\$OUT/done" 2> /dev/null) != 0 ]]; then
    delete=()
    filters=(--exclude 'run/statepoint.*')
    all=\$(statepoints "\$WORK/\$OUT/run")
    for f in \$(head -n -1 <<< "\$all"); do filters+=(--include "run/\$f"); done
    f=\$(tail -n 1 <<< "\$all")
    if [[ -n \$f ]] && mtime=\$(stat -c %Y "\$WORK/\$OUT/run/\$f" 2> /dev/null) && ((\$(date +%s) - mtime > 120)); then
      filters+=(--include "run/\$f")
    fi
  fi
  # --delete のときは、先に消さずにコピーしてから消す（途中で失敗しても、S3 から statepoint が消えて最後のものもない、にならない）
  "\$AWS" s3 sync --only-show-errors --exclude done --exclude 'boot-*.log' "\${filters[@]}" "\$WORK/\$OUT" "\$PREFIX"
  if ((\${#delete[@]})); then
    "\$AWS" s3 sync --only-show-errors --delete --exclude done --exclude 'boot-*.log' "\$WORK/\$OUT" "\$PREFIX"
  fi
}
# 終わるときに、結果とログを S3 に置いてから電源を切る。done は、ジョブが終わったとき（run_remote.sh が書く）と、
# READY のあとでこの user-data が失敗したとき（run_remote.sh の done を上書きして終了コードを書く。完了の後の整理の失敗も
# 含む）に置く。SIGTERM（スポットの回収、最長時間の shutdown、手元からの終了）と READY の前の失敗（準備。S3 のジョブには
# 触れていない）では書かず、中断として aws_run.sh で起動し直せるようにする（原因は boot-*.log に残る）。
# SIGTERM のときの同期とログは、シャットダウンの途中で間に合わないことがある（確かめた範囲では置かれなかった）。中断で残るのは
# 決まった間隔の同期の分。同期に 3 回失敗したら、S3 の結果がそろわないので done を置かない
finish() {
  local code=\$? synced=0 try
  # ここで止まると電源が切れないので、途中の失敗では止めない。シャットダウンで tee が先に消えても、
  # 出力の SIGPIPE で止まらない。finish の途中の SIGTERM（最長時間、手元からの終了）で exit 143 に飛ばない
  set +e
  trap '' PIPE TERM
  echo "xlc: exit \$code \$(date -u +%FT%TZ)"
  # SIGTERM を受けた（中断）ときは done を置かない。シャットダウンでジョブが殺されて書いた done も消し、同期が書いている途中の
  # statepoint を除くようにする。ジョブ自身がシグナルで落ちた（OOM の 137 など）のは失敗として done に残す
  if ((INTERRUPTED)); then
    rm -f "\$WORK/\$OUT/done"
  fi
  if ((READY && code != 0 && !INTERRUPTED)); then
    mkdir -p "\$WORK/\$OUT"
    echo "\$code" > "\$WORK/\$OUT/done"
  fi
  if [[ -x \$AWS ]]; then
    for try in 1 2 3; do
      sync_out && { synced=1; break; }
      sleep 10
    done
    # tee が書き終えてからログを置く（最後の行が欠けないように）
    exec > /dev/null 2>&1
    wait "\$TEE" 2> /dev/null || sleep 1
    "\$AWS" s3 cp --only-show-errors /var/log/xlc-job.log "\$PREFIX/\$BOOTLOG" || true
    if ((synced && !INTERRUPTED)) && [[ -f \$WORK/\$OUT/done ]]; then
      "\$AWS" s3 cp --only-show-errors "\$WORK/\$OUT/done" "\$PREFIX/done" || true
    fi
  fi
  shutdown -h now
}
trap finish EXIT
INTERRUPTED=0
trap 'INTERRUPTED=1; exit 143' TERM

# installer の版を固定する（起動し直しの VM でも同じ Nix になる）。v3.23.1 は 2026-10-13 の確認で入った版
curl -fsSL https://install.determinate.systems/nix/tag/v3.23.1 | sh -s -- install --no-confirm
# AWS CLI は cache.nixos.org から（自前の cache を取るのに要る）
"\$NIX/nix-store" --realise $awscli > /dev/null

mkdir -p "\$WORK"
cd "\$WORK"
"\$AWS" s3 sync --only-show-errors "s3://\$BUCKET/store/nix/$closure" cache
"\$NIX/nix" copy --no-check-sigs --from "file://\$WORK/cache" $openmc $python
rm -rf cache
"\$AWS" s3 sync --only-show-errors "s3://\$BUCKET/store/openmc-library/$library" examples/vera/openmc-library
"\$AWS" s3 sync --only-show-errors "s3://\$BUCKET/runs/$run/src" .
printf 'OPENMC_BIN=%s\nPYTHON=%s\n' $openmc/bin $python/bin/python3 > tools.env

# S3 のジョブ（前の起動のログなど）を戻してから同期を始める。statepoint は最新の 2 個だけ戻す（再開に要るのは最新の読めるもの。
# 走っている間の同期は S3 の statepoint を消さないので、全部を戻すと EBS を使う）
"\$AWS" s3 sync --only-show-errors --exclude 'boot-*.log' --exclude done --exclude 'run/statepoint.*' "\$PREFIX" "\$OUT"
keys=\$("\$AWS" s3api list-objects-v2 --bucket "\$BUCKET" --prefix "runs/$run/$job/run/statepoint." --query 'Contents[].Key' --output text)
if [[ \$keys != None ]]; then
  for key in \$(tr '\t' '\n' <<< "\$keys" | sort -t. -k2 -n | tail -2); do
    "\$AWS" s3 cp --only-show-errors "s3://\$BUCKET/\$key" "\$OUT/run/\${key##*/}"
  done
fi
READY=1
USERDATA
  local flags=
  ((restart)) && flags=" --restart"
  [[ -z $interval ]] || flags="$flags --statepoint-interval $interval"
  echo "bash examples/vera/run_remote.sh$flags $input"

  cat <<USERDATA
# run_remote.sh は done を 1 回の echo で書くので、空のうちに読まないよう数字になるまで待つ
# 10 秒ごとに done とジョブを見て、$seconds 秒ごとに同期する。run_remote.sh が切り離したジョブ（bash -c job ... results/<job>
# と python ... results/<job>/run）が done を書かずに消えたら失敗とする（最長の時間まで待たない）
waited=0
until [[ -f \$OUT/done && \$(cat "\$OUT/done") =~ ^[0-9]+\$ ]]; do
  # 2 度見る（done を書いた直後に消えた場合と、起動の直後でまだ exec していない場合）
  if ! pgrep -f -- " \$OUT(/run)?( |\$)" > /dev/null; then
    sleep 5
    if [[ ! -f \$OUT/done ]] && ! pgrep -f -- " \$OUT(/run)?( |\$)" > /dev/null; then
      echo "xlc: the job ended without done \$(date -u +%FT%TZ)"
      exit 1
    fi
  fi
  sleep 10
  waited=\$((waited + 10))
  if ((waited >= $seconds)); then
    waited=0
    sync_out || echo "xlc: sync failed \$(date -u +%FT%TZ)"
  fi
done
code=\$(cat "\$OUT/done")
if [[ \$code == 0 ]]; then
  # 完了: 最後の世代のほかの statepoint を消し、manifest を書く
  last=\$(statepoints "\$OUT/run" | tail -1)
  find "\$OUT/run" -name 'statepoint.*.h5' ! -name "\$last" -delete
  (cd "\$OUT" && find . -type f ! -name SHA256SUMS ! -name done -printf '%P\0' | sort -z | xargs -0 sha256sum > SHA256SUMS)
fi
exit "\$code"
USERDATA
}
