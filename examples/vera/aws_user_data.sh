# shellcheck shell=bash
# aws_run.sh が source する。1 つのジョブを自己完結に流す VM の user-data（cloud-init のシェルスクリプト）を標準出力に書く
# （ADR-0155）。AWS を呼ばないので、単体の検査（test_aws_user_data.sh）でも source できる。
#
# aws_user_data <bucket> <run> <job> <入力> <再開 0|1> <最長の分> <同期の秒> <ライブラリの鍵> <閉包の鍵> <openmc> <python> <awscli>
#   <job>      results/ の下の名前（入力の basename から .toml を除いたもの。run_remote.sh と同じ）
#   <入力>     作業ディレクトリ（S3 の runs/<run>/src/ を写したもの）からの相対
#   <ライブラリの鍵>、<閉包の鍵>  S3 の store/openmc-library/<鍵>/、store/nix/<鍵>/
aws_user_data() {
  (($# == 12)) || { echo "aws_user_data: 引数は 12 個（$# 個）" >&2; return 2; }
  local bucket=$1 run=$2 job=$3 input=$4 restart=$5 minutes=$6 seconds=$7 library=$8 closure=$9
  local openmc=${10} python=${11} awscli=${12} v
  [[ $restart == 0 || $restart == 1 ]] || { echo "aws_user_data: 再開は 0 か 1（$restart）" >&2; return 2; }
  [[ $minutes =~ ^[1-9][0-9]*$ && $seconds =~ ^[1-9][0-9]*$ ]] ||
    { echo "aws_user_data: 最長の分と同期の秒は正の整数（$minutes、$seconds）" >&2; return 2; }
  for v in "$openmc" "$python" "$awscli"; do
    [[ $v =~ ^/nix/store/[a-z0-9]{32}-[^/[:space:]]+$ ]] || { echo "aws_user_data: store のパスではない: $v" >&2; return 2; }
  done
  for v in "$bucket" "$run" "$job" "$input" "$library" "$closure"; do
    [[ $v =~ ^[A-Za-z0-9._/-]+$ ]] || { echo "aws_user_data: 使えない文字: $v" >&2; return 2; }
  done
  cat <<USERDATA
#!/bin/bash
# aws_run.sh が作った user-data（ADR-0155）。cloud-init が root で 1 度だけ走らせる
set -euo pipefail
# 安全装置: 何が起きても最長 $minutes 分で電源を切る（InstanceInitiatedShutdownBehavior=terminate なので終了する）
shutdown -h +$minutes
exec > >(tee -a /var/log/xlc-job.log) 2>&1
echo "xlc: start \$(date -u +%FT%TZ)"

export HOME=/root AWS_DEFAULT_REGION=us-east-1
BUCKET=$bucket
PREFIX=s3://\$BUCKET/runs/$run/$job
AWS=$awscli/bin/aws
NIX=/nix/var/nix/profiles/default/bin
WORK=/var/lib/xlc
OUT=results/$job

# 同期: done と boot.log は finish が最後に置く（S3 に done があれば、ほかのファイルはすべて置いてある）
sync_out() {
  [[ -d \$WORK/\$OUT ]] || return 0
  "\$AWS" s3 sync --only-show-errors --delete --exclude done --exclude boot.log "\$WORK/\$OUT" "\$PREFIX"
}
# どこで終わっても、結果とログを S3 に置いてから電源を切る。done は、ジョブが終わったとき（run_remote.sh が書く）と、
# この user-data が失敗したとき（run_remote.sh の done を上書きして終了コードを書く。完了の後の整理の失敗も含む）に置く。SIGTERM（スポットの回収、最長時間の shutdown、手元からの終了）
# では書かず、中断として aws_run.sh で起動し直せるようにする。同期に失敗したら、S3 の結果がそろわないので done を置かない
finish() {
  local code=\$? synced=0
  echo "xlc: exit \$code \$(date -u +%FT%TZ)"
  if ((code != 0 && code != 143)); then
    mkdir -p "\$WORK/\$OUT"
    echo "\$code" > "\$WORK/\$OUT/done"
  fi
  if [[ -x \$AWS ]]; then
    sync_out && synced=1
    "\$AWS" s3 cp --only-show-errors /var/log/xlc-job.log "\$PREFIX/boot.log" || true
    if ((synced && code != 143)); then
      "\$AWS" s3 cp --only-show-errors "\$WORK/\$OUT/done" "\$PREFIX/done" || true
    fi
  fi
  shutdown -h now
}
trap finish EXIT
trap 'exit 143' TERM

curl -fsSL https://install.determinate.systems/nix | sh -s -- install --no-confirm
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

USERDATA
  if ((restart)); then
    cat <<USERDATA
# 再開: S3 の statepoint とログを戻す
"\$AWS" s3 sync --only-show-errors --exclude boot.log --exclude done "\$PREFIX" "\$OUT"
bash examples/vera/run_remote.sh --restart $input
USERDATA
  else
    cat <<USERDATA
bash examples/vera/run_remote.sh $input
USERDATA
  fi
  cat <<USERDATA
# run_remote.sh は done を 1 回の echo で書くので、空のうちに読まないよう数字になるまで待つ
until [[ -f \$OUT/done && \$(cat "\$OUT/done") =~ ^[0-9]+\$ ]]; do
  sleep $seconds
  sync_out || echo "xlc: sync failed \$(date -u +%FT%TZ)"
done
code=\$(cat "\$OUT/done")
if [[ \$code == 0 ]]; then
  # 完了: 最後の世代のほかの statepoint を消し、manifest を書く
  last=\$(ls "\$OUT/run" | grep -E '^statepoint\.[0-9]+\.h5\$' | sort -t. -k2 -n | tail -1)
  find "\$OUT/run" -name 'statepoint.*.h5' ! -name "\$last" -delete
  (cd "\$OUT" && find . -type f ! -name SHA256SUMS ! -name done -printf '%P\0' | sort -z | xargs -0 sha256sum > SHA256SUMS)
fi
exit "\$code"
USERDATA
}
