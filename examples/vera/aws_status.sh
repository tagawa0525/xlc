#!/usr/bin/env bash
# 走っている（終了前の）インスタンスと起動からの時間、スポット価格の目安を出す（docs/environment.md の「外部の計算機（AWS）」）。
# 続けて、S3 の runs/ にある自己完結な実行（aws_run.sh。ADR-0155）のジョブごとの done と最新の statepoint を出す。
#
# 使い方: examples/vera/aws_status.sh
set -euo pipefail
# shellcheck source=aws_common.sh
source "$(dirname "$0")/aws_common.sh"
[[ ${1:-} != --help && ${1:-} != -h ]] || { sed -n '2,5p' "$0" | sed 's/^# \{0,1\}//'; exit 0; }

now=$(date +%s)
found=0
while read -r id name type state ip launch az; do
  found=1
  hours=$(awk -v n="$now" -v l="$(date -d "$launch" +%s)" 'BEGIN {printf "%.2f", (n - l) / 3600}')
  price=$(aws ec2 describe-spot-price-history --instance-types "$type" --availability-zone "$az" --product-descriptions Linux/UNIX \
    --start-time "$(date -u +%FT%TZ)" --query 'SpotPriceHistory[0].SpotPrice' --output text)
  [[ $price != None ]] || { echo "$name $id $type $state $ip スポット価格を取得できません（$az）"; continue; }
  echo "$name $id $type $state $ip 起動から ${hours} 時間、スポット価格 \$${price}/時間（直近の市場価格（その AZ）。1 台あたり約 \$$(awk -v h="$hours" -v p="$price" 'BEGIN {printf "%.3f", h * p}')）"
done < <(aws ec2 describe-instances --filters "Name=instance-state-name,Values=pending,running,stopping,stopped" \
  --query 'Reservations[].Instances[].[InstanceId,Tags[?Key==`Name`]|[0].Value,InstanceType,State.Name,PublicIpAddress,LaunchTime,Placement.AvailabilityZone]' --output text)
((found)) || echo "走っているインスタンスはありません"

echo "== S3 の実行（s3://$XLC_BUCKET/runs/）"
mapfile -t runs < <(aws s3 ls "s3://$XLC_BUCKET/runs/" | sed -n 's#^ *PRE \(.*\)/$#\1#p')
((${#runs[@]})) || echo "ありません"
for run in "${runs[@]}"; do
  prefix=s3://$XLC_BUCKET/runs/$run
  env=$(aws s3 cp "$prefix/launch.env" - 2> /dev/null) || { echo "$run: launch.env がありません（置いている途中か失敗）"; continue; }
  read -r -a inputs <<< "$(sed -n 's/^INPUTS=//p' <<< "$env" | tr -d '"')"
  for input in "${inputs[@]}"; do
    job=$(basename "$input" .toml)
    code=$(aws s3 cp "$prefix/$job/done" - 2> /dev/null || true)
    last=$(aws s3 ls "$prefix/$job/run/" 2> /dev/null | sed -n 's/.* \(statepoint\.[0-9]*\.h5\)$/\1/p' | sort -t. -k2 -n | tail -1)
    echo "$run/$job: done ${code:-なし}、最新の statepoint ${last:-なし}"
  done
done
