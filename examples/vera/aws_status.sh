#!/usr/bin/env bash
# 走っている（終了前の）インスタンスと起動からの時間、スポット価格の目安を出す（docs/environment.md の「外部の計算機（AWS）」）。
#
# 使い方: examples/vera/aws_status.sh
set -euo pipefail
# shellcheck source=aws_common.sh
source "$(dirname "$0")/aws_common.sh"
[[ ${1:-} != --help && ${1:-} != -h ]] || { sed -n '2,4p' "$0" | sed 's/^# \{0,1\}//'; exit 0; }

now=$(date +%s)
found=0
while read -r id name type state ip launch; do
  found=1
  hours=$(awk -v n="$now" -v l="$(date -d "$launch" +%s)" 'BEGIN {printf "%.2f", (n - l) / 3600}')
  price=$(aws ec2 describe-spot-price-history --instance-types "$type" --product-descriptions Linux/UNIX \
    --start-time "$(date -u +%FT%TZ)" --query 'SpotPriceHistory[0].SpotPrice' --output text)
  echo "$name $id $type $state $ip 起動から ${hours} 時間、スポット価格 \$${price}/時間（直近の市場価格（AZ のうち先頭の 1 つ）。1 台あたり約 \$$(awk -v h="$hours" -v p="$price" 'BEGIN {printf "%.3f", h * p}')）"
done < <(aws ec2 describe-instances --filters "Name=instance-state-name,Values=pending,running,stopping,stopped" \
  --query 'Reservations[].Instances[].[InstanceId,Tags[?Key==`Name`]|[0].Value,InstanceType,State.Name,PublicIpAddress,LaunchTime]' --output text)
((found)) || echo "走っているインスタンスはありません"
