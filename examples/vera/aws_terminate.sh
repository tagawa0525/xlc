#!/usr/bin/env bash
# タグ Name が <name>（または <name>-N）のインスタンスを終了し、terminated を待って ~/.ssh/config の Host を消す
# （docs/environment.md の「外部の計算機（AWS）」）。EBS は DeleteOnTermination で消える。
#
# 使い方: examples/vera/aws_terminate.sh [--dry-run] <name>
set -euo pipefail
# shellcheck source=aws_common.sh
source "$(dirname "$0")/aws_common.sh"

dry=0
name=""
while (($#)); do
  case $1 in
    --help | -h) sed -n '2,6p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    --dry-run) dry=1 ;;
    -*) echo "知らない引数 $1" >&2; exit 2 ;;
    *) name=$1 ;;
  esac
  shift
done
[[ -n $name ]] || { sed -n '2,6p' "$0" | sed 's/^# \{0,1\}//'; exit 2; }

mapfile -t ids < <(aws ec2 describe-instances --filters "Name=tag:Name,Values=$name,$name-*" \
  "Name=instance-state-name,Values=pending,running,stopping,stopped" --query 'Reservations[].Instances[].InstanceId' --output text | tr '\t' '\n' | sed '/^$/d')
mapfile -t hosts < <(aws ec2 describe-instances --filters "Name=tag:Name,Values=$name,$name-*" \
  "Name=instance-state-name,Values=pending,running,stopping,stopped" \
  --query 'Reservations[].Instances[].Tags[?Key==`Name`].Value' --output text | tr '\t' '\n' | sed '/^$/d')
echo "対象: ${ids[*]:-なし}（${hosts[*]:-}）"
((dry)) && exit 0
if ((${#ids[@]})); then
  aws ec2 terminate-instances --instance-ids "${ids[@]}" > /dev/null
  aws ec2 wait instance-terminated --instance-ids "${ids[@]}"
  aws ec2 describe-instances --instance-ids "${ids[@]}" --query 'Reservations[].Instances[].[InstanceId,State.Name]' --output text
fi
ssh_config_remove "$name"
for h in "${hosts[@]}"; do ssh_config_remove "$h"; done
echo "終了しました"
