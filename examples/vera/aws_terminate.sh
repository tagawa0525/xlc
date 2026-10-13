#!/usr/bin/env bash
# タグ Name が <name>（または <name>-N）のインスタンスと、aws_run.sh の実行 <name> のインスタンス（タグ XlcRun が <name>）を終了し、terminated を待って ~/.ssh/config の Host を消す
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
    --help | -h) sed -n '2,5p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    --dry-run) dry=1 ;;
    -*) echo "知らない引数 $1" >&2; exit 2 ;;
    *) name=$1 ;;
  esac
  shift
done
[[ -n $name ]] || { sed -n '2,5p' "$0" | sed 's/^# \{0,1\}//'; exit 2; }

# <name> または <name>-<数字> に厳密に一致するものと、XlcRun が <name> のものだけ（xlc で xlc-ssh-test を巻き込まない）
ids=()
hosts=()
# 失敗で止まるよう、一覧を先に取る（< <(...) の失敗は set -e で止まらない）
listing=$(aws ec2 describe-instances --filters "Name=tag:Name,Values=$name,$name-*" \
  "Name=instance-state-name,Values=pending,running,stopping,stopped" \
  --query 'Reservations[].Instances[].[InstanceId,Tags[?Key==`Name`]|[0].Value,Tags[?Key==`XlcRun`]|[0].Value]' --output text)
while read -r id host run; do
  [[ $host == "$name" || $host =~ ^$name-[0-9]+$ || $run == "$name" ]] || continue
  ids+=("$id")
  hosts+=("$host")
done <<< "$listing"
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
