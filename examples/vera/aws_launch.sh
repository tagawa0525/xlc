#!/usr/bin/env bash
# AWS EC2 のスポットで Ubuntu 24.04 の VM を起動し、~/.ssh/config に Host を書く（docs/environment.md の「外部の計算機（AWS）」）。
# 起動の前に aws_allow_ip.sh を呼ぶ（新しい VM に入れなければ意味がないため）。
#
# 使い方: examples/vera/aws_launch.sh [--dry-run] <name> [--count N] [--type c7a.4xlarge]
#   <name>   タグ Name と ssh の Host。--count が 2 以上なら <name>-1 ... <name>-N
#   --count  台数（既定 1）
#   --type   インスタンスタイプ（既定 c7a.4xlarge。16 物理コア、SMT なし、32 GB）
#   --dry-run  AMI と起動の引数を検証（aws の --dry-run）するだけで、起動しない
# スポット（one-time）、EBS は gp3 40 GB、キーペア xlc-aws、セキュリティグループ xlc-ssh。
# 起動後に lscpu の物理コア数と free -g を出す。
set -euo pipefail
# shellcheck source=aws_common.sh
source "$(dirname "$0")/aws_common.sh"

dry=0
name=""
count=1
type=c7a.4xlarge
while (($#)); do
  case $1 in
    --help | -h) sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    --dry-run) dry=1 ;;
    --count) count=$2; shift ;;
    --type) type=$2; shift ;;
    -*) echo "知らない引数 $1" >&2; exit 2 ;;
    *) [[ -z $name ]] || { echo "name は 1 つ" >&2; exit 2; }; name=$1 ;;
  esac
  shift
done
[[ -n $name ]] || { sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; exit 2; }
[[ $count =~ ^[1-9][0-9]*$ ]] || { echo "--count は正の整数" >&2; exit 2; }

ami=$(aws ec2 describe-images --owners 099720109477 \
  --filters "Name=name,Values=ubuntu/images/hvm-ssd-gp3/ubuntu-noble-24.04-amd64-server-*" \
  --query 'sort_by(Images,&CreationDate)[-1].ImageId' --output text)
echo "AMI: $ami、タイプ: $type、台数: $count"

names=()
for ((i = 1; i <= count; i++)); do
  if ((count == 1)); then names+=("$name"); else names+=("$name-$i"); fi
done

# 同名が走っていたら止める（ssh の Host を取り違えるため）
for n in "${names[@]}"; do
  live=$(aws ec2 describe-instances --filters "Name=tag:Name,Values=$n" \
    "Name=instance-state-name,Values=pending,running,stopping,stopped" --query 'Reservations[].Instances[].InstanceId' --output text)
  [[ -z $live ]] || { echo "$n はすでにあります: $live（aws_terminate.sh）" >&2; exit 1; }
done

run_args=(--image-id "$ami" --instance-type "$type" --count 1 --key-name xlc-aws --security-group-ids "$SG_ID"
  --instance-market-options 'MarketType=spot,SpotOptions={SpotInstanceType=one-time}'
  --block-device-mappings 'DeviceName=/dev/sda1,Ebs={VolumeSize=40,VolumeType=gp3,DeleteOnTermination=true}'
  --query 'Instances[0].InstanceId' --output text)
if ((dry)); then
  out=$(aws ec2 run-instances "${run_args[@]}" --dry-run 2>&1 || true)
  [[ $out == *DryRunOperation* ]] || { echo "$out" >&2; exit 1; }
  echo "dry-run: 起動の引数は有効"
  exit 0
fi

"$(dirname "$0")/aws_allow_ip.sh"
ids=()
for n in "${names[@]}"; do
  ids+=("$(aws ec2 run-instances "${run_args[@]}" --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$n}]")")
  echo "$n: ${ids[-1]} を起動"
done
aws ec2 wait instance-running --instance-ids "${ids[@]}"
for i in "${!names[@]}"; do
  ip=$(aws ec2 describe-instances --instance-ids "${ids[$i]}" --query 'Reservations[0].Instances[0].PublicIpAddress' --output text)
  ssh_config_set "${names[$i]}" "$ip"
  echo "${names[$i]}: $ip"
  for _ in $(seq 60); do
    ssh -o ConnectTimeout=5 -o BatchMode=yes "${names[$i]}" true 2> /dev/null && break
    sleep 3
  done
  ssh "${names[$i]}" "echo '-- ${names[$i]}'; lscpu | grep -E '^(Model name|Thread\(s\) per core|Core\(s\) per socket|Socket\(s\)|CPU\(s\))'; free -g | head -2"
done
