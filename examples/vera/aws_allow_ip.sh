#!/usr/bin/env bash
# セキュリティグループ xlc-ssh の 22/tcp を、手元の今のグローバル IP の /32 だけにする（docs/environment.md の「外部の計算機（AWS）」）。
# 古い規則は消す。すでに一致なら何もしない。
#
# 使い方: examples/vera/aws_allow_ip.sh [--dry-run] [--ip <IPv4>]
#   --dry-run  変更せず、何をするかだけ出す
#   --ip       checkip.amazonaws.com の代わりにこの IP を使う（IP の変更の確認用）
set -euo pipefail
# shellcheck source=aws_common.sh
source "$(dirname "$0")/aws_common.sh"

dry=0
ip=""
while (($#)); do
  case $1 in
    --help | -h) sed -n '2,8p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    --dry-run) dry=1 ;;
    --ip) ip=$2; shift ;;
    *) echo "知らない引数 $1" >&2; exit 2 ;;
  esac
  shift
done

[[ -n $ip ]] || ip=$(curl -fsS https://checkip.amazonaws.com | tr -d '[:space:]')
[[ $ip =~ ^[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+$ ]] || { echo "IPv4 として読めません: $ip" >&2; exit 1; }
want="$ip/32"

mapfile -t have < <(aws ec2 describe-security-groups --group-ids "$SG_ID" \
  --query 'SecurityGroups[0].IpPermissions[?FromPort==`22` && ToPort==`22` && IpProtocol==`tcp`].IpRanges[].CidrIp' --output text | tr '\t' '\n' | sed '/^$/d')

if [[ ${#have[@]} -eq 1 && ${have[0]} == "$want" ]]; then
  echo "xlc-ssh: 22/tcp は $want のまま（変更なし）"
  exit 0
fi
echo "xlc-ssh: 22/tcp ${have[*]:-（規則なし）} -> $want"
((dry)) && exit 0
for old in "${have[@]}"; do
  [[ $old == "$want" ]] ||
    aws ec2 revoke-security-group-ingress --group-id "$SG_ID" --protocol tcp --port 22 --cidr "$old" > /dev/null
done
if ! printf '%s\n' "${have[@]}" | grep -qx "$want"; then
  aws ec2 authorize-security-group-ingress --group-id "$SG_ID" --protocol tcp --port 22 --cidr "$want" > /dev/null
fi
echo "更新しました"
