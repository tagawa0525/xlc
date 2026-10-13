# shellcheck shell=bash
# aws_*.sh が source する共通の定義（docs/environment.md の「外部の計算機（AWS）」）。単体では実行しない。
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-east-1}
SSH_CONFIG=${SSH_CONFIG:-$HOME/.ssh/config}
# 自己完結な実行（aws_run.sh、aws_fetch.sh。ADR-0155）の S3 のバケット
XLC_BUCKET=${XLC_BUCKET:-xlc-openmc-118645514640}
SG_ID=$(aws ec2 describe-security-groups --filters Name=group-name,Values=xlc-ssh \
  --query 'SecurityGroups[0].GroupId' --output text)
[[ $SG_ID == sg-* ]] || { echo "セキュリティグループ xlc-ssh が $AWS_DEFAULT_REGION にありません" >&2; exit 1; }

# ~/.ssh/config の「# BEGIN xlc-aws <host>」から「# END xlc-aws <host>」までを消す（なければ何もしない）
ssh_config_remove() {
  local host=$1 tmp
  [[ -f $SSH_CONFIG ]] || return 0
  ssh-keygen -R "$host" > /dev/null 2>&1 || true
  grep -qxF "# BEGIN xlc-aws $host" "$SSH_CONFIG" || return 0
  grep -qxF "# END xlc-aws $host" "$SSH_CONFIG" ||
    { echo "$SSH_CONFIG に「# END xlc-aws $host」がありません。手で直してください" >&2; return 1; }
  tmp=$(mktemp)
  awk -v b="# BEGIN xlc-aws $host" -v e="# END xlc-aws $host" \
    '$0==b {skip=1} !skip {print} $0==e {skip=0}' "$SSH_CONFIG" > "$tmp"
  cat "$tmp" > "$SSH_CONFIG"
  rm -f "$tmp"
}

# Host <host> を書く（既存の同名は置き換え）
ssh_config_set() {
  local host=$1 ip=$2
  mkdir -p "$(dirname "$SSH_CONFIG")"
  ssh_config_remove "$host"
  cat >> "$SSH_CONFIG" <<EOT
# BEGIN xlc-aws $host
Host $host
    HostName $ip
    User ubuntu
    IdentityFile ~/.ssh/xlc-aws.pem
    IdentitiesOnly yes
    StrictHostKeyChecking accept-new
    HostKeyAlias $host
# END xlc-aws $host
EOT
}
