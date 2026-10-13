# aws_*.sh が source する共通の定義（docs/environment.md の「外部の計算機（AWS）」）。単体では実行しない。
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-east-1}
SSH_CONFIG=${SSH_CONFIG:-$HOME/.ssh/config}
SG_ID=sg-07ed76965b2bc0539   # xlc-ssh

# ~/.ssh/config の「# BEGIN xlc-aws <host>」から「# END xlc-aws <host>」までを消す（なければ何もしない）
ssh_config_remove() {
  local host=$1 tmp
  [[ -f $SSH_CONFIG ]] || return 0
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
# END xlc-aws $host
EOT
}
