# shellcheck shell=bash
# aws_*.sh が source する共通の定義（docs/environment.md の「外部の計算機（AWS）」）。単体では走らせない。
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-east-1}
SSH_CONFIG=${SSH_CONFIG:-$HOME/.ssh/config}
# 自己完結な実行（aws_run.sh、aws_fetch.sh。ADR-0155）の S3 のバケット
XLC_BUCKET=${XLC_BUCKET:-xlc-openmc-118645514640}

# SSH で入る道具（aws_launch.sh、aws_allow_ip.sh）が使うセキュリティグループ xlc-ssh の ID を SG_ID に入れる
# （自己完結な実行は使わないので、source しただけでは引かない）
ssh_sg_id() {
  SG_ID=$(aws ec2 describe-security-groups --filters Name=group-name,Values=xlc-ssh \
    --query 'SecurityGroups[0].GroupId' --output text)
  [[ $SG_ID == sg-* ]] || { echo "セキュリティグループ xlc-ssh が $AWS_DEFAULT_REGION にありません" >&2; exit 1; }
}

# 自己完結な実行 <name> の S3 の launch.env を LAUNCH_ENV に、入力を RUN_INPUTS（配列）に読む。なければ止まる
read_launch_env() {
  LAUNCH_ENV=$(aws s3 cp "s3://$XLC_BUCKET/runs/$1/launch.env" -) ||
    { echo "s3://$XLC_BUCKET/runs/$1/launch.env がありません" >&2; exit 1; }
  # shellcheck disable=SC2034 # source した道具が使う
  read -r -a RUN_INPUTS <<< "$(launch_env_value INPUTS)"
}

# LAUNCH_ENV の <key> の値（launch.env は aws_run.sh が KEY=値 の行で書く。INPUTS だけ引用符つき）
launch_env_value() {
  sed -n "s/^$1=//p" <<< "$LAUNCH_ENV" | tr -d '"'
}

# 入力のパスからジョブの名前（results/ と S3 の runs/<name>/ の下の名前。run_remote.sh と同じ）
job_name() {
  basename "$1" .toml
}

# S3 の runs/<name>/<job>/done の中身（run と extract の終了コード）。なければ空
job_done() {
  aws s3 cp "s3://$XLC_BUCKET/runs/$1/$2/done" - 2> /dev/null || true
}

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
