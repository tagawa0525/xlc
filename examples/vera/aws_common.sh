# shellcheck shell=bash
# aws_*.sh が source する共通の定義（docs/environment.md の「外部の計算機（AWS）」）。単体では実行しない。
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-east-1}
SSH_CONFIG=${SSH_CONFIG:-$HOME/.ssh/config}
# 自己完結な実行（aws_run.sh、aws_fetch.sh。ADR-0155）の S3 のバケット
XLC_BUCKET=${XLC_BUCKET:-xlc-openmc-118645514640}

# SSH で入るツール（aws_launch.sh、aws_allow_ip.sh）が使うセキュリティグループ xlc-ssh の ID を SG_ID に入れる
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
  # shellcheck disable=SC2034 # source したツールが使う
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

# バケットの <接頭辞> で始まるキーを 1 行ずつ。ないときは何も出さず、S3 に届かないなどの失敗では止まる
# （aws s3 ls は「ない」と失敗のどちらも 1 で終わり、区別できない）
s3_keys() {
  local out
  out=$(aws s3api list-objects-v2 --bucket "$XLC_BUCKET" --prefix "$1" --query 'Contents[].Key' --output text) ||
    { echo "s3://$XLC_BUCKET/$1 を一覧できません" >&2; exit 1; }
  [[ $out == None ]] || tr '\t' '\n' <<< "$out"
}

# S3 の runs/<name>/<job>/done の中身（run と extract の終了コード）。なければ空。読めなければ止まる
# （呼ぶ側は code=$(job_done ...) の代入で受け、失敗を set -e で止める）
job_done() {
  local key=runs/$1/$2/done keys
  keys=$(s3_keys "$key")
  if grep -qxF "$key" <<< "$keys"; then
    aws s3 cp "s3://$XLC_BUCKET/$key" - || { echo "s3://$XLC_BUCKET/$key を読めません" >&2; exit 1; }
  fi
}

# S3 の runs/<name>/<job>/run/ の statepoint のうち最新のファイル名。なければ空
job_latest_statepoint() {
  local keys
  keys=$(s3_keys "runs/$1/$2/run/statepoint.")
  sed -n 's#.*/\(statepoint\.[0-9]*\.h5\)$#\1#p' <<< "$keys" | sort -t. -k2 -n | tail -1
}

# ジョブの名前から、run と extract を行うスクリプト（run_remote.sh の case と同じ対応。VM の上では run_remote.sh が選ぶ）。
# 知らない名前なら空
job_script() {
  case $1 in
    openmc-3?) echo examples/vera/openmc_model_3d.py ;;
    openmc-5?-2D) echo examples/vera/openmc_core_2d.py ;;
  esac
}

# Ubuntu 24.04 の最新の AMI（このアカウントは SCP で SSM の公開パラメータが引けないので describe-images）
ubuntu_ami() {
  aws ec2 describe-images --owners 099720109477 \
    --filters "Name=name,Values=ubuntu/images/hvm-ssd-gp3/ubuntu-noble-24.04-amd64-server-*" \
    --query 'sort_by(Images,&CreationDate)[-1].ImageId' --output text
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
