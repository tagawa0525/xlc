# shellcheck shell=bash
# aws_*.sh が source する共通の定義（docs/environment.md の「外部の計算機（AWS）」）。単体では走らせない。
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-east-1}
# $(...) の中の失敗でも止める（既定では command substitution の中で set -e が切れ、S3 の失敗が空の結果に化ける）
shopt -s inherit_errexit
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

# 自己完結な実行 <name> の S3 の launch.env を LAUNCH_ENV に、入力を RUN_INPUTS（配列）に、定義のパスを RUN_DEFINITION に読み、
# 定義（S3 の src/ のもの）を def_parse で DEF_* に入れる。なければ止まる。--allow-no-definition を付けると、DEFINITION のない
# （定義に対応する前の）実行は RUN_DEFINITION を空にして戻る（一覧を出す aws_status.sh が、古い実行 1 つで止まらないため）
read_launch_env() {
  LAUNCH_ENV=$(aws s3 cp "s3://$XLC_BUCKET/runs/$1/launch.env" -) ||
    { echo "s3://$XLC_BUCKET/runs/$1/launch.env がありません" >&2; exit 1; }
  # shellcheck disable=SC2034 # source した道具が使う
  read -r -a RUN_INPUTS <<< "$(launch_env_value INPUTS)"
  # 実行の定義（ADR-0174）は、VM が使った S3 の src/ のもの。なければ（定義の前の実行）読み込みの互換は作らず止まる
  RUN_DEFINITION=$(launch_env_value DEFINITION)
  if [[ -z $RUN_DEFINITION && ${2:-} == --allow-no-definition ]]; then return 0; fi
  [[ -n $RUN_DEFINITION ]] ||
    { echo "s3://$XLC_BUCKET/runs/$1/launch.env に DEFINITION がありません（実行の定義に対応する前の実行は扱えない。ADR-0174）" >&2; exit 1; }
  def_parse "$(aws s3 cp "s3://$XLC_BUCKET/runs/$1/src/$RUN_DEFINITION" -)"
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
job_done() {
  local key=runs/$1/$2/done keys
  keys=$(s3_keys "$key") || exit 1
  if grep -qxF "$key" <<< "$keys"; then
    aws s3 cp "s3://$XLC_BUCKET/$key" - || { echo "s3://$XLC_BUCKET/$key を読めません" >&2; exit 1; }
  fi
}

# S3 の runs/ の下の実行の名前を 1 行ずつ。失敗では止まる
s3_runs() {
  local out
  out=$(aws s3api list-objects-v2 --bucket "$XLC_BUCKET" --prefix runs/ --delimiter / \
    --query 'CommonPrefixes[].Prefix' --output text) || { echo "s3://$XLC_BUCKET/runs/ を一覧できません" >&2; exit 1; }
  [[ $out == None ]] || tr '\t' '\n' <<< "$out" | sed 's#^runs/##; s#/$##'
}

# 標準入力の TOML（入力）から、<key> = <値> の行の値を読む。値は <正規表現> の形（と行末のコメント）に限り、それ以外の書き方
# （2_500、複数行など。TOML として正しくても）は読み違えないよう止まる。読むのは起動の前の確認だけで、VM は tomllib で読む
input_value() { # <名前（表示用）> <key> <値の正規表現>
  local line
  line=$(grep -E "^$2[[:space:]]*=" || true)
  [[ $(wc -l <<< "$line") == 1 && $line =~ ^$2[[:space:]]*=[[:space:]]*($3)[[:space:]]*(#.*)?$ ]] ||
    { echo "$1: $2 の行を読めません（$line）" >&2; exit 1; }
  echo "${BASH_REMATCH[1]}"
}

# S3 の runs/<name>/<job>/run/ の進み具合のうち最新のもの（定義の progress による。read_launch_env のあとに呼ぶ）。なければ空
#   statepoint: statepoint.<世代>.h5（ADR-0168）   steps: 完了の印 complete のある区間 step-<k>（ADR-0174）
job_latest_progress() {
  local keys
  case $DEF_PROGRESS in
    statepoint)
      keys=$(s3_keys "runs/$1/$2/run/statepoint.") || exit 1
      sed -n 's#.*/\(statepoint\.[0-9]*\.h5\)$#\1#p' <<< "$keys" | sort -t. -k2 -n | tail -1
      ;;
    steps)
      keys=$(s3_keys "runs/$1/$2/run/step-") || exit 1
      sed -n 's#.*/\(step-[0-9]*\)/complete$#\1#p' <<< "$keys" | sort -t- -k2 -n | tail -1
      ;;
  esac
}

# 最新の進み具合 <latest>（job_latest_progress の値）が、再開できないジョブの終わりを指すなら 0。<input> は S3 の src/ の入力。
# statepoint は入力の batches と最後の世代が一致するとき（世代が残っていないので再開できない）。steps は、最後の区間まで
# 終わっていても再開できる（区間を流さずに extract を流して done を置く）ので、常に 1
job_progress_finished() { # <name> <job> <latest> <input>
  local batches
  case $DEF_PROGRESS in
    statepoint)
      # 条件の中で呼ばれると set -e が効かないので、失敗は明示して止める
      batches=$(aws s3 cp "s3://$XLC_BUCKET/runs/$1/src/$4" - | input_value "s3://$XLC_BUCKET/runs/$1/src/$4" batches '[1-9][0-9]*') || exit 1
      [[ $3 == "statepoint.$batches.h5" ]]
      ;;
    steps) return 1 ;;
  esac
}

# 実行の定義（examples/<問題>/aws.toml。ADR-0174）の本文 <text> を読み、次の変数に入れる。キーの綴りや形が違えば止まる。
#   DEF_TOOLS（配列）  VM の src/ に置く道具（リポジトリの直下からの相対。入口 DEF_ENTRY と定義自身は含めない）
#   DEF_ENTRY          VM が流す入口のスクリプト
#   DEF_PROGRESS       進み具合の種類（statepoint か steps）
#   DEF_LIBRARY        OpenMC のライブラリのディレクトリ（cross_sections.xml と全 path= のファイルがある）
#   DEF_DATA（配列）   そのほかの入力ファイル（省略できる）
#   DEF_JOBS           受け付けるジョブ名（入力の basename から .toml を除いたもの）の正規表現
# TOML の部分集合だけを読む: 「キー = 値」の行、値は "文字列"（英数字と ._/-。絶対パスと .. は不可）、その配列（複数行にできる）、
# jobs だけは '正規表現'。# 以降はコメントなので、値に # は使えない。戻りの値でなく exit で止まるので、条件の中で呼ばない
def_parse() {
  local text=$1 line buf='' key value inner part seen=' ' re rc
  # shellcheck disable=SC2034 # source した道具が使う
  DEF_TOOLS=() DEF_ENTRY='' DEF_PROGRESS='' DEF_LIBRARY='' DEF_DATA=() DEF_JOBS=''
  while IFS= read -r line; do
    line=${line%%#*}
    [[ -n $buf || ! $line =~ ^[[:space:]]*$ ]] || continue
    buf+="${buf:+ }$line"
    # 配列は ] で閉じるまで次の行を足す
    if [[ $buf =~ ^[[:space:]]*[a-z_]+[[:space:]]*=[[:space:]]*\[ && ! $buf =~ \][[:space:]]*$ ]]; then continue; fi
    [[ $buf =~ ^[[:space:]]*([a-z_]+)[[:space:]]*=[[:space:]]*(.*[^[:space:]])[[:space:]]*$ ]] ||
      { echo "実行の定義: 読めない行です: $buf" >&2; exit 1; }
    key=${BASH_REMATCH[1]} value=${BASH_REMATCH[2]} buf=''
    [[ $seen != *" $key "* ]] || { echo "実行の定義: $key が重複しています" >&2; exit 1; }
    seen+="$key "
    case $key in
      tools | data)
        [[ $value =~ ^\[[[:space:]]*(.*[^,[:space:]])?[[:space:]]*,?[[:space:]]*\]$ ]] ||
          { echo "実行の定義: $key は文字列の配列です: $value" >&2; exit 1; }
        inner=${BASH_REMATCH[1]}
        local -a items=()
        if [[ -n $inner ]]; then
          while IFS= read -r part; do
            def_path_string "$key" "$part"
            items+=("$DEF_STR")
          done < <(tr ',' '\n' <<< "$inner")
        fi
        # shellcheck disable=SC2034 # source した道具が使う
        if [[ $key == tools ]]; then DEF_TOOLS=("${items[@]}"); else DEF_DATA=("${items[@]}"); fi
        ;;
      entry) def_path_string entry "$value"; DEF_ENTRY=$DEF_STR ;;
      library) def_path_string library "$value"; DEF_LIBRARY=$DEF_STR ;;
      progress)
        def_path_string progress "$value"
        [[ $DEF_STR == statepoint || $DEF_STR == steps ]] ||
          { echo "実行の定義: progress は \"statepoint\" か \"steps\": $value" >&2; exit 1; }
        DEF_PROGRESS=$DEF_STR
        ;;
      jobs)
        [[ $value =~ ^\'([^\']+)\'$ ]] || { echo "実行の定義: jobs は '正規表現': $value" >&2; exit 1; }
        re=${BASH_REMATCH[1]} rc=0
        # [[ =~ ]] は一致しなければ 1、正規表現が不正なら 2
        # shellcheck disable=SC2319
        [[ x =~ $re ]] || rc=$?
        ((rc != 2)) || { echo "実行の定義: jobs が正規表現として読めません: $re" >&2; exit 1; }
        DEF_JOBS=$re
        ;;
      *) echo "実行の定義: 知らないキー $key" >&2; exit 1 ;;
    esac
  done <<< "$text"
  [[ -z $buf ]] || { echo "実行の定義: 配列が閉じていません: $buf" >&2; exit 1; }
  ((${#DEF_TOOLS[@]})) && [[ -n $DEF_ENTRY && -n $DEF_PROGRESS && -n $DEF_LIBRARY && -n $DEF_JOBS ]] ||
    { echo "実行の定義: tools（1 つ以上）、entry、progress、library、jobs が要ります" >&2; exit 1; }
}

# 定義の値 <値> を "文字列" として読み、リポジトリの直下からの相対のパスなら DEF_STR に入れる。<key> はエラーの表示用
def_path_string() {
  local value=${2#"${2%%[![:space:]]*}"}
  value=${value%"${value##*[![:space:]]}"}
  [[ $value =~ ^\"([A-Za-z0-9._/-]+)\"$ ]] || { echo "実行の定義: $1 の値は \"文字列\"（英数字と ._/-）: $2" >&2; exit 1; }
  DEF_STR=${BASH_REMATCH[1]}
  [[ $DEF_STR != /* && /$DEF_STR/ != */../* ]] ||
    { echo "実行の定義: $1 はリポジトリの直下からの相対のパスにしてください（絶対パスと .. は不可）: $DEF_STR" >&2; exit 1; }
}

# 入力 <path> から上のディレクトリへたどって、aws.toml のある最初のものを出す（リポジトリの直下から走らせる）。なければ止まる
def_find() {
  local d
  d=$(dirname "$1")
  while [[ $d != . && $d != / ]]; do
    if [[ -f $d/aws.toml ]]; then echo "$d/aws.toml"; return 0; fi
    d=$(dirname "$d")
  done
  echo "$1: 実行の定義（入力の親のディレクトリにある aws.toml。ADR-0174）がありません" >&2
  exit 1
}

# AZ <az> の default サブネットの ID を SUBNET_ID に入れる（現在のリージョンの AZ でなければ、なければ止まる）。
# AZ を固定して起動するための --subnet-id に渡す
az_subnet_id() {
  [[ $1 =~ ^[a-z]{2}(-[a-z]+)+-[0-9][a-z]$ ]] || { echo "AZ の名前でない: $1" >&2; exit 2; }
  SUBNET_ID=$(aws ec2 describe-subnets --filters Name=default-for-az,Values=true "Name=availability-zone,Values=$1" \
    --query 'Subnets[0].SubnetId' --output text)
  [[ $SUBNET_ID == subnet-* ]] ||
    { echo "AZ $1 の default サブネットが ${AWS_REGION:-$AWS_DEFAULT_REGION} にありません（AZ のリージョンは AWS_REGION、なければ AWS_DEFAULT_REGION と同じ）" >&2; exit 1; }
}

# タグ Name が <host> の終了前（終了中を含む）のインスタンスの ID（なければ空）
live_instance() {
  aws ec2 describe-instances --filters "Name=tag:Name,Values=$1" \
    "Name=instance-state-name,Values=pending,running,shutting-down,stopping,stopped" --query 'Reservations[].Instances[].InstanceId' --output text
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

# 機種の設定（aws_types.toml）の [speed] の "タイプ" = 速度 を「タイプ 速度」で 1 行ずつ出す。書き方が違う行は止まる
# （速度は c7a.4xlarge を 1 とした見積もりの値。ADR-0171）
auto_types_read() {
  local file=$1 line type speed in_speed=0 in_vcpus=0 n=0 seen=' '
  [[ -f $file ]] || { echo "機種の設定 $file がありません" >&2; exit 1; }
  while IFS= read -r line; do
    [[ $line =~ ^[[:space:]]*(#.*)?$ ]] && continue
    if [[ $line =~ ^\[speed\][[:space:]]*(#.*)?$ ]]; then in_speed=1; in_vcpus=0; continue; fi
    if [[ $line =~ ^\[vcpus\][[:space:]]*(#.*)?$ ]]; then in_speed=0; in_vcpus=1; continue; fi
    # [vcpus] は auto_vcpus_read が読む。ここでは形だけ確かめる
    if ((in_vcpus)) && [[ $line =~ ^\"[a-z0-9-]+\.[a-z0-9-]+\"[[:space:]]*=[[:space:]]*[1-9][0-9]*[[:space:]]*(#.*)?$ ]]; then continue; fi
    if ((in_speed)) && [[ $line =~ ^\"([a-z0-9-]+\.[a-z0-9-]+)\"[[:space:]]*=[[:space:]]*([0-9]+(\.[0-9]+)?)[[:space:]]*(#.*)?$ ]]; then
      type=${BASH_REMATCH[1]} speed=${BASH_REMATCH[2]}
      # 0 や 0.0 は速度にならない
      if [[ $speed =~ [1-9] ]]; then
        [[ $seen != *" $type "* ]] || { echo "$file: $type が重複しています" >&2; exit 1; }
        seen="$seen$type "
        echo "$type $speed"
        n=$((n + 1))
        continue
      fi
    fi
    echo "$file: 読めない行です（[speed] の下に \"タイプ\" = 速度、[vcpus] の下に \"タイプ\" = vCPU 数）: $line" >&2
    exit 1
  done < "$file"
  ((n)) || { echo "$file: 候補の機種がありません" >&2; exit 1; }
}

# 機種の設定 <file> の候補について、直近 1 時間のスポット価格を AZ ごとに引き（aws_spot_prices.sh と同じく AZ ごとの最新）、
# 1 本の費用 = 価格 ÷ 速度 の昇順に「費用 タイプ AZ 価格 速度」を 1 行ずつ出す。<az> があればその AZ だけ。価格のない
# （その AZ に提供のない）組み合わせは出ない
auto_candidates() {
  local file=$1 only_az=${2:-} speeds types prices start
  speeds=$(auto_types_read "$file")
  types=$(awk '{print $1}' <<< "$speeds" | tr '\n' ' ')
  start=$(date -u -d '1 hours ago' +%FT%TZ)
  # shellcheck disable=SC2086 # types は空白区切りの複数引数
  prices=$(aws ec2 describe-spot-price-history --instance-types $types --product-descriptions Linux/UNIX --start-time "$start" \
    --query 'SpotPriceHistory[].[InstanceType,AvailabilityZone,Timestamp,SpotPrice]' --output text) ||
    { echo "スポット価格を引けません" >&2; exit 1; }
  [[ -n $prices ]] || return 0
  sort -k1,1 -k2,2 -k3,3r <<< "$prices" | awk -v speeds="$speeds" -v only="$only_az" '
    BEGIN { n = split(speeds, l, "\n"); for (i = 1; i <= n; i++) { split(l[i], f, " "); sp[f[1]] = f[2] } }
    !seen[$1 " " $2]++ && (only == "" || $2 == only) && sp[$1] > 0 { printf "%.6f %s %s %s %s\n", $4 / sp[$1], $1, $2, $4, sp[$1] }' |
    LC_ALL=C sort -k1,1g -k2,2 -k3,3
}

# 機種の設定 <file> の [vcpus]（機種の vCPU 数。クォータに収まる機種の選別に使う）から「タイプ vCPU」を 1 行ずつ出す。
# 候補の機種（[speed]）に vCPU がなければ止まる
auto_vcpus_read() {
  local file=$1 line in_vcpus=0 type vcpus
  vcpus=$(while IFS= read -r line; do
    if [[ $line =~ ^\[vcpus\][[:space:]]*(#.*)?$ ]]; then in_vcpus=1; continue; fi
    if [[ $line =~ ^\[ ]]; then in_vcpus=0; continue; fi
    if ((in_vcpus)) && [[ $line =~ ^\"([a-z0-9-]+\.[a-z0-9-]+)\"[[:space:]]*=[[:space:]]*([1-9][0-9]*)[[:space:]]*(#.*)?$ ]]; then
      echo "${BASH_REMATCH[1]} ${BASH_REMATCH[2]}"
    fi
  done < "$file")
  while read -r type _; do
    grep -q "^$type " <<< "$vcpus" || { echo "$file: $type の vCPU 数が [vcpus] にありません" >&2; exit 1; }
  done < <(auto_types_read "$file")
  echo "$vcpus"
}

# スポットの vCPU のクォータ（L-34B43A08）から、この account で今走っている（pending/running）スポットのインスタンスの
# vCPU の合計を引いた残りを出す。タグ XlcRun に限らず全部を数える。取れなければ 1 を返す（呼び出し側が警告する）
spot_vcpu_remaining() {
  local quota used
  quota=$(aws service-quotas get-service-quota --service-code ec2 --quota-code L-34B43A08 --query Quota.Value --output text) || return 1
  [[ $quota =~ ^[0-9]+(\.[0-9]+)?$ ]] || return 1
  used=$(aws ec2 describe-instances --filters Name=instance-lifecycle,Values=spot Name=instance-state-name,Values=pending,running \
    --query 'Reservations[].Instances[].[CpuOptions.CoreCount,CpuOptions.ThreadsPerCore]' --output text |
    awk '{ s += $1 * $2 } END { print s + 0 }') || return 1
  echo "$(( ${quota%.*} - used )) ${quota%.*} $used"
}
