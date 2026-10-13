#!/usr/bin/env bash
# OpenMC の参照解を、1 つの入力に 1 台の EC2 のスポットで自己完結に流す（ADR-0155、docs/environment.md の「自己完結な実行」）。
# VM は user-data（aws_user_data.sh）で Nix とクロージャ、核データ、入力を S3 から取り、run_remote.sh で流し、S3 の
# runs/<name>/<ジョブ>/ に同期し、終わると自分を終了する。手元はこのコマンドのあと要らない。
#
# 使い方: examples/vera/aws_run.sh [--dry-run] [--type T] [--max-hours H] [--sync-minutes M] <name> [<入力>...]
#   <name>     実行の名前（S3 の runs/<name>/、タグ XlcRun、インスタンスの Name は <name>-<ジョブ>）
#   <入力>     examples/vera/p3/openmc-3A.toml など（リポジトリの直下からの相対。.. を含まない）。1 つにつき 1 台。
#              入力を渡すと新しい実行（runs/<name>/ がすでにあれば止まる）。クロージャ、核データ、スクリプトと入力を S3 に置く
#              入力を省くと起動し直し: S3 の runs/<name>/launch.env と src/ を使い（手元の作業ツリーは読まない）、
#              done のないジョブを起動する。statepoint があれば --restart で続ける。VM が走っているジョブは飛ばす。
#              done が 0 でないジョブ（失敗）があれば、何も起動せずに止まる
#   --type     インスタンスタイプ（既定 c7a.4xlarge）
#   --max-hours  最長の稼働時間（既定 14。過ぎると VM は電源を切る。statepoint は S3 にあるので起動し直せば続く）
#   --sync-minutes  S3 に同期する間隔（既定 5）
#   --dry-run  S3 に置かず起動もせず、ジョブごとの起動の引数（aws の --dry-run で検証）と user-data の要約を出す
# バケットは XLC_BUCKET（既定 xlc-openmc-118645514640）、インスタンスプロファイルは xlc-openmc-vm。
set -euo pipefail
here=$(dirname "$0")

usage() { sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'; }
dry=0
type=c7a.4xlarge
max_hours=14
sync_minutes=5
name=""
inputs=()
while (($#)); do
  case $1 in
    --help | -h) usage; exit 0 ;;
    --dry-run) dry=1 ;;
    --type) type=$2; shift ;;
    --max-hours) max_hours=$2; shift ;;
    --sync-minutes) sync_minutes=$2; shift ;;
    -*) echo "知らない引数 $1" >&2; exit 2 ;;
    *) if [[ -z $name ]]; then name=$1; else inputs+=("$1"); fi ;;
  esac
  shift
done
[[ -n $name ]] || { usage; exit 2; }
[[ $name =~ ^[a-z0-9][a-z0-9-]*$ ]] || { echo "name は英小文字、数字、- だけ: $name" >&2; exit 2; }
[[ $max_hours =~ ^[1-9][0-9]*$ && $sync_minutes =~ ^[1-9][0-9]*$ ]] ||
  { echo "--max-hours と --sync-minutes は正の整数" >&2; exit 2; }
# shellcheck source=aws_common.sh
source "$here/aws_common.sh"
# shellcheck source=aws_user_data.sh
source "$here/aws_user_data.sh"

run_prefix=s3://$XLC_BUCKET/runs/$name
# VM はバケットと同じリージョンで起動する（VM と S3 の間の転送を無料にし、user-data の AWS CLI のリージョンを合わせる）
bucket_region=$(aws s3api get-bucket-location --bucket "$XLC_BUCKET" --query LocationConstraint --output text)
[[ $bucket_region != None ]] || bucket_region=us-east-1
[[ $bucket_region == "$AWS_DEFAULT_REGION" ]] ||
  { echo "バケット $XLC_BUCKET は $bucket_region、起動のリージョンは $AWS_DEFAULT_REGION（AWS_DEFAULT_REGION を合わせる）" >&2; exit 1; }

# タグ Name が <host> の終了前のインスタンスの ID（なければ空）
live_instance() {
  aws ec2 describe-instances --filters "Name=tag:Name,Values=$1" \
    "Name=instance-state-name,Values=pending,running,shutting-down,stopping,stopped" --query 'Reservations[].Instances[].InstanceId' --output text
}

# VM は AWS CLI を cache.nixos.org から取る（自前の cache は AWS CLI がないと取れない）。なければ VM は S3 に何も置けずに止まる
require_public_awscli() {
  nix path-info --store https://cache.nixos.org "$1" > /dev/null 2>&1 ||
    { echo "$1 が cache.nixos.org にありません（VM が AWS CLI を取れない）" >&2; exit 1; }
}

if ((${#inputs[@]})); then
  # 新しい実行
  [[ -f flake.nix && -d examples/vera ]] || { echo "リポジトリの直下から実行してください" >&2; exit 1; }
  [[ -f examples/vera/openmc-library/cross_sections.xml ]] ||
    { echo "examples/vera/openmc-library がありません（dvc checkout examples/vera/dvc.yaml）" >&2; exit 1; }
  keys=$(s3_keys "runs/$name/")
  if [[ -n $keys ]]; then
    echo "$run_prefix/ はすでにあります。起動し直すなら入力を省く。やり直すなら aws_fetch.sh $name --delete か別の名前" >&2
    exit 1
  fi
  jobs=()
  for input in "${inputs[@]}"; do
    # VM の作業ディレクトリに同じ相対パスで置くので、絶対パスと .. は受け付けない
    [[ $input != /* && /$input/ != */../* ]] || { echo "$input: リポジトリの直下からの相対パスにしてください" >&2; exit 1; }
    [[ $input =~ ^[A-Za-z0-9._/-]+$ ]] || { echo "$input: 使えない文字があります（英数字と ._/-）" >&2; exit 1; }
    [[ -f $input ]] || { echo "$input がありません" >&2; exit 1; }
    # run_remote.sh が受け付ける名前（VM に上げてから断られないように、置く前に確かめる）
    [[ -n $(job_script "$(job_name "$input")") ]] ||
      { echo "$input: Problem 3（openmc-3?）か 5-2D（openmc-5?-2D）の入力ではありません" >&2; exit 1; }
    jobs+=("$(job_name "$input")")
    [[ -z $(live_instance "$name-${jobs[-1]}") ]] || { echo "$name-${jobs[-1]} の VM が走っています" >&2; exit 1; }
  done
  dup=$(printf '%s\n' "${jobs[@]}" | sort | uniq -d)
  [[ -z $dup ]] || { echo "ジョブの名前（入力の basename）が重なります: $dup" >&2; exit 1; }

  mapfile -t outs < <(nix build .#openmc .#openmc-python --no-link --print-out-paths)
  awscli=$(nix build --inputs-from . nixpkgs#awscli2 --no-link --print-out-paths)
  require_public_awscli "$awscli"
  openmc=${outs[0]}
  python=${outs[1]}
  # 鍵は内容のハッシュ。同じクロージャと核データは 1 度だけ置く
  closure=$(printf '%s\n' "$openmc" "$python" | sha256sum | cut -c1-16)
  library=$(cd examples/vera/openmc-library && find . -type f -print0 | sort -z | xargs -0 sha256sum | sha256sum | cut -c1-16)
  commit=$(git rev-parse HEAD)
  # 追跡していない入力（trials/ など。git status は無視したファイルを出さない）もコミットから再現できないので dirty
  if [[ -n $(git status --porcelain) ]] || ! git ls-files --error-unmatch "${inputs[@]}" > /dev/null 2>&1; then
    commit="$commit-dirty"
  fi

  src=$(mktemp -d)
  trap 'rm -rf "$src"' EXIT
  files=(
    examples/openmc_common.py
    examples/openmc_run.py
    examples/openmc_resume.py
    examples/vera/openmc_tallies.py
    examples/vera/openmc_model.py
    examples/vera/openmc_model_3d.py
    examples/vera/openmc_core_2d.py
    examples/vera/run_remote.sh
    "${inputs[@]}"
  )
  for f in "${files[@]}"; do [[ -f $f ]] || { echo "$f がありません" >&2; exit 1; }; done
  cp --parents "${files[@]}" "$src/"
  cat > "$src/launch.env" <<ENV
COMMIT=$commit
OPENMC=$openmc
PYTHON=$python
AWSCLI=$awscli
CLOSURE=$closure
LIBRARY=$library
INPUTS="${inputs[*]}"
ENV
  echo "== 実行 $name（$commit）"
  cat "$src/launch.env"

  if ((!dry)); then
    # 毎回作って sync する（手元で約 2 秒）。途中で切れた置き場を「ある」と見なさない
    cache=$(mktemp -d)
    nix copy --quiet --to "file://$cache?compression=zstd" "$openmc" "$python"
    # ファイル名は内容で決まるので大きさだけで比べる（作り直した cache は mtime が新しく、既定の比較では全部を送り直す）
    aws s3 sync --only-show-errors --size-only "$cache" "s3://$XLC_BUCKET/store/nix/$closure"
    rm -rf "$cache"
    aws s3 sync --only-show-errors --size-only examples/vera/openmc-library "s3://$XLC_BUCKET/store/openmc-library/$library"
    aws s3 sync --only-show-errors --exclude launch.env "$src" "$run_prefix/src"
    # launch.env は最後に置く（あれば src/ と store/ はそろっている）
    aws s3 cp --only-show-errors "$src/launch.env" "$run_prefix/launch.env"
  fi
  todo=("${inputs[@]}")
  restart_flags=()
  for _ in "${todo[@]}"; do restart_flags+=(0); done
else
  # 起動し直し
  read_launch_env "$name"
  openmc=$(launch_env_value OPENMC) python=$(launch_env_value PYTHON) awscli=$(launch_env_value AWSCLI)
  closure=$(launch_env_value CLOSURE) library=$(launch_env_value LIBRARY)
  require_public_awscli "$awscli"
  echo "== 起動し直し $name（$(launch_env_value COMMIT)）"
  todo=()
  restart_flags=()
  for input in "${RUN_INPUTS[@]}"; do
    job=$(job_name "$input")
    code=$(job_done "$name" "$job")
    if [[ $code == 0 ]]; then
      echo "$job: 完了"
      continue
    fi
    [[ -z $code ]] || { echo "$job: 失敗しています（done: $code）。$run_prefix/$job/boot.log を見てください" >&2; exit 1; }
    if [[ -n $(live_instance "$name-$job") ]]; then
      echo "$job: VM が走っているので飛ばす"
      continue
    fi
    latest=$(job_latest_statepoint "$name" "$job")
    if [[ -n $latest ]]; then
      echo "$job: $latest から再開"
      restart_flags+=(1)
    else
      echo "$job: 最初から"
      restart_flags+=(0)
    fi
    todo+=("$input")
  done
  ((${#todo[@]})) || { echo "起動するジョブはありません（aws_status.sh、aws_fetch.sh $name）"; exit 0; }
fi

ami=$(ubuntu_ami)
echo "AMI: $ami、タイプ: $type、最長 $max_hours 時間、同期 $sync_minutes 分ごと"

for i in "${!todo[@]}"; do
  input=${todo[$i]}
  job=$(job_name "$input")
  host=$name-$job
  user_data=$(aws_user_data "$bucket_region" "$XLC_BUCKET" "$name" "$job" "$input" "${restart_flags[$i]}" $((max_hours * 60)) \
    $((sync_minutes * 60)) "$library" "$closure" "$openmc" "$python" "$awscli")
  args=(--image-id "$ami" --instance-type "$type" --count 1
    --iam-instance-profile Name=xlc-openmc-vm
    --instance-initiated-shutdown-behavior terminate
    --metadata-options 'HttpTokens=required,HttpEndpoint=enabled'
    --instance-market-options 'MarketType=spot,SpotOptions={SpotInstanceType=one-time}'
    --block-device-mappings 'DeviceName=/dev/sda1,Ebs={VolumeSize=40,VolumeType=gp3,DeleteOnTermination=true}'
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$host},{Key=Project,Value=xlc},{Key=XlcRun,Value=$name}]"
    --user-data "$user_data"
    --query 'Instances[0].InstanceId' --output text)
  if ((dry)); then
    out=$(aws ec2 run-instances "${args[@]}" --dry-run 2>&1 || true)
    [[ $out == *DryRunOperation* ]] || { echo "$out" >&2; exit 1; }
    echo "$host: 起動の引数は有効（再開 ${restart_flags[$i]}、user-data $(wc -c <<< "$user_data") バイト）"
  else
    id=$(aws ec2 run-instances "${args[@]}")
    echo "$host: $id を起動（再開 ${restart_flags[$i]}）"
  fi
done
((dry)) || echo "進み具合: examples/vera/aws_status.sh、終わったら examples/vera/aws_fetch.sh $name"
