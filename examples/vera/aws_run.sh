#!/usr/bin/env bash
# OpenMC の参照解を、1 つの入力に 1 台の EC2 のスポットで自己完結に流す（ADR-0155、docs/environment.md の「自己完結な実行」）。
# VM は user-data（aws_user_data.sh）で Nix と閉包、核データ、入力を S3 から取り、run_remote.sh で流し、S3 の
# runs/<name>/<ジョブ>/ に同期し、終わると自分を終了する。手元はこのコマンドのあと要らない。
#
# 使い方: examples/vera/aws_run.sh [--dry-run] [--type T] [--max-hours H] [--sync-minutes M] <name> [<入力>...]
#   <name>     実行の名前（S3 の runs/<name>/、タグ XlcRun、インスタンスの Name は <name>-<ジョブ>）
#   <入力>     examples/vera/p3/openmc-3A.toml など（リポジトリの直下からの相対。.. を含まない）。1 つにつき 1 台。
#              入力を渡すと新しい実行（runs/<name>/ がすでにあれば止まる）。閉包、核データ、スクリプトと入力を S3 に置く
#              入力を省くと起動し直し: S3 の runs/<name>/launch.env と src/ を使い（手元の作業ツリーは読まない）、
#              done のないジョブを起動する。statepoint があれば --restart で続ける。VM が走っているジョブは飛ばす。
#              done が 0 でないジョブ（失敗）と、最後の世代の statepoint があって done のないジョブは飛ばして示し、1 で終わる
#   --type     インスタンスタイプ（既定 c7a.4xlarge）
#   --max-hours  最長の稼働時間（既定 14。過ぎると VM は電源を切る。statepoint は S3 にあるので起動し直せば続く）
#   --sync-minutes  S3 に同期する間隔（既定 5）。起動し直しで 3 つの指定がなければ、最後の起動の記録の値を使う
#   --dry-run  S3 に置かず起動もせず、ジョブごとの起動の引数（aws の --dry-run で検証）と user-data の要約を出す
# バケットは XLC_BUCKET（既定 xlc-openmc-118645514640）、インスタンスプロファイルは xlc-openmc-vm。
set -euo pipefail
here=$(dirname "$0")

usage() { sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'; }
dry=0
type=
max_hours=
sync_minutes=
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
[[ $max_hours =~ ^([1-9][0-9]*)?$ && $sync_minutes =~ ^([1-9][0-9]*)?$ ]] ||
  { echo "--max-hours と --sync-minutes は正の整数" >&2; exit 2; }
# shellcheck source=aws_common.sh
source "$here/aws_common.sh"
# shellcheck source=aws_user_data.sh
source "$here/aws_user_data.sh"

run_prefix=s3://$XLC_BUCKET/runs/$name
# VM はバケットと同じリージョンで起動する（VM と S3 の間の転送を無料にし、user-data の AWS CLI のリージョンを合わせる）
bucket_region=$(aws s3api get-bucket-location --bucket "$XLC_BUCKET" --query LocationConstraint --output text)
[[ $bucket_region != None ]] || bucket_region=us-east-1
# AWS CLI は AWS_REGION を AWS_DEFAULT_REGION より先に読む
launch_region=${AWS_REGION:-$AWS_DEFAULT_REGION}
[[ $bucket_region == "$launch_region" ]] ||
  { echo "バケット $XLC_BUCKET は $bucket_region、起動のリージョンは $launch_region（AWS_REGION か AWS_DEFAULT_REGION を合わせる）" >&2; exit 1; }

# VM は AWS CLI を cache.nixos.org から取る（自前の cache は AWS CLI がないと取れない）。なければ VM は S3 に何も置けずに止まる
# 標準入力の TOML（入力）から、<key> = <値> の行の値を読む。値は <正規表現> の形（と行末のコメント）に限り、それ以外の書き方
# （2_500、複数行など。TOML として正しくても）は読み違えないよう止まる。読むのは起動の前の確認だけで、VM は tomllib で読む
input_value() { # <名前（表示用）> <key> <値の正規表現>
  local line
  line=$(grep -E "^$2[[:space:]]*=" || true)
  [[ $(wc -l <<< "$line") == 1 && $line =~ ^$2[[:space:]]*=[[:space:]]*($3)[[:space:]]*(#.*)?$ ]] ||
    { echo "$1: $2 の行を読めません（$line）" >&2; exit 1; }
  echo "${BASH_REMATCH[1]}"
}

require_public_awscli() {
  nix path-info --store https://cache.nixos.org "$1" > /dev/null 2>&1 ||
    { echo "$1 が cache.nixos.org にありません（VM が AWS CLI を取れない）" >&2; exit 1; }
}

if ((${#inputs[@]})); then
  # 新しい実行
  [[ -f flake.nix && -d examples/vera ]] || { echo "リポジトリの直下から走らせてください" >&2; exit 1; }
  [[ -f examples/vera/openmc-library/cross_sections.xml ]] ||
    { echo "examples/vera/openmc-library がありません（dvc checkout examples/vera/dvc.yaml）" >&2; exit 1; }
  # 一部だけ checkout されたライブラリを置かない（VM で OpenMC が核データを読むときに初めて失敗する）
  while read -r f; do
    [[ -f examples/vera/openmc-library/$f ]] || { echo "examples/vera/openmc-library/$f がありません（dvc checkout examples/vera/dvc.yaml）" >&2; exit 1; }
  done < <(sed -n 's/.* path="\([^"]*\)".*/\1/p' examples/vera/openmc-library/cross_sections.xml)
  keys=$(s3_keys "runs/$name/")
  if [[ -n $keys ]]; then
    echo "$run_prefix/ はすでにあります。起動し直すなら入力を省く。やり直すなら別の名前にする（古い実行を消すなら" \
      "取り出してから aws_fetch.sh --verify --delete $name、完了していなければ aws s3 rm --recursive $run_prefix/）" >&2
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
    # VM に置くのは examples/vera/openmc-library だけなので、入力の cross-sections はそこを指すこと
    xs=$(input_value "$input" cross-sections '"[^"\\]*"' < "$input")
    xs=${xs//\"/}
    [[ $(realpath -m "$(dirname "$input")/$xs") == "$(realpath examples/vera/openmc-library/cross_sections.xml)" ]] ||
      { echo "$input: cross-sections（$xs）が examples/vera/openmc-library/cross_sections.xml を指しません" >&2; exit 1; }
    jobs+=("$(job_name "$input")")
    live=$(live_instance "$name-${jobs[-1]}")
    [[ -z $live ]] || { echo "$name-${jobs[-1]} の VM が走っています" >&2; exit 1; }
  done
  dup=$(printf '%s\n' "${jobs[@]}" | sort | uniq -d)
  [[ -z $dup ]] || { echo "ジョブの名前（入力の basename）が重なります: $dup" >&2; exit 1; }

  # パスだけが要る（VM が cache.nixos.org から取る）ので、手元に 422 MB の閉包を取らない
  awscli=$(nix eval --raw --inputs-from . nixpkgs#awscli2.outPath)
  require_public_awscli "$awscli"
  # 1 つずつ取る（1 回の build の出力の順に頼らない。$(...) なので失敗で止まる）
  openmc=$(nix build .#openmc --no-link --print-out-paths)
  python=$(nix build .#openmc-python --no-link --print-out-paths)
  [[ $openmc =~ ^/nix/store/[^[:space:]]+$ && $python =~ ^/nix/store/[^[:space:]]+$ ]] ||
    { echo "閉包の出力が 1 つの store のパスでない: $openmc $python" >&2; exit 1; }
  # 鍵は内容のハッシュ。同じ閉包と核データは 1 度だけ置く
  closure=$(printf '%s\n' "$openmc" "$python" | sha256sum | cut -c1-16)
  library=$(cd examples/vera/openmc-library && find . -type f -print0 | sort -z | xargs -0 sha256sum | sha256sum | cut -c1-16)
  commit=$(git rev-parse HEAD)
  # 追跡していない入力（trials/ など。git status は無視したファイルを出さない）もコミットから再現できないので dirty
  if [[ -n $(git status --porcelain) ]] || ! git ls-files --error-unmatch "${inputs[@]}" > /dev/null 2>&1; then
    commit="$commit-dirty"
  fi

  src=$(mktemp -d)
  cache=$(mktemp -d)
  trap 'rm -rf "$src" "$cache"' EXIT
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
    nix copy --quiet --to "file://$cache?compression=zstd" "$openmc" "$python"
    # ファイル名は内容で決まるので大きさだけで比べる（作り直した cache は mtime が新しく、既定の比較では全部を送り直す）
    aws s3 sync --only-show-errors --size-only "$cache" "s3://$XLC_BUCKET/store/nix/$closure"
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
  # 指定のない --type、--max-hours、--sync-minutes は、最後の起動の記録の値を使う（起動し直しで黙って設定を変えない）
  last_key=$(s3_keys "runs/$name/launches/" | sort | tail -1)
  if [[ -n $last_key ]]; then
    last=$(aws s3 cp "s3://$XLC_BUCKET/$last_key" -)
    type=${type:-$(sed -n 's/^TYPE=//p' <<< "$last")}
    max_hours=${max_hours:-$(sed -n 's/^MAX_HOURS=//p' <<< "$last")}
    sync_minutes=${sync_minutes:-$(sed -n 's/^SYNC_MINUTES=//p' <<< "$last")}
    echo "最後の起動（${last_key##*/}）の設定を引き継ぐ（指定がなければ）"
  fi
  todo=()
  restart_flags=()
  stuck=0
  for input in "${RUN_INPUTS[@]}"; do
    job=$(job_name "$input")
    code=$(job_done "$name" "$job")
    if [[ $code == 0 ]]; then
      echo "$job: 完了"
      continue
    fi
    if [[ -n $code ]]; then
      # 失敗したジョブは飛ばし、ほかのジョブは起動する（影響を 1 本に閉じる）。終わりの状態は 1
      echo "$job: 失敗しています（done: $code）。$run_prefix/$job/boot-*.log を見る。続けるなら原因を直して done を消す" >&2
      stuck=1
      continue
    fi
    live=$(live_instance "$name-$job")
    if [[ -n $live ]]; then
      echo "$job: VM が走っているので飛ばす"
      continue
    fi
    latest=$(job_latest_statepoint "$name" "$job")
    # 最後の世代の statepoint があって done がないのは、完了したあと done を置く前に止まったとき。OpenMC はそこから再開
    # できない（世代が残っていない）ので、起動しない
    batches=
    if [[ -n $latest ]]; then
      src_input=$(aws s3 cp "$run_prefix/src/$input" -)
      batches=$(input_value "$run_prefix/src/$input" batches '[1-9][0-9]*' <<< "$src_input")
    fi
    if [[ -n $batches && $latest == "statepoint.$batches.h5" ]]; then
      echo "$job: 最後の世代の statepoint があり done がない（完了のあと done を置く前に止まった）。$run_prefix/$job/ の" \
        "SHA256SUMS と boot-*.log を確かめ、そろっていれば echo 0 | aws s3 cp - $run_prefix/$job/done で完了にする" >&2
      stuck=1
      continue
    fi
    if [[ -n $latest ]]; then
      echo "$job: $latest から再開"
      restart_flags+=(1)
    else
      echo "$job: 最初から"
      restart_flags+=(0)
    fi
    todo+=("$input")
  done
  ((${#todo[@]})) || { echo "起動するジョブはありません（aws_status.sh、aws_fetch.sh $name）"; exit "$stuck"; }
fi

# 既定（ADR-0155）。起動し直しでは、指定も記録もないときだけ
type=${type:-c7a.4xlarge}
max_hours=${max_hours:-14}
sync_minutes=${sync_minutes:-5}
ami=$(ubuntu_ami)
echo "AMI: $ami、タイプ: $type、最長 $max_hours 時間、同期 $sync_minutes 分ごと"

if ((!dry)); then
  # user-data を作ったこのスクリプトのチェックアウト（作業ディレクトリでなく）
  tools=$(git -C "$here" rev-parse HEAD)
  [[ -z $(git -C "$here" status --porcelain) ]] || tools="$tools-dirty"
fi

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
    # 起動ごとの記録（起動し直しは S3 の src/ を使うが、user-data を作る道具、AMI、機種はこの起動のもの）。起動できたものだけ書く
    printf 'TOOLS_COMMIT=%s\nAMI=%s\nTYPE=%s\nMAX_HOURS=%s\nSYNC_MINUTES=%s\nJOB=%s\nRESTART=%s\nINSTANCE=%s\n' "$tools" "$ami" \
      "$type" "$max_hours" "$sync_minutes" "$job" "${restart_flags[$i]}" "$id" |
      aws s3 cp --only-show-errors - "$run_prefix/launches/$(date -u +%Y%m%dT%H%M%S.%NZ)-$job.env"
  fi
done
((dry)) || echo "進み具合: examples/vera/aws_status.sh、終わったら examples/vera/aws_fetch.sh $name"
exit "${stuck:-0}"
