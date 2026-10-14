#!/usr/bin/env bash
# OpenMC の参照解を、1 つの入力に 1 台の EC2 のスポットで自己完結に流す（ADR-0155、docs/environment.md の「自己完結な実行」）。
# VM は user-data（aws_user_data.sh）で Nix と閉包、核データ、入力を S3 から取り、run_remote.sh で流し、S3 の
# runs/<name>/<ジョブ>/ に同期し、終わると自分を終了する。手元はこのコマンドのあと要らない。
#
# 使い方: examples/vera/aws_run.sh [--dry-run] [--type T] [--max-hours H] [--sync-minutes M] [--statepoint-interval N] [--az AZ] [--on-demand | --spot] [--job J]... [--retry-failed] <name> [<入力>...]
#   <name>     実行の名前（S3 の runs/<name>/、タグ XlcRun、インスタンスの Name は <name>-<ジョブ>）
#   <入力>     examples/vera/p3/openmc-3A.toml など（リポジトリの直下からの相対。.. を含まない）。1 つにつき 1 台。
#              入力を渡すと新しい実行（runs/<name>/ がすでにあれば止まる）。閉包、核データ、スクリプトと入力を S3 に置く
#              入力を省くと起動し直し: S3 の runs/<name>/launch.env と入力を使い、道具（run_remote.sh、openmc_*.py。新しい実行が
#              src/ に置くものと同じ）は手元の作業ツリーの版を src/ に置き直してから、done のないジョブを起動する（入力と
#              launch.env は書き換えない。道具の版は launches/ の TOOLS_COMMIT）。statepoint があれば --restart で続ける。
#              VM が走っているジョブは飛ばす。done が 0 でないジョブ（失敗）と、最後の世代の statepoint があって done のない
#              ジョブは飛ばして示し、1 で終わる
#   --type     インスタンスタイプ（既定 c7a.4xlarge）。auto なら、aws_types.toml の候補の機種と直近 1 時間のスポット価格から、
#              1 本の費用（価格 ÷ 速度）の安い (機種, AZ) の順に起動を試す（ADR-0171。空き、クォータ、価格で失敗したら次へ。
#              全部失敗なら 1 で終わり、起動し直しで再試行できる）。--az と併用すると AZ を固定して機種だけ選ぶ。--on-demand とは
#              併用できない。スポットの vCPU の残り（クォータ − 走行中）に収まらない機種は外す（欠けたジョブが複数なら残りを数で割る）。起動し直しでは、最後の起動が auto の選択なら --type 省略も auto になり、auto のとき記録の AZ を引き継がない
#   --max-hours  最長の稼働時間（既定 14。過ぎると VM は電源を切る。statepoint は S3 にあるので起動し直せば続く）
#   --sync-minutes  S3 に同期する間隔（既定 5）。起動し直しで 3 つの指定がなければ、最後の起動の記録の値を使う
#   --statepoint-interval  起動し直しで、その起動の OpenMC の走行だけ、入力の statepoint-interval を N 世代に上書きする（中断が続くとき、
#              失う分を縮める。S3 の src/ の入力は書き換えない。物理の結果には影響しない）。毎回明示する（launches/ の記録から引き継がず、
#              指定がなければ入力の値）。記録には STATEPOINT_INTERVAL= を残す。新しい実行では指定できない（入力に書く）
#   --az       AZ を固定する（例 us-east-1b）。その AZ の default サブネットで起動する。無指定は AWS が選ぶ。起動し直しでは、
#              指定がなければ最後の起動の記録の値（--job があれば、そのジョブの記録。固定を外すのは --az any）
#   --on-demand  スポットでなくオンデマンドで起動する（中断が続くときの逃げ道。価格は約 2.5 倍）。起動し直しでは、
#              指定がなければ最後の起動の記録の値。スポットに戻すのは --spot
#   --job      起動し直しで、このジョブだけを対象にする（例 --job openmc-3A。複数回指定できる。無指定は全部）。名前は launch.env の
#              INPUTS のジョブ名で、なければエラー。新しい実行では指定できない。done のあるジョブと VM が走っているジョブは今までどおり飛ばす。
#              目的: 監視のスクリプトが、中断されたジョブごとに別の AZ やオンデマンドを選んで起動し直すため
#   --retry-failed  起動し直しで、--job で指定したジョブの done が 0 でないとき、その done を S3 から消してから起動し直す
#              （--job の指定が要る。指定のないジョブは今までどおり飛ばす）。失敗の原因を S3 の runs/<name>/<ジョブ>/log と
#              boot-*.log で確かめ、直ったと分かってから使う（原因が残っていれば同じ失敗を繰り返す）
#   --dry-run  S3 に置かず起動もせず、ジョブごとの起動の引数（aws の --dry-run で検証）と user-data の要約を出す。--type auto では
#              候補の表（費用順）を出し、最初の候補の引数を検証する
# バケットは XLC_BUCKET（既定 xlc-openmc-118645514640）、インスタンスプロファイルは xlc-openmc-vm。
set -euo pipefail
here=$(dirname "$0")

usage() { sed -n '2,36p' "$0" | sed 's/^# \{0,1\}//'; }
dry=0
type=
max_hours=
sync_minutes=
interval=
az=
market=
retry=0
name=""
inputs=()
jobs=()
while (($#)); do
  case $1 in
    --help | -h) usage; exit 0 ;;
    --az) az=${2:?--az に値が要ります}; shift ;;
    --on-demand) [[ $market != spot ]] || { echo "--on-demand と --spot は同時に指定できない" >&2; exit 2; }; market=on-demand ;;
    --spot) [[ $market != on-demand ]] || { echo "--on-demand と --spot は同時に指定できない" >&2; exit 2; }; market=spot ;;
    --job) jobs+=("${2:?--job に値が要ります}"); shift ;;
    --dry-run) dry=1 ;;
    --retry-failed) retry=1 ;;
    --type) type=$2; shift ;;
    --max-hours) max_hours=$2; shift ;;
    --sync-minutes) sync_minutes=$2; shift ;;
    --statepoint-interval) interval=${2:?--statepoint-interval に値が要ります}; shift ;;
    -*) echo "知らない引数 $1" >&2; exit 2 ;;
    *) if [[ -z $name ]]; then name=$1; else inputs+=("$1"); fi ;;
  esac
  shift
done
[[ -n $name ]] || { usage; exit 2; }
[[ $name =~ ^[a-z0-9][a-z0-9-]*$ ]] || { echo "name は英小文字、数字、- だけ: $name" >&2; exit 2; }
((${#inputs[@]} == 0 || ${#jobs[@]} == 0)) || { echo "--job は起動し直し（入力を省く形）でだけ指定できる" >&2; exit 2; }
[[ $max_hours =~ ^([1-9][0-9]*)?$ && $sync_minutes =~ ^([1-9][0-9]*)?$ && $interval =~ ^([1-9][0-9]*)?$ ]] ||
  { echo "--max-hours、--sync-minutes、--statepoint-interval は正の整数" >&2; exit 2; }
((!retry || (${#jobs[@]} > 0 && ${#inputs[@]} == 0))) ||
  { echo "--retry-failed は起動し直し（入力を省く形）で、--job と一緒にだけ指定できる" >&2; exit 2; }
((${#inputs[@]} == 0)) || [[ -z $interval ]] ||
  { echo "--statepoint-interval は起動し直し（入力を省く形）でだけ指定できる（新しい実行は入力の statepoint-interval）" >&2; exit 2; }
[[ $type != auto || $market != on-demand ]] || { echo "--type auto はスポットの価格で並べるので --on-demand と併用できない" >&2; exit 2; }
# shellcheck source=aws_common.sh
source "$here/aws_common.sh"
# shellcheck source=aws_user_data.sh
source "$here/aws_user_data.sh"
# 置く前に確かめる（誤った AZ で S3 に半端な実行を残さない。起動し直しの値は後で確かめる）
[[ -z $az || $az == any ]] || az_subnet_id "$az"

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

clear_done=()
# VM が実行する道具（src/ に置く。新しい実行も起動し直しも同じ集合）
tool_files=(
  examples/openmc_common.py
  examples/openmc_input.py
  examples/openmc_run.py
  examples/openmc_resume.py
  examples/vera/openmc_tallies.py
  examples/vera/openmc_model.py
  examples/vera/openmc_model_3d.py
  examples/vera/openmc_core_2d.py
  examples/vera/run_remote.sh
)

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
  files=("${tool_files[@]}" "${inputs[@]}")
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
  all_jobs=()
  for input in "${RUN_INPUTS[@]}"; do all_jobs+=("$(job_name "$input")"); done
  for want in "${jobs[@]}"; do
    grep -Fxq -- "$want" <<< "$(printf '%s\n' "${all_jobs[@]}")" || { echo "--job $want というジョブはない（${all_jobs[*]}）" >&2; exit 2; }
  done
  require_public_awscli "$awscli"
  echo "== 起動し直し $name（$(launch_env_value COMMIT)）"
  # 指定のない --type、--max-hours、--sync-minutes は、最後の起動の記録の値を使う（起動し直しで黙って設定を変えない）
  # --job があれば、そのジョブの記録だけから引き継ぐ（別のジョブの AZ や市場を黙って持ち込まない）。複数なら、その中で最後のもの
  last_keys=$(s3_keys "runs/$name/launches/" | sort)
  if ((${#jobs[@]})); then
    pattern=$(printf '%s\n' "${jobs[@]}" | sed 's/[][\.*^$/]/\\&/g' | paste -sd'|')
    last_keys=$(grep -E -- "-($pattern)\.env\$" <<< "$last_keys" || true)
  fi
  last_key=$(tail -1 <<< "$last_keys")
  [[ -n $last_key ]] || ((${#jobs[@]} == 0)) || echo "警告: --job のジョブに起動の記録がないので、既定の設定で起動する（--type などで指定する）" >&2
  if [[ -n $last_key ]]; then
    last=$(aws s3 cp "s3://$XLC_BUCKET/$last_key" -)
    # 最後の起動が auto の選択なら、起動し直しも auto にする（同じ機種と AZ に固定すると、中断された空きの池を試し直すだけになる）
    if [[ -z $type && $market != on-demand ]] && grep -Fxq SELECTION=auto <<< "$last"; then type=auto; fi
    type=${type:-$(sed -n 's/^TYPE=//p' <<< "$last")}
    max_hours=${max_hours:-$(sed -n 's/^MAX_HOURS=//p' <<< "$last")}
    sync_minutes=${sync_minutes:-$(sed -n 's/^SYNC_MINUTES=//p' <<< "$last")}
    # auto は AZ も選ぶので、記録の AZ（前に選んだもの）で固定しない
    [[ $type == auto ]] || az=${az:-$(sed -n 's/^AZ=//p' <<< "$last")}
    market=${market:-$(sed -n 's/^MARKET=//p' <<< "$last")}
    [[ $market != on-demand ]] || echo "警告: 最後の起動はオンデマンド（約 2.5 倍）。スポットに戻すなら --spot" >&2
    echo "最後の起動（${last_key##*/}）の設定を引き継ぐ（指定がなければ）"
  fi
  todo=()
  restart_flags=()
  stuck=0
  for input in "${RUN_INPUTS[@]}"; do
    job=$(job_name "$input")
    # --job があれば、それ以外のジョブは対象にしない（状態も示さない）
    ((${#jobs[@]} == 0)) || grep -Fxq -- "$job" <<< "$(printf '%s\n' "${jobs[@]}")" || continue
    code=$(job_done "$name" "$job")
    retry_code=
    if [[ $code == 0 ]]; then
      echo "$job: 完了"
      continue
    fi
    if [[ -n $code ]]; then
      # 失敗したジョブは飛ばし、ほかのジョブは起動する（影響を 1 本に閉じる）。終わりの状態は 1
      if ((retry)) && grep -Fxq -- "$job" <<< "$(printf '%s\n' "${jobs[@]}")"; then
        retry_code=$code
      else
        echo "$job: 失敗しています（done: $code）。$run_prefix/$job/log と boot-*.log で原因を確かめる。直ったなら" \
          "--job $job --retry-failed で done を消して起動し直す" >&2
        stuck=1
        continue
      fi
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
    if [[ -n $retry_code ]]; then
      echo "$job: done=$retry_code を消して起動する（--retry-failed${dry:+、--dry-run なので消さない}。原因は $run_prefix/$job/log と boot-*.log で確かめたもの）"
      clear_done+=("$job")
    fi
    todo+=("$input")
  done
  ((${#todo[@]})) || { echo "起動するジョブはありません（aws_status.sh、aws_fetch.sh $name）"; exit "$stuck"; }
fi

# 既定（ADR-0155）。起動し直しでは、指定も記録もないときだけ
type=${type:-c7a.4xlarge}
max_hours=${max_hours:-14}
sync_minutes=${sync_minutes:-5}
market=${market:-spot}
[[ $market == spot || $market == on-demand ]] || { echo "起動の記録の MARKET が spot でも on-demand でもない: $market" >&2; exit 1; }
[[ $az != any ]] || az=
# 起動を試す候補を「順位 費用 タイプ AZ 価格 速度」で 1 行ずつ（cands）。auto でなければ指定の 1 つ（AZ は空なら AWS が選ぶ）
cands=()
if [[ $type == auto ]]; then
  [[ $market == spot ]] || { echo "--type auto はスポットの価格で並べるので --on-demand と併用できない（最後の起動の記録がオンデマンドなら --spot を付ける）" >&2; exit 2; }
  types_file=${XLC_AWS_TYPES:-$here/aws_types.toml}
  [[ -z $az ]] || az_subnet_id "$az"
  # 先に変数へ取る（< <(...) の中の失敗は止まらず、設定や認証の誤りが「空きがない」に化ける）
  cand_text=$(auto_candidates "$types_file" "$az" | awk '{print NR, $0}')
  [[ -n $cand_text ]] || { echo "候補のスポット価格がありません（$types_file の機種が ${az:-どの AZ}にない）。空きがない" >&2; exit 1; }
  # スポットの vCPU の残り（クォータ − この account の走行中）に収まらない機種を外す。起動し直しで欠けたジョブが複数あるときは、
  # 1 本が使い切って他が入らなくならないよう、残りをジョブの数で割った値を 1 本の上限にする（--job 1 本の指定ならしない）
  if remaining=$(spot_vcpu_remaining); then
    read -r remaining quota used <<< "$remaining"
    share=$remaining
    ((${#inputs[@]} || ${#jobs[@]} == 1 || ${#todo[@]} < 2)) || share=$((remaining / ${#todo[@]}))
    vcpus=$(auto_vcpus_read "$types_file")
    dropped=$(awk -v lim="$share" -v v="$vcpus" 'BEGIN { n = split(v, l, "\n"); for (i = 1; i <= n; i++) { split(l[i], f, " "); vc[f[1]] = f[2] } }
      vc[$3] > lim && !d[$3]++ { printf "%s%s(%s)", s, $3, vc[$3]; s = " " }' <<< "$cand_text")
    cand_text=$(awk -v lim="$share" -v v="$vcpus" 'BEGIN { n = split(v, l, "\n"); for (i = 1; i <= n; i++) { split(l[i], f, " "); vc[f[1]] = f[2] } }
      vc[$3] <= lim' <<< "$cand_text" | awk '{ $1 = NR; print }')
    echo "== スポットの vCPU: クォータ $quota、走行中 $used、残り $remaining、このジョブの上限 $share。外した候補: ${dropped:-なし}"
    [[ -n $cand_text ]] || { echo "上限 $share vCPU に収まる候補がない（クォータの増量か、走行中の VM の終了を待つ。起動し直しで再試行: aws_run.sh --type auto $name）" >&2; exit 1; }
  else
    echo "警告: スポットの vCPU のクォータか走行中の vCPU を取れないので、候補を絞らない（権限: servicequotas:GetServiceQuota、ec2:DescribeInstances）" >&2
  fi
  mapfile -t cands <<< "$cand_text"
  echo "== 候補（1 本の費用 = 価格 ÷ 速度の昇順。速度は $types_file の見積もりで測っていない）"
  printf '%-4s %12s  %-14s %-12s %10s %6s\n' 順位 '費用$/時間' タイプ AZ '価格$/時間' 速度
  for c in "${cands[@]}"; do
    read -r rank cost ctype caz cprice cspeed <<< "$c"
    printf '%-4s %12s  %-14s %-12s %10s %6s\n' "$rank" "$cost" "$ctype" "$caz" "$cprice" "$cspeed"
  done
else
  cands=("1 - $type ${az:--} - -")
fi
ami=$(ubuntu_ami)
echo "AMI: $ami、タイプ: $type、最長 $max_hours 時間、同期 $sync_minutes 分ごと、市場: $market、AZ: ${az:-指定なし}、statepoint の間隔: ${interval:-入力のまま}"

if ((!dry)); then
  # user-data を作ったこのスクリプトのチェックアウト（作業ディレクトリでなく）
  tools=$(git -C "$here" rev-parse HEAD)
  [[ -z $(git -C "$here" status --porcelain) ]] || tools="$tools-dirty"
fi

if ((${#inputs[@]} == 0)); then
  # 道具は手元の版に合わせる（S3 の src/ は最初の起動の版で、新しい引数を知らない。入力と launch.env は実行の定義なので書き換えない）
  root=$(cd "$here/../.." && pwd)
  [[ -f $root/flake.nix && -d $root/examples/vera ]] || { echo "リポジトリの中の aws_run.sh を使ってください（$root）" >&2; exit 1; }
  for f in "${tool_files[@]}"; do [[ -f $root/$f ]] || { echo "$root/$f がありません" >&2; exit 1; }; done
  if ((dry)); then
    echo "道具 ${#tool_files[@]} 個を $run_prefix/src/ に置き直す（--dry-run なので置かない）"
  else
    tools_src=$(mktemp -d)
    trap 'rm -rf "$tools_src"' EXIT
    (cd "$root" && cp --parents "${tool_files[@]}" "$tools_src/")
    aws s3 sync --only-show-errors "$tools_src" "$run_prefix/src"
  fi
fi

no_capacity=0
for i in "${!todo[@]}"; do
  input=${todo[$i]}
  job=$(job_name "$input")
  host=$name-$job
  user_data=$(aws_user_data "$bucket_region" "$XLC_BUCKET" "$name" "$job" "$input" "${restart_flags[$i]}" $((max_hours * 60)) \
    $((sync_minutes * 60)) "$library" "$closure" "$openmc" "$python" "$awscli" "$interval")
  launched=0
  for c in "${cands[@]}"; do
    read -r rank cost ctype caz cprice cspeed <<< "$c"
    [[ $caz != - ]] || caz=
    subnet_args=()
    if [[ -n $caz ]]; then
      az_subnet_id "$caz"
      subnet_args=(--subnet-id "$SUBNET_ID")
    fi
    args=(--image-id "$ami" --instance-type "$ctype" --count 1
      --iam-instance-profile Name=xlc-openmc-vm
      --instance-initiated-shutdown-behavior terminate
      --metadata-options 'HttpTokens=required,HttpEndpoint=enabled'
      --block-device-mappings 'DeviceName=/dev/sda1,Ebs={VolumeSize=40,VolumeType=gp3,DeleteOnTermination=true}'
      --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$host},{Key=Project,Value=xlc},{Key=XlcRun,Value=$name}]"
      --user-data "$user_data"
      "${subnet_args[@]}"
      --query 'Instances[0].InstanceId' --output text)
    if [[ $market == spot ]]; then args+=(--instance-market-options 'MarketType=spot,SpotOptions={SpotInstanceType=one-time}'); fi
    if ((dry)); then
      out=$(aws ec2 run-instances "${args[@]}" --dry-run 2>&1 || true)
      [[ $out == *DryRunOperation* ]] || { echo "$out" >&2; exit 1; }
      echo "$host: 起動の引数は有効（$ctype ${caz:-AZ 指定なし}、再開 ${restart_flags[$i]}、user-data $(wc -c <<< "$user_data") バイト）"
      launched=1
      break
    fi
    if out=$(aws ec2 run-instances "${args[@]}" 2>&1); then
      id=$(grep -Eo '^i-[0-9a-f]+$' <<< "$out" | head -1 || true)
      [[ -n $id ]] || { echo "run-instances の出力にインスタンス ID がない: $out" >&2; exit 1; }
      echo "$host: $id を起動（$ctype ${caz:-AZ 指定なし}、再開 ${restart_flags[$i]}）"
      auto_record=
      if [[ $type == auto ]]; then
        auto_record=$(printf 'SELECTION=auto\nSELECTION_RANK=%s/%s\nSPOT_PRICE=%s\nSPEED=%s\nCOST_PER_SPEED=%s\n' "$rank" "${#cands[@]}" "$cprice" "$cspeed" "$cost")
      fi
      # 起動ごとの記録（user-data を作る道具と、起動し直しで src/ に置き直す道具は TOOLS_COMMIT のもの。AMI、機種もこの起動のもの）。起動できたものだけ書く。
      # auto は選んだ機種と AZ、その価格と順位を残す（ADR-0171）
      printf 'TOOLS_COMMIT=%s\nAMI=%s\nTYPE=%s\nMAX_HOURS=%s\nSYNC_MINUTES=%s\nJOB=%s\nRESTART=%s\nINSTANCE=%s\nAZ=%s\nMARKET=%s\nSTATEPOINT_INTERVAL=%s\n%s' "$tools" "$ami" \
        "$ctype" "$max_hours" "$sync_minutes" "$job" "${restart_flags[$i]}" "$id" "$caz" "$market" "$interval" "${auto_record:+$auto_record$'\n'}" |
        aws s3 cp --only-show-errors - "$run_prefix/launches/$(date -u +%Y%m%dT%H%M%S.%NZ)-$job.env"
      # 起動の記録を書いてから消す（起動に失敗したときに、失敗の印を失わない）
      if grep -Fxq -- "$job" <<< "$(printf '%s\n' "${clear_done[@]:-}")"; then
        aws s3 rm --only-show-errors "$run_prefix/$job/done"
      fi
      launched=1
      break
    fi
    # auto の空き、クォータ、価格の失敗は次の候補へ。それ以外（と auto でない起動）の失敗は止まる
    if [[ $type == auto ]] && grep -Eq 'InsufficientInstanceCapacity|MaxSpotInstanceCountExceeded|VcpuLimitExceeded|SpotMaxPriceTooLow' <<< "$out"; then
      echo "$host: $ctype $caz は起動できない（$(grep -Eo 'InsufficientInstanceCapacity|MaxSpotInstanceCountExceeded|VcpuLimitExceeded|SpotMaxPriceTooLow' <<< "$out" | head -1)）。次の候補へ" >&2
      continue
    fi
    echo "$out" >&2
    exit 1
  done
  if ((!launched)); then
    echo "$host: 候補 ${#cands[@]} 個のどれにも空きがない（起動し直しで再試行: aws_run.sh --type auto --job $job $name）" >&2
    no_capacity=1
  fi
done
((dry)) || echo "進み具合: examples/vera/aws_status.sh、終わったら examples/vera/aws_fetch.sh $name"
((!no_capacity)) || exit 1
exit "${stuck:-0}"
