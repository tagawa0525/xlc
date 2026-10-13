#!/usr/bin/env bash
# aws_run.sh で流した実行の結果を S3 から取り出し、検証し、S3 から消す（ADR-0155、docs/environment.md の「自己完結な実行」）。
# リポジトリの直下から実行する。
#
# 使い方: examples/vera/aws_fetch.sh [--dry-run] [--verify] [--delete] [--dest <置き場>] <name>
#   取り出す: S3 の runs/<name>/ の src/（VM が使ったスクリプトと入力）と、完了したジョブの結果（最後の世代の statepoint、
#            summary.h5、k.json など、log、boot.log）を <置き場>/ に写し、ジョブごとの SHA256SUMS で検証する
#   --verify  さらに、手元（flake の .#openmc-python）で取り出した statepoint から extract し、VM の結果とバイト一致を確かめる
#   --delete  すべてのジョブが完了し、--verify が通ったら、S3 の runs/<name>/ を消す（--verify が要る）
#   --dest    置き場（既定 references/openmc/aws/<name>）
#   --dry-run  何もコピーせず、ジョブの状態（done）だけ出す
# 完了していないジョブ（done がないか 0 でない）があれば、取り出せる分を取り出したうえで 1 で終わる。
set -euo pipefail
# shellcheck source=aws_common.sh
source "$(dirname "$0")/aws_common.sh"

usage() { sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; }
dry=0
verify=0
delete=0
dest=""
name=""
while (($#)); do
  case $1 in
    --help | -h) usage; exit 0 ;;
    --dry-run) dry=1 ;;
    --verify) verify=1 ;;
    --delete) delete=1 ;;
    --dest) dest=$2; shift ;;
    -*) echo "知らない引数 $1" >&2; exit 2 ;;
    *) [[ -z $name ]] || { echo "name は 1 つ" >&2; exit 2; }; name=$1 ;;
  esac
  shift
done
[[ -n $name ]] || { usage; exit 2; }
((!delete || verify)) || { echo "--delete には --verify が要ります（検証してから消す）" >&2; exit 2; }
[[ -f flake.nix && -d examples/vera ]] || { echo "リポジトリの直下から実行してください" >&2; exit 1; }
dest=${dest:-references/openmc/aws/$name}
run_prefix=s3://$XLC_BUCKET/runs/$name

env=$(aws s3 cp "$run_prefix/launch.env" -) || { echo "$run_prefix/launch.env がありません" >&2; exit 1; }
read -r -a inputs <<< "$(sed -n 's/^INPUTS=//p' <<< "$env" | tr -d '"')"
echo "== $name（$(sed -n 's/^COMMIT=//p' <<< "$env")）"

done_inputs=()
status=0
for input in "${inputs[@]}"; do
  job=$(basename "$input" .toml)
  code=$(aws s3 cp "$run_prefix/$job/done" - 2> /dev/null || true)
  echo "$job: done ${code:-なし}"
  if [[ $code == 0 ]]; then done_inputs+=("$input"); else status=1; fi
done
((!dry)) || exit "$status"

aws s3 sync --only-show-errors "$run_prefix/src" "$dest/src"
for input in "${done_inputs[@]}"; do
  job=$(basename "$input" .toml)
  aws s3 sync --only-show-errors "$run_prefix/$job" "$dest/$job"
  (cd "$dest/$job" && sha256sum --quiet -c SHA256SUMS) || { echo "$job: SHA256SUMS が合いません" >&2; exit 1; }
  echo "$job: 取り出して SHA256SUMS が一致（$dest/$job）"
done

if ((verify)); then
  python=$(nix build .#openmc-python --no-link --print-out-paths)/bin/python3
  for input in "${done_inputs[@]}"; do
    job=$(basename "$input" .toml)
    case $job in
      openmc-3?) script=examples/vera/openmc_model_3d.py ;;
      openmc-5?-2D) script=examples/vera/openmc_core_2d.py ;;
      *) echo "$job: extract のツールが分かりません" >&2; exit 1 ;;
    esac
    check=$(mktemp -d)
    "$python" "$dest/src/$script" extract "$dest/src/$input" "$dest/$job/run" "$check" > /dev/null
    n=0
    while IFS= read -r -d '' f; do
      cmp -s "$check/$f" "$dest/$job/$f" || { echo "$job: $f が VM の結果と違います" >&2; exit 1; }
      n=$((n + 1))
    done < <(cd "$check" && find . -type f -printf '%P\0')
    rm -rf "$check"
    ((n)) || { echo "$job: extract が何も書きませんでした" >&2; exit 1; }
    echo "$job: 手元の extract の $n 個のファイルが VM の結果とバイト一致"
  done
fi

if ((delete)); then
  ((status == 0)) || { echo "完了していないジョブがあるので、S3 の $run_prefix/ は消しません" >&2; exit 1; }
  aws s3 rm --only-show-errors --recursive "$run_prefix/"
  echo "S3 の $run_prefix/ を消しました"
fi
exit "$status"
