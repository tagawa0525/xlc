#!/usr/bin/env bash
# aws_run.sh で流した実行の結果を S3 から取り出し、検証し、S3 から消す（ADR-0155、docs/environment.md の「自己完結な実行」）。
# リポジトリの直下から実行する。
#
# 使い方: examples/vera/aws_fetch.sh [--dry-run] [--verify] [--delete] [--dest <置き場>] <name>
#   取り出す: S3 の runs/<name>/ の src/（VM が使ったスクリプトと入力）と、完了したジョブの結果（最後の世代の statepoint、
#            summary.h5、k.json など、log、boot-*.log）を <置き場>/ に写し、ジョブごとの SHA256SUMS で検証する
#   --verify  さらに、VM と同じ Python（launch.env の PYTHON。手元の store にあること）で取り出した statepoint から extract し、
#             VM の結果とバイト一致を確かめる
#   --delete  すべてのジョブが完了し、--verify が通ったら、S3 の runs/<name>/ を消す（--verify が要る）
#   --dest    置き場（既定 references/openmc/aws/<name>）
#   --dry-run  何もコピーせず、ジョブの状態（done）だけ出す
# 完了していないジョブ（done がないか 0 でない）があれば、取り出せる分を取り出したうえで 1 で終わる。
set -euo pipefail

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
# shellcheck source=aws_common.sh
source "$(dirname "$0")/aws_common.sh"
dest=${dest:-references/openmc/aws/$name}
run_prefix=s3://$XLC_BUCKET/runs/$name

read_launch_env "$name"
echo "== $name（$(launch_env_value COMMIT)）"

done_inputs=()
status=0
for input in "${RUN_INPUTS[@]}"; do
  job=$(job_name "$input")
  code=$(job_done "$name" "$job")
  echo "$job: done ${code:-なし}"
  if [[ $code == 0 ]]; then done_inputs+=("$input"); else status=1; fi
done
((!dry)) || exit "$status"

# --delete: 同じ名前の前の取り出しの残りを SHA256SUMS の検証に混ぜない（置き場は S3 のコピーにする）
aws s3 sync --only-show-errors --delete "$run_prefix/src" "$dest/src"
for input in "${done_inputs[@]}"; do
  job=$(job_name "$input")
  aws s3 sync --only-show-errors --delete "$run_prefix/$job" "$dest/$job"
  (cd "$dest/$job" && sha256sum --quiet -c SHA256SUMS) || { echo "$job: SHA256SUMS が合いません" >&2; exit 1; }
  echo "$job: 取り出して SHA256SUMS が一致（$dest/$job）"
done

if ((verify)); then
  # VM と同じ Python（launch.env の store のパス）で extract する。手元になければ、そのコミットで nix build .#openmc-python する
  python=$(launch_env_value PYTHON)/bin/python3
  [[ -x $python ]] || { echo "$python が手元にありません（git checkout $(launch_env_value COMMIT) で nix build .#openmc-python）" >&2; exit 1; }
  check=$(mktemp -d)
  trap 'rm -rf "$check"' EXIT
  for input in "${done_inputs[@]}"; do
    job=$(job_name "$input")
    script=$(job_script "$job")
    [[ -n $script ]] || { echo "$job: extract のツールが分かりません" >&2; exit 1; }
    rm -rf "${check:?}"/*
    "$python" "$dest/src/$script" extract "$dest/src/$input" "$dest/$job/run" "$check" > /dev/null
    n=0
    while IFS= read -r -d '' f; do
      cmp -s "$check/$f" "$dest/$job/$f" || { echo "$job: $f が VM の結果と違います" >&2; exit 1; }
      n=$((n + 1))
    done < <(cd "$check" && find . -type f -printf '%P\0')
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
