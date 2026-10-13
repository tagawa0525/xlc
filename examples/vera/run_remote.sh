#!/usr/bin/env bash
# VM 上で、Problem 3 と 5-2D の OpenMC の参照解を並列に流す（Slurm は使わない。docs/environment.md の「外部の計算機での参照解」）。
# remote_send.sh が送った作業ディレクトリの直下から走らせる。
#
# 使い方: run_remote.sh [--restart] [--cores <物理コア数>] [--threads <1 本あたりのスレッド数>] [--statepoint-interval <世代>] <入力>...
#   <入力>   examples/vera/p3/openmc-3A.toml など（作業ディレクトリからの相対）。1 つにつき 1 本のジョブを並列に流す
#   --cores  VM の物理コア数。1 本あたりのスレッド数 = 物理コア数 / ジョブの数（切り捨て、最小 1）。省略すると lscpu の物理コア数。
#            c7a は SMT がなく vCPU が物理コア。SMT のある機種では lscpu の物理コア数が vCPU の半分。--threads があればそちらを使う
#   --threads  1 本あたりのスレッド数（入力の `threads` を上書きする。openmc_resume.py の `run --threads`）
#   --statepoint-interval  入力の `statepoint-interval` を上書きする（openmc_resume.py の `run --statepoint-interval`。入力のファイルは
#            変えない。全ジョブに同じ値。ADR-0168）
#   --restart  作業ディレクトリ（results/<入力の名前>/run）の最新の statepoint から再開する（スポットの中断のあと）
#
# 各ジョブは setsid と nohup で切り離し、接続が切れても続く。ログは results/<名前>/log、完了すると results/<名前>/done
# （成功なら 0、失敗なら終了コード）を書く。run が終わると続けて extract も走らせ、k.json などは results/<名前>/ に書く。
# 進み具合: tail -n 3 results/*/log、ls results/*/done。戻すのは remote_fetch.sh。
set -euo pipefail

restart=()
interval=()
cores=""
threads=""
inputs=()
while (($#)); do
  case $1 in
    --help | -h) sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    --restart) restart=(--restart) ;;
    --statepoint-interval) interval=(--statepoint-interval "${2:?--statepoint-interval に値が要ります}"); shift ;;
    --cores) cores=$2; shift ;;
    --threads) threads=$2; shift ;;
    -*) echo "知らない引数 $1" >&2; exit 2 ;;
    *) inputs+=("$1") ;;
  esac
  shift
done
((${#inputs[@]})) || { echo "入力を 1 つ以上渡してください（--help）" >&2; exit 2; }

source ./tools.env
export PATH="$OPENMC_BIN:$PATH"

if [[ -z $threads ]]; then
  if [[ -z $cores ]]; then
    cores=$(lscpu -p=core,socket | grep -v '^#' | sort -u | wc -l)
  fi
  threads=$((cores / ${#inputs[@]}))
  ((threads >= 1)) || threads=1
fi
echo "ジョブ ${#inputs[@]} 本、1 本あたり $threads スレッド（物理コア ${cores:-指定なし}）"

job() { # <script> <入力> <出力ディレクトリ> <スレッド数> [--restart]
  local script=$1 input=$2 out=$3 threads=$4 code
  shift 4
  "$PYTHON" "$script" run "$@" --threads "$threads" "$input" "$out/run" && \
    "$PYTHON" "$script" extract "$input" "$out/run" "$out"
  code=$?
  echo "$code" > "$out/done"
}
export -f job
export PYTHON

for input in "${inputs[@]}"; do
  [[ -f $input ]] || { echo "$input がありません" >&2; exit 1; }
  name=$(basename "$input" .toml)
  case $name in
    openmc-3?) script=examples/vera/openmc_model_3d.py ;;
    openmc-5?-2D) script=examples/vera/openmc_core_2d.py ;;
    *) echo "$input は Problem 3 か 5-2D の入力ではありません" >&2; exit 1 ;;
  esac
  out=results/$name
  mkdir -p "$out"
  rm -f "$out/done"
  # 再開のときは、ログに追記する。statepoint は $out/run にある
  setsid nohup bash -c 'job "$@"' _ "$script" "$input" "$out" "$threads" "${restart[@]}" "${interval[@]}" >>"$out/log" 2>&1 < /dev/null &
  echo "$name: $out/log"
done
