#!/usr/bin/env bash
# VM 上で、OpenMC の燃焼の参照解を 1 つ流す（AWS の実行の定義 aws.toml の入口。ADR-0174）。作業ディレクトリ（送ったツールと入力がある
# ところ）の直下から実行する。
#
# 使い方: run_remote.sh [--restart] [--threads <スレッド数>] [--extract-script] <入力>
#   <入力>   examples/depletion/openmc-1C-rings1.toml など（作業ディレクトリからの相対）
#   --restart  results/<入力の名前>/run の最後の完了した区間から再開する（スポットの中断のあと）
#   --threads  スレッド数（入力の `threads` を上書きする）。省略すると lscpu の物理コア数
#   --extract-script  流さずに、<入力>の run と extract を行うスクリプト（リポジトリの直下からの相対）を出す
#            （aws_fetch.sh --verify が手元で extract をやり直すために問い合わせる）
#
# ジョブは setsid と nohup で切り離し、接続が切れても続く。ログは results/<名前>/log、完了すると results/<名前>/done
# （成功なら 0、失敗なら終了コード）を書く。run が終わると続けて extract も実行し、<名前>.msgpack は results/<名前>/ に書く。
set -euo pipefail

script=examples/depletion/openmc_depletion_resume.py
restart=()
threads=""
query=0
inputs=()
while (($#)); do
  case $1 in
    --help | -h) sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    --restart) restart=(--restart) ;;
    --threads) threads=${2:?--threads に値が要ります}; shift ;;
    --extract-script) query=1 ;;
    -*) echo "知らない引数 $1" >&2; exit 2 ;;
    *) inputs+=("$1") ;;
  esac
  shift
done
((${#inputs[@]} == 1)) || { echo "入力を 1 つ渡してください（--help）" >&2; exit 2; }
input=${inputs[0]}
if ((query)); then echo "$script"; exit 0; fi

source ./tools.env
export PATH="$OPENMC_BIN:$PATH"
[[ -f $input ]] || { echo "$input がありません" >&2; exit 1; }
if [[ -z $threads ]]; then
  threads=$(lscpu -p=core,socket | grep -v '^#' | sort -u | wc -l)
fi
echo "ジョブ 1 本、$threads スレッド"

name=$(basename "$input" .toml)
out=results/$name
mkdir -p "$out"
rm -f "$out/done"
# OMP_NUM_THREADS は OpenMC の読み込みの前に要る（Python の中で設定しても、読み込み済みのライブラリには効かない）
job() { # <入力> <出力ディレクトリ> <スレッド数> [--restart]
  local input=$1 out=$2 threads=$3 code
  shift 3
  OMP_NUM_THREADS=$threads "$PYTHON" "$script" run "$@" --threads "$threads" "$input" "$out/run" && \
    "$PYTHON" "$script" extract "$input" "$out/run" "$out"
  code=$?
  echo "$code" > "$out/done"
}
export -f job
export PYTHON script
# 再開のときは、ログに追記する。区間は $out/run にある
setsid nohup bash -c 'job "$@"' _ "$input" "$out" "$threads" "${restart[@]}" >>"$out/log" 2>&1 < /dev/null &
echo "$name: $out/log"
