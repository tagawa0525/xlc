#!/usr/bin/env bash
# 自己完結な実行（aws_run.sh。ADR-0155）のジョブを一定の間隔で見て、欠けたジョブを 1 回に 1 つ、スポットで起動し直す
# （docs/environment.md の「自己完結な実行」）。S3 の done の中身が 0 のジョブは完了、VM が走っていないジョブは欠けたジョブ。
# 起動し直しは aws_run.sh --job J --spot --type auto --statepoint-interval N [--retry-failed] <name>
# （N はジョブ名が openmc-3? なら 50、それ以外は 25。done が 0 以外なら --retry-failed）。
# 終了コード: 0 全部 done、2 起動し直しの上限、3 想定外の失敗、4 認証切れ（aws login をやり直す。VM は動き続ける）
#
# 使い方: examples/vera/aws_watch.sh [--interval 秒] [--limit 回] <name>
#   --interval  見る間隔（既定 300 秒）
#   --limit     起動し直しの回数の上限（既定 60。実際に起動した回だけ数える）
set -euo pipefail
here=$(dirname "$0")
# shellcheck source=aws_common.sh
source "$here/aws_common.sh"

usage() { sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'; }
interval=300 limit=60 name=
while (($#)); do
  case $1 in
    --help | -h) usage; exit 0 ;;
    --interval) interval=${2:?--interval に値が要ります}; shift ;;
    --limit) limit=${2:?--limit に値が要ります}; shift ;;
    -*) echo "知らない引数 $1" >&2; exit 2 ;;
    *) [[ -z $name ]] || { echo "name は 1 つだけ" >&2; exit 2; }; name=$1 ;;
  esac
  shift
done
[[ -n $name ]] || { usage; exit 2; }
[[ $interval =~ ^[1-9][0-9]*$ && $limit =~ ^[0-9]+$ ]] || { echo "--interval は正の整数、--limit は整数" >&2; exit 2; }

expired='ExpiredToken|RequestExpired|session has expired'
# 失敗の出力を $1 に見て、認証切れなら 4 で終わる
check_expired() {
  if grep -qE "$expired" <<< "$1"; then echo "ATTENTION: aws login expired"; exit 4; fi
}

read_out=$(read_launch_env "$name" 2>&1) || { check_expired "$read_out"; echo "$read_out" >&2; exit 3; }
read_launch_env "$name"
jobs=()
for input in "${RUN_INPUTS[@]}"; do jobs+=("$(job_name "$input")"); done
total=${#jobs[@]}

relaunches=0
while :; do
  rc=0
  st=$("$here/aws_status.sh" 2>&1) || rc=$?
  if ((rc)); then check_expired "$st"; echo "ATTENTION: aws_status.sh failed (rc=$rc)"; echo "$st"; exit 3; fi
  done_n=0 missing=() failed=()
  for j in "${jobs[@]}"; do
    rc=0
    code=$(job_done "$name" "$j" 2>&1) || rc=$?
    if ((rc)); then check_expired "$code"; echo "ATTENTION: reading done failed (rc=$rc)"; echo "$code"; exit 3; fi
    code=${code//[[:space:]]/}
    if [[ $code == 0 ]]; then
      done_n=$((done_n + 1))
    elif grep -q "$name-$j .* running " <<< "$st"; then
      :
    else
      missing+=("$j")
      [[ -z $code ]] || failed+=("$j")
    fi
  done
  echo "$(date '+%m-%d %H:%M') done=$done_n/$total missing=[${missing[*]:-}] failed=[${failed[*]:-}] relaunches=$relaunches"
  if ((done_n == total)); then echo "ALL DONE"; exit 0; fi
  if ((${#missing[@]})); then
    if ((relaunches >= limit)); then
      echo "ATTENTION: relaunch limit reached for ${missing[*]}"; echo "$st"; exit 2
    fi
    job=${missing[0]}
    interval_n=25
    [[ $job != openmc-3? ]] || interval_n=50
    args=(--job "$job" --spot --type auto --statepoint-interval "$interval_n")
    for f in "${failed[@]:-}"; do [[ $f != "$job" ]] || args+=(--retry-failed); done
    echo "relaunching $job"
    rc=0
    out=$("$here/aws_run.sh" "${args[@]}" "$name" 2>&1) || rc=$?
    tail -6 <<< "$out"
    if ((rc == 0)); then
      relaunches=$((relaunches + 1))
    else
      check_expired "$out"
      if ((rc == 1)) && grep -qE '空きがない|候補がない' <<< "$out"; then
        echo "no capacity for $job; retry next round"
      else
        echo "ATTENTION: aws_run.sh failed (rc=$rc)"; echo "$st"; exit 3
      fi
    fi
  fi
  sleep "$interval"
done
