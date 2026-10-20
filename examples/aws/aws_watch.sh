#!/usr/bin/env bash
# 自己完結な実行（aws_run.sh。ADR-0155）のジョブを一定の間隔で見て、欠けたジョブを 1 回にまとめて、スポットで起動し直す
# （docs/environment.md の「自己完結な実行」）。S3 の done の中身が 0 のジョブは完了、VM が起動中・実行中でないジョブは欠けたジョブ。
# 起動し直しは aws_run.sh --job J1 --job J2 ... --spot --type auto [--statepoint-interval N] [--retry-failed] <name>
# （欠けたジョブの全部を 1 回の呼び出しに渡す。aws_run.sh が 1 本の vCPU の上限を残り ÷ 起動する数で割る。ADR-0171）
# （N は定義の progress が statepoint のときだけで、1 回の呼び出しに 1 つ。欠けたジョブがすべて openmc-3? なら 50、ほかが混じれば 25。
# 欠けたジョブに done が 0 以外のものがあれば --retry-failed。起動の上限を超える分は次の回に回す）。
# 実行名ごとに flock のロックを取る（同じ <name> の監視は 1 つだけ）。取れなければ「すでに走っている」と出して 1 で終わる。
# 終了コード: 0 全部 done、1 すでに走っている（--status では走っていない）、2 起動し直しの上限、3 想定外の失敗、
#            4 認証切れ（aws login をやり直す。VM は動き続ける）
#
# 使い方: examples/aws/aws_watch.sh [--interval 秒] [--limit 回] [--status] <name>
#   --interval  見る間隔（既定 300 秒）
#   --limit     起動し直しの回数の上限（既定 60。実際に起動した台数を数える）
#   --status    監視が走っているかをロックで確かめる（走っていれば 0、走っていなければ 1）
set -euo pipefail
here=$(dirname "$0")
# shellcheck source=aws_common.sh
source "$here/aws_common.sh"

usage() { sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'; }
interval=300 limit=60 name='' status=0
while (($#)); do
  case $1 in
    --help | -h) usage; exit 0 ;;
    --interval) interval=${2:?--interval に値が要ります}; shift ;;
    --limit) limit=${2:?--limit に値が要ります}; shift ;;
    --status) status=1 ;;
    -*) echo "知らない引数 $1" >&2; exit 2 ;;
    *) [[ -z $name ]] || { echo "name は 1 つだけ" >&2; exit 2; }; name=$1 ;;
  esac
  shift
done
[[ -n $name ]] || { usage; exit 2; }
[[ $interval =~ ^[1-9][0-9]*$ && $limit =~ ^[0-9]+$ ]] || { echo "--interval は正の整数、--limit は整数" >&2; exit 2; }

# 走っているかを pgrep -f で調べない。コマンド文字列は別の待ちループや自分のシェルにも当たる（2026-10-14 に誤認で監視が 2 時間止まり、
# pkill は自分のシェルを殺して exit 144 になった）。ロックはプロセスが終われば OS が放す
# ロックの場所は環境（XDG_RUNTIME_DIR など）で変えない。シェルによって別のファイルを取ると、二重起動を防げない
[[ $name =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || { echo "name は英数字と . _ - だけ（先頭は英数字）" >&2; exit 2; }
lock_dir=$HOME/.cache/xlc
mkdir -p "$lock_dir"
exec 9>> "$lock_dir/aws-watch-$name.lock"
if ((status)); then
  if flock -n 9; then echo "$name の監視は走っていない"; exit 1; fi
  echo "$name の監視は走っている"; exit 0
fi
# 2 秒待つ: --status の確認が一瞬ロックを取る間に、起動した監視が取り損ねない。本当の二重起動は待っても取れない
flock -w 2 9 || { echo "$name の監視はすでに走っている（確認は aws_watch.sh --status $name）" >&2; exit 1; }

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
  st=$("$here/aws_status.sh" 2>&1 9>&-) || rc=$?
  if ((rc)); then check_expired "$st"; echo "ATTENTION: aws_status.sh failed (rc=$rc)"; echo "$st"; exit 3; fi
  done_n=0 missing=() failed=()
  for j in "${jobs[@]}"; do
    rc=0
    code=$(job_done "$name" "$j" 2>&1) || rc=$?
    if ((rc)); then check_expired "$code"; echo "ATTENTION: reading done failed (rc=$rc)"; echo "$code"; exit 3; fi
    code=${code//[[:space:]]/}
    if [[ $code == 0 ]]; then
      done_n=$((done_n + 1))
    elif grep -qE "^$name-$j [^ ]+ [^ ]+ (pending|running) " <<< "$st"; then
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
    # 1 回で上限を超えないよう、残りの回数までのジョブを起動する
    missing=("${missing[@]:0:limit - relaunches}")
    args=()
    for j in "${missing[@]}"; do args+=(--job "$j"); done
    args+=(--spot --type auto)
    # statepoint の実行だけ、失う分を縮める間隔を渡す（steps は区間の境界で再開するので、間隔がない）。
    # 間隔は 1 回の呼び出しに 1 つしか渡せないので、openmc-3? だけが欠けたときだけ 50、ほかが混じれば 25
    if [[ $DEF_PROGRESS == statepoint ]]; then
      interval_n=50
      for j in "${missing[@]}"; do [[ $j == openmc-3? ]] || interval_n=25; done
      args+=(--statepoint-interval "$interval_n")
    fi
    # --retry-failed は --job のジョブのうち done が 0 でないものにだけ効く
    for j in "${missing[@]}"; do
      for f in "${failed[@]:-}"; do [[ $f != "$j" ]] || { args+=(--retry-failed); break 2; }; done
    done
    echo "relaunching ${missing[*]}"
    rc=0
    # ロックの fd 9 は子に渡さない（監視が止まったあとも走る aws_run.sh がロックを持ち続けない）
    out=$("$here/aws_run.sh" "${args[@]}" "$name" 2>&1 9>&-) || rc=$?
    tail -"$((${#missing[@]} + 6))" <<< "$out"
    # 起動した台数を数える（aws_run.sh の「<host>: i-... を起動」の行。文言を変えたら、ここと test_aws_watch.sh の偽も変える）。
    # 1 回の呼び出しで一部だけ起動することもある。起動するジョブがなかった回（VM が終了中など）は 0
    launched=$(grep -c '[^ ]*: i-[0-9a-f]* を起動' <<< "$out" || true)
    relaunches=$((relaunches + launched))
    if ((rc == 0)); then
      ((launched > 0)) || echo "nothing launched for ${missing[*]}"
    else
      check_expired "$out"
      if ((rc == 1)) && grep -qE '空きがない|候補がない' <<< "$out"; then
        echo "no capacity for ${missing[*]}; retry next round"
      else
        echo "ATTENTION: aws_run.sh failed (rc=$rc)"; echo "$st"; exit 3
      fi
    fi
  fi
  sleep "$interval" 9>&-
done
