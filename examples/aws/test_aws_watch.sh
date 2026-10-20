#!/usr/bin/env bash
# aws_watch.sh の検査: 同じ実行の二重起動は 1 で終わる、--status がロックで走っているかを答える、欠けたジョブは 1 回の
# aws_run.sh にまとめて渡る。偽の aws_common.sh、aws_status.sh、aws_run.sh を aws_watch.sh の隣に置いて実行する。
# flake の check（aws-watch）が実行する。
#
# 使い方: examples/aws/test_aws_watch.sh [<aws_watch.sh のパス>]
set -euo pipefail
watch=$(realpath "${1:-$(dirname "$0")/aws_watch.sh}")

failures=0
fail() {
  echo "FAIL: $*" >&2
  failures=$((failures + 1))
}
expect() { # <名前> <期待> <実際>
  [[ $2 == "$3" ]] || fail "$1: 期待 '$2'、実際 '$3'"
}

tmp=$(mktemp -d)
trap 'kill $(jobs -p) 2> /dev/null || true; rm -rf "$tmp"' EXIT
mkdir "$tmp/bin" "$tmp/home"
cp "$watch" "$tmp/bin/aws_watch.sh"
export HOME=$tmp/home STATE=$tmp/state
mkdir "$STATE"

# 偽の aws_common.sh: 入力は 4 本。done は $STATE/done-<job> の中身（なければ空）
cat > "$tmp/bin/aws_common.sh" <<'COMMON'
DEF_PROGRESS=statepoint
read_launch_env() { RUN_INPUTS=(a/openmc-3A.toml a/openmc-3B.toml a/openmc-5B.toml a/openmc-5C.toml); }
job_name() { basename "$1" .toml; }
job_done() { if [[ -f $STATE/done-$2 ]]; then cat "$STATE/done-$2"; fi; }
COMMON
# 偽の aws_status.sh: t3 のジョブは実行中（欠けないので、監視は待ち続ける）、ほかは VM なし
cat > "$tmp/bin/aws_status.sh" <<'STATUS'
for j in openmc-3A openmc-3B openmc-5B openmc-5C; do echo "t3-$j i-0 c8a.4xlarge running 0"; done
STATUS
# 偽の aws_run.sh: 引数を記録し、--job のジョブを done にして（NO_DONE があれば done にしない）、ジョブごとに起動の行を出す
cat > "$tmp/bin/aws_run.sh" <<'RUN'
echo "$*" >> "$STATE/calls"
args=("$@")
for ((i = 0; i < ${#args[@]}; i++)); do
  if [[ ${args[i]} == --job ]]; then
    [[ -n ${NO_DONE:-} ]] || echo 0 > "$STATE/done-${args[i + 1]}"
    echo "${args[i + 1]}: i-0123abcd を起動（c8a.4xlarge us-east-1a、再開 なし）"
  fi
done
RUN
chmod +x "$tmp/bin/"*.sh

# 欠けた 4 本は 1 回の呼び出しにまとまり、起動した 4 台が回数に数えられる
out=$(timeout 20 bash "$tmp/bin/aws_watch.sh" --interval 1 --limit 60 t1)
expect "呼び出しは 1 回" 1 "$(wc -l < "$STATE/calls")"
expect "欠けた 4 本が 1 回で渡る" \
  "--job openmc-3A --job openmc-3B --job openmc-5B --job openmc-5C --spot --type auto --statepoint-interval 25 t1" "$(cat "$STATE/calls")"
grep -q 'relaunches=0$' <<< "$out" || fail "最初の回の relaunches が 0 でない: $out"
grep -q 'ALL DONE' <<< "$out" || fail "ALL DONE が出ない: $out"

# 起動の上限は台数で数える: 上限 4 のとき、4 台起動した次の回で止まる
rm -f "$STATE"/done-* "$STATE/calls"
rc=0
# 起動のあとも done にならない偽（2 回目の回で上限に当たる）
out=$(NO_DONE=1 timeout 20 bash "$tmp/bin/aws_watch.sh" --interval 1 --limit 4 t2) || rc=$?
expect "上限は終了コード 2" 2 "$rc"
expect "上限の前に呼び出しは 1 回" 1 "$(wc -l < "$STATE/calls")"
grep -q 'relaunches=4$' <<< "$out" || fail "起動した台数 4 が relaunches に数えられない: $out"

# 3A、3B だけが欠けたときは間隔 50、done が 1 のジョブがあれば --retry-failed
rm -f "$STATE"/done-* "$STATE/calls"
echo 0 > "$STATE/done-openmc-5B"
echo 0 > "$STATE/done-openmc-5C"
echo 1 > "$STATE/done-openmc-3A"
timeout 20 bash "$tmp/bin/aws_watch.sh" --interval 1 --limit 60 t4 > /dev/null
expect "3? だけなら間隔 50 と --retry-failed" \
  "--job openmc-3A --job openmc-3B --spot --type auto --statepoint-interval 50 --retry-failed t4" "$(cat "$STATE/calls")"

# 上限の残りが少ないときは、残りの数までのジョブだけを渡す（上限を超えて起動しない）
rm -f "$STATE"/done-* "$STATE/calls"
rc=0
out=$(NO_DONE=1 timeout 20 bash "$tmp/bin/aws_watch.sh" --interval 1 --limit 2 t5) || rc=$?
expect "残り 2 本までのジョブ" \
  "--job openmc-3A --job openmc-3B --spot --type auto --statepoint-interval 50 t5" "$(cat "$STATE/calls")"
expect "上限で終了コード 2" 2 "$rc"

# 二重起動: 走っている間は同じ実行の 2 つ目が 1 で終わり、--status は走っていると答える
rm -f "$STATE/calls"
bash "$tmp/bin/aws_watch.sh" --interval 30 --limit 0 t3 > "$tmp/t3.log" &
# 最初の回の行が出れば、ロックは取れている（--status の確認とロックの取得が重なって取り損ねるのを避ける）
for _ in {1..100}; do grep -q 'done=' "$tmp/t3.log" && break; sleep 0.1; done
grep -q 'done=' "$tmp/t3.log" || fail "t3 の監視が最初の回の行を出さない"
rc=0
out=$(timeout 20 bash "$tmp/bin/aws_watch.sh" t3 2>&1) || rc=$?
expect "二重起動は 1" 1 "$rc"
grep -q 'すでに走っている' <<< "$out" || fail "二重起動のメッセージがない: $out"
rc=0
timeout 20 bash "$tmp/bin/aws_watch.sh" --status t3 > /dev/null || rc=$?
expect "--status は走っていれば 0" 0 "$rc"
rc=0
timeout 20 bash "$tmp/bin/aws_watch.sh" --status other > /dev/null || rc=$?
expect "別の実行の --status は走っていない（1）" 1 "$rc"
kill %1
wait || true
rc=0
timeout 20 bash "$tmp/bin/aws_watch.sh" --status t3 > /dev/null || rc=$?
expect "終わればロックは放される" 1 "$rc"

((failures == 0)) || { echo "$failures 件失敗" >&2; exit 1; }
echo "aws_watch: ok"
